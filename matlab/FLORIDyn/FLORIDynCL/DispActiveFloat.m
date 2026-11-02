function [T, x_disp, y_disp, x_dot, y_dot, theta, thetad, Faero] = DispActiveFloat(T, tmpM, Sim, paramFLORIDyn, paramFLORIS, it, Ct)

%==========================================================================
% OLD FLOATER IN "FLORIDynCL/Dispfloater"
%==========================================================================

% DispActiveFloat
% -------------------------------------------------------------------------
% ActiveFloat platform dynamics using a LINEAR 2-DOF model (surge + pitch):
%
%   M2 * qddot + Beff2 * qdot + Keff2 * q = F2
%
% where:
%   q    = [surge; pitch]
%   qdot = [surge_rate; pitch_rate]
%   F2   = [Fs; Mtheta]   (generalized external loads)
%
% IMPORTANT MODELING NOTE:
% - The provided matrices (M, AddBlin, K, AddCLin) typically already include:
%     * total mass + hydrodynamic added mass  -> M
%     * linear radiation damping (+ extra lin. damping) -> AddBlin
%     * hydrostatic + gravity + (linearized) mooring stiffness -> K (+ AddCLin)
% - Therefore, we start with ONLY AERODYNAMIC EXCITATION on the RHS (F2).
%   Do not add MooringForce/HydroForce here unless you are explicitly adding
%   nonlinear residual forces not already represented by the matrices.
%
% Integration:
% - Fixed sub-stepping with semi-implicit Euler (dt_sub = dt/50) to remain
%   synchronous with the FLORIDyn macro time step.
%
% Outputs:
% - Updates T.floater.{x,vx,theta,theta_dot}
% - Updates T.posBase using surge + pitch-induced nacelle x-shift
% -------------------------------------------------------------------------





    % -------------------- Current states (vectors, one per turbine) -------
    s      = T.floater.x(:,1);          % surge [m]
    sd     = T.floater.vx(:,1);         % surge rate [m/s]
    theta  = T.floater.theta(:,1);      % pitch [rad]
    thetad = T.floater.theta_dot(:,1);  % pitch rate [rad/s]

    % Keep sway variables unchanged in this 2-DOF implementation
    y_disp = T.floater.y(:,1);
    y_dot  = T.floater.vy(:,1);

    % -------------------- Extract 2x2 matrices (surge=1, pitch=5) --------
    % DOF order assumed: [surge sway heave roll pitch yaw] = [1 2 3 4 5 6]
    idx = [1 5];

    % M2    = T.floater.MatFloat.M(idx, idx);
    % K2    = T.floater.MatFloat.K(idx, idx);
    % Cadd2 = T.floater.MatFloat.AddCLin(idx, idx);
    % B2    = T.floater.MatFloat.AddBLin(idx, idx);
    % 
    % % trying to reach out 5° of pitch instead of 7°
    % % Cadd2(2,2) = Cadd2(2,2) + 9*Cadd2(2,2);
    % % M2(2,2) = M2(2,2) + 0.55*M2(2,2);
    % 
    % Keff2 = K2 + Cadd2;
    % Beff2 = B2;
    M2 = T.floater.VolturnUS.M;
    Beff2 = T.floater.VolturnUS.B;
    Keff2 = T.floater.VolturnUS.K;

    % Added mass a frequenza infinita (Table 4, Allen et al.)
    A_inf = [9.64e6   -1.01e8;
            -1.01e8    1.16e10];
    
    % Rigidezza idrostatica WAMIT
    Chs55 = 2.182166e5 * 1025 * 9.8066;  % ≈ 2.193e9
    
    % Contributo gravitazionale
    m_tot = 20706630;    % massa sistema completo [kg]
    g     = 9.8066;
    z_cg  = -1.9;        % CG sistema sotto SWL [m]
    C_grav = -m_tot * g * z_cg;  % ≈ +3.86e8
    
    % Correzioni
    M2 = M2 + A_inf;                          % aggiungi added mass
    Keff2(2,2) = Keff2(2,2) + Chs55 + C_grav; % aggiungi idrostatica + gravitazionale

    % -------------------- Rotor area ------------------------------------
    A = pi*(T.D/2).^2;

    % -------------------- Effective wind speed --------------------------
    Ue_now = tmpM(:,3);

    % -------------------- Surge damping: hydro quadratic + aero stab ---
    %
    % 1) QUADRATIC HYDRODYNAMIC DRAG (linearized around current sd)
    %    ActiveFloat COREWIND (Mahfouz et al. 2021, Wind Energ. Sci.)
    %    Semi-submersible: 3 external columns + 1 central, D=12m, draft=26.5m
    %    F_drag = -0.5*rho*Cd*A_wet*sd*|sd|  ->  B_quad ≈ rho*Cd*A_wet*|sd|
    rho_water    = 1025;                      % seawater density [kg/m³]
    Cd_surge     = 0.7;                       % transverse Cd, COREWIND value [-]
    D_col        = 12.0;                      % external column diameter [m]
    draft        = 26.5;                      % platform draft [m]
    n_col        = 4;                         % 3 external + 1 central column
    A_wet_total  = n_col * D_col * draft;     % projected wetted area [m²]

    sd_now       = T.floater.vx(:,1);
    B_quad_surge = 0.5*rho_water * Cd_surge * A_wet_total * abs(sd_now);
    B_quad_surge = max(B_quad_surge, 1e4);    % floor when sd≈0 at startup

    % 2) AERODYNAMIC STABILIZING DAMPING (linearization of dFT/d(sd))
    %    FT = 0.5*rho_air*A*Ct*(Ue-sd)^2, linearized around sd=0:
    %    FT ≈ 0.5*rho_air*A*Ct*Ue² - (rho_air*A*Ct*Ue)*sd

    idxT = T.StartI;         % turbine indices at current time step
    yaw = deg2rad(T.States_T(idxT(:),2));

    dFdU = paramFLORIS.airDen .* A .* Ct .* Ue_now .* cos(yaw).^paramFLORIS.p_t;  % per turbine [Ns/m]
    B_aero_stab = dFdU;

    % Apply both terms to surge DOF only (1,1)
    % Beff2(1,1) = Beff2(1,1) + B_quad_surge + B_aero_stab;

    debug = true;
    if debug == true

        % ========== DEBUG DAMPING ==========================================
        fprintf('\n--- FLOATER DEBUG (it=%d) ---\n', it);
        fprintf('Beff2 after damping additions:\n');
        disp(Beff2);
        fprintf('Keff2:\n');
        disp(Keff2);
    
        [~, Deig] = eig(Keff2, M2);
        fn = sqrt(diag(Deig)) / (2*pi);
        fprintf('Natural frequencies surge/pitch [Hz]: %.4f  %.4f\n', fn(1), fn(2));
    
        fprintf('B_quad_surge   = %.4e Ns/m\n', B_quad_surge(1));
        fprintf('B_aero_stab    = %.4e Ns/m\n', B_aero_stab(1));
        fprintf('B_aero_neg     = %.4e Ns/m\n', dFdU(1));
        fprintf('B_net_surge    = %.4e Ns/m  <-- negative = unstable\n', ...
            (Beff2(1,1) + B_quad_surge(1) + B_aero_stab(1)) - dFdU(1));
    
        Fs_static        = 0.5 * paramFLORIS.airDen * A(1) * Ct(1) * Ue_now(1)^2;
        Mth_static       = T.posNac(1,3) * Fs_static;
        theta_static_deg = rad2deg(Mth_static / Keff2(2,2));
        surge_static     = Fs_static / Keff2(1,1);
        fprintf('Expected static pitch turbine 1: %.2f deg\n', theta_static_deg);
        fprintf('Expected static surge turbine 1: %.2f m\n',   surge_static);
        fprintf('------------------------------------------\n');
        % ========== END DEBUG ==============================================

    end

    % -------------------- Time stepping ---------------------------------
    dt     = Sim.TimeStep;   % FLORIDyn macro step [s]
    Nsub   = 50;
    dt_sub = dt / Nsub;


        
    % -------------------- Aerodynamic thrust (from Ct) ------------------
    % Rotor area and local effective wind speed at rotor (Ue)

    
    Ue  = tmpM(:,3);

    a = T.States_T(idxT(:),1);
    % Ct   = 4*a.*(1-a);
    rot_tilt = ones(size(a)) .* deg2rad(paramFLORIS.rotor_tilt);
    eff_tilt = rot_tilt + theta;


    % eff_tilt = eff_tilt.*0;

    if T.floater_use == true
        FT = 0.5 * paramFLORIS.airDen * A .* Ct .* Ue.^2 .* cos(yaw).^paramFLORIS.p_t;        
    else
        eff_angle = acos(cos(yaw) .* cos(eff_tilt));
        skew_angle = (0.6 * a + 1) .* eff_angle; 

        Ct_loc = 4 * a .* (cos(eff_angle) + tan(skew_angle / 2).*...
            sin(eff_angle) - a .* (1 ./ cos(skew_angle / 2)) .^ 2);
        FT = 0.5 * paramFLORIS.airDen * A .* Ct_loc .* Ue.^2;
    end

    % -------------------- External forces (build generalized load) ------
    % Compute aerodynamic force components (Faero: [Fx Fy ...] from your helper)
    
    yaw_iT = deg2rad(T.States_T(idxT(:),2));
    Faero  = AeroForce(FT, yaw_iT);
    
    % -------------------- Build generalized loads F2 = [Fs; Mtheta] ------
    % Surge generalized force:
    Fs = Faero(:,1);   

    % Pitch generalized moment:
    z_hub = T.posNac(:,3);          
    Mtheta = z_hub .* Fs;           
    

    % -------------------- Integrate platform ODE for each turbine --------
    for iT = 1:T.nT

        % State vectors for this turbine
        q    = [ s(iT);    theta(iT)  ];   % [surge; pitch]
        qdot = [ sd(iT);   thetad(iT) ];   % [surge_rate; pitch_rate]

        % Hold generalized loads constant during the macro-step
        F2 = [ Fs(iT); Mtheta(iT) ];

        
        % Fixed-step semi-implicit Euler (stable and synchronous)
        for k = 1:Nsub
            Beff2_iT = Beff2;
            Beff2_iT(1,1) = Beff2_iT(1,1)+ B_aero_stab(iT); % + B_quad_surge(iT) ;
            qddot = M2 \ (F2 - Beff2_iT*qdot - Keff2*q);
            qdot  = qdot + dt_sub * qddot;
            q     = q    + dt_sub * qdot;
        end

        % Store back
        s(iT)      = q(1);
        theta(iT)  = q(2);
        sd(iT)     = qdot(1);
        thetad(iT) = qdot(2);
    end


    % update:
    % surge motions
    T.floater.x(:,1) = s;
    T.floater.vx(:,1) = sd;
    
    % sway %motion
    T.floater.y(:,1) = y_disp;
    T.floater.vy(:,1) = y_dot;

    % Pitch motion
    T.floater.theta(:,1) = theta;
    T.floater.theta_dot(:,1) = thetad;
    
    % disp calculation taking in account pitch:
    % Update postitions:
    x_disp_pitch = z_hub.*sin(theta);
    x_dot_pitch = (z_hub .* thetad) .* cos(theta);
    % for plotting the pitch relative wind speed component
    T.floater.pitch_ws = x_dot_pitch;
    
    % Update platform base position used by FLORIDyn/FLORIS geometry and
    % wind direction
    wd    = T.States_WF(T.StartI,2);
    alpha = mod(wd - 270, 360);

    T.posBase(:,1) = T.pos0(:,1) + (T.floater.x(:,1) + x_disp_pitch) .* cosd(alpha) - T.floater.y(:,1) .* sind(alpha);
    T.posBase(:,2) = T.pos0(:,2) + (T.floater.y(:,1)) .* cosd(alpha) + (T.floater.x(:,1) + x_disp_pitch) .* sind(alpha) ;

    x_disp = T.floater.x(:,1) + x_disp_pitch;
    x_dot  = T.floater.vx(:,1) + x_dot_pitch;

    y_disp = T.floater.y(:,1) ;
    y_dot  = T.floater.vy(:,1);

    % Optional OP update (keep your previous behavior)
    T.States_OP(T.StartI,3) = z_hub .* cos(theta);

end
