function [T, x_disp, y_disp, x_dot, y_dot, theta_disp, theta_dot, Faero] = Dispfloater(T, tmpM, Sim, paramFLORIDyn, paramFLORIS, it, Ct)

   % 
   % nT = T.nT;
    
    % Prealloc outputs (vectors  nT x 1)
    % x_disp = zeros(nT,1);
    % y_disp = zeros(nT,1);
    % x_dot  = zeros(nT,1);
    % y_dot  = zeros(nT,1);

    x_disp     = T.floater.x(:,1);
    y_disp     = T.floater.y(:,1);
    x_dot      = T.floater.vx(:,1);
    y_dot      = T.floater.vy(:,1);
    theta_disp = T.floater.theta(:,1);
    theta_dot  = T.floater.theta_dot(:,1);

    
    dt   = Sim.TimeStep;
    idxT = T.StartI;
    
    % Thrust aerodinamica (usa Ct)
    A  = pi*(T.D/2).^2;
    Ue  = tmpM(:,3);
    a = T.States_T(idxT(:),1);
    % Ct   = 4*a.*(1-a);
    rot_tilt = ones(size(a)) .* deg2rad(paramFLORIS.rotor_tilt);
    yaw= T.States_T(T.StartI,2);
    eff_tilt = rot_tilt + theta_disp;
    % eff_tilt = eff_tilt.*0;
    eff_angle = acos(cos(yaw) .* cos(eff_tilt));
    skew_angle = (0.6 * a + 1) .* eff_angle; 
    if T.floater_use == true
        FT = 0.5 * paramFLORIS.airDen * A .* Ct .* cos(yaw).* Ue.^2;        
    else
        Ct = 4 * a .* (cos(eff_angle) + tan(skew_angle / 2).*...
            sin(eff_angle) - a .* (1 ./ cos(skew_angle / 2)) .^ 2);
        FT = 0.5 * paramFLORIS.airDen * A .* Ct .* Ue.^2;
    end

    % yaw relativo rispetto al vento locale
    yaw_iT = T.States_T(idxT(:),2);
    
    anchor_radius = 837.46;
    
    % Forces
    Faero = AeroForce(FT, yaw_iT, eff_tilt);
    Fhydro = HydroForce(x_dot, y_dot);
    Fmoor = MooringForce(T.floater.x(:,1), T.floater.y(:,1), anchor_radius);
    % Moments
    Maero_pitch = T.posNac(:,3).*Faero(:,1);



    % Accelerazione
    resistance_forces_x = (Fmoor(:,1) + Fhydro(:,1));
    resistance_forces_y = (Fmoor(:,2) + Fhydro(:,2));
    Forces_x = Faero(:,1) + resistance_forces_x;
    Forces_y = Faero(:,2) + resistance_forces_y;
    ax= (Forces_x) / T.floater.M_eff;
    ay= (Forces_y) / T.floater.M_eff;

    % Angular acceleration
    I = T.floater.I_pitch;
    C = T.floater.C_pitch;
    K = T.floater.K_pitch;
    
    theta_k    = theta_disp;     % (nT x 1) rad
    thetad_k   = theta_dot;      % (nT x 1) rad/s
    M_k        = Maero_pitch;    % (nT x 1) N*m
    
    den = (1 + dt*(C/I) + dt^2*(K/I));     % scalare
    

    % Integrazione (Eulero semi-implicito)
    vx_new = x_dot + dt * ax;
    vy_new = y_dot + dt * ay;
    x_new = x_disp + dt * vx_new;
    y_new = y_disp + dt * vy_new;

    % integrazione pitch

    thetad_new = ( thetad_k + dt*(M_k./I) - dt*(K/I)*theta_k ) / den;
    theta_new  = theta_k + dt*thetad_new;
    theta_deg = rad2deg(theta_new);


    % update:
    T.floater.x = x_new;
    T.floater.y = y_new;
    T.floater.vx = vx_new;
    T.floater.vy = vy_new;
    T.floater.theta = theta_new;
    T.floater.theta_dot = thetad_new;
    
    % disp calculation taking in account pitch:
    % Update postitions:
    z = T.posNac(:,3);
    x_disp_pitch = z.*sin(T.floater.theta(:));
    x_dot_pitch = (z .* thetad_new) .* cos(theta_new);
    % for plotting the pitch relative wind speed component
    T.floater.pitch_ws = x_dot_pitch;

    T.posBase(:,1) = T.pos0(:,1) + T.floater.x(:,1) + x_disp_pitch;
    T.posBase(:,2) = T.pos0(:,2) + T.floater.y(:,1);

    x_disp = T.floater.x(:,1) + x_disp_pitch;
    y_disp = T.floater.y(:,1) ;
    x_dot  = T.floater.vx(:,1) + x_dot_pitch;
    y_dot  = T.floater.vy(:,1);
    T.States_OP(T.StartI,3) = z.*cos(T.floater.theta(:));
    

    % %% for visualization of platform motion
    % if k == 1
    %     floater.xHist = zeros(Sim.NoTimeSteps, nT);
    %     floater.vHist = zeros(Sim.NoTimeSteps, nT);
    % end

    % floater.xHist(k,:) = floater.x.';  % [1 x nT]
    % floater.vHist(k,:) = floater.v.';  % [1 x nT]

    end
    %% =========================================================== %
