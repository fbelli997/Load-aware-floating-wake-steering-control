function [T, x_disp, y_disp, x_dot, y_dot, theta, thetad, Faero] = ...
    DispActiveFloat_ss(T, tmpM, Sim, paramFLORIDyn, paramFLORIS, it, Ct)


%% Current platform states

s      = T.floater.x(:,1);
h      = T.floater.z(:,1);
theta  = T.floater.theta(:,1);
tfa    = T.floater.tfa(:,1);
sd     = T.floater.vx(:,1);
hd     = T.floater.vz(:,1);
thetad = T.floater.theta_dot(:,1);
tfad   = T.floater.tfa_dot(:,1);

% Sway remains unchanged in this 2-DOF model

y_disp = T.floater.y(:,1);
y_dot  = T.floater.vy(:,1);


Ass = T.floater.StateSpace.A;
Bss = T.floater.StateSpace.B;

addbquad = [9.23e5    0        -8.92e6;
            0         2.30e6    0;
           -8.92e6    0         1.68e10];

%% Debug: verify eigenvalues
debug = true;
if debug
    % % M-R-K natural frequencies
    % [~, Deig] = eig(Keff2, M2);
    % fn_mrk = sqrt(diag(Deig)) / (2*pi);
    lambda = eig(Ass);
    fn_ss  = abs(imag(lambda)) / (2*pi);
    zeta   = -real(lambda) ./ abs(lambda);  % damping ratio
    
    fprintf('\n--- STATE-SPACE DEBUG (senza generatore) ---\n');
    fprintf('Dimensione A: %dx%d\n', size(Ass));
    fprintf('Dimensione B: %dx%d\n', size(Bss));
    fprintf('\nAutovalori:\n');
    for i = 1:length(lambda)
        fprintf('  lambda_%d = %+.4e %+.4ei   f=%.4f Hz   zeta=%.4f\n', ...
            i, real(lambda(i)), imag(lambda(i)), fn_ss(i), zeta(i));
    end
    fprintf('\nFrequenze ordinate [Hz]:');
    fprintf(' %.4f', sort(fn_ss));
    fprintf('\n');
    fprintf('Target: surge ~0.0077 Hz, pitch ~0.036 Hz\n');
    fprintf('------------------------------------------\n');
end



%% Rotor area

A = pi*(T.D/2).^2;



%% Hydrodynamic quadratic drag parameters

rho_water   = 1025;      % water density [kg/m3]
Cd_surge    = 0.7;       % drag coefficient [-]

D_col       = 12.0;      % column diameter [m]
draft       = 26.5;      % draft [m]
n_col       = 4;

A_wet_total = n_col*D_col*draft;



%% Time integration parameters

dt     = Sim.TimeStep;
Nsub   = 1;
dt_sub = dt/Nsub;


idxT = T.StartI;



%% Wind speed and hub height

Ue = tmpM(:,3);

z_hub = T.posNac(:,3);



%% Allocate aerodynamic force output

Faero = zeros(T.nT,2);



%% Integrate nonlinear platform dynamics

for iT = 1:T.nT


    % Initial state

    x = [s(iT);
         h(iT);
         theta(iT);
         tfa(iT);
         sd(iT);
         hd(iT);
         thetad(iT);
         tfad(iT)];


    for k = 1:Nsub


        % Current surge velocity

        surge_speed = x(5);



        % Relative wind velocity seen by rotor

        Urel = Ue(iT)-surge_speed;

        % Avoid negative relative wind speed
        Urel = max(Urel,0);



        % Nonlinear aerodynamic thrust

        FT_iT = 0.5 * paramFLORIS.airDen ...
             * A(iT) * Ct(iT) * Urel^2;



        % Project thrust according to yaw angle

        yaw_iT = deg2rad(T.States_T(idxT(iT),2));


        Faero_iT = AeroForce(FT_iT,yaw_iT);


        Faero(iT,:) = Faero_iT;


        surge_speed = x(5);
        heave_speed = x(6);
        pitch_speed = x(7);
        
        vel = [surge_speed; heave_speed; pitch_speed];
        F_damp = -addbquad * (abs(vel) .* vel);

        % u = [FT_iT;
        %     F_damp(1);
        %     F_damp(2);
            % F_damp(3)];
        u = [2e6;
            0;
            0;
            0];

        x = Ass*x + Bss*u;

    end



    % Store updated states
    
    s(iT)      = x(1);    % surge
    h(iT)      = x(2);    % heave
    theta(iT)  = x(3);    % pitch
    tfa(iT)    = x(4);    % tower
    sd(iT)     = x(5);    % surge vel
    hd(iT)     = x(6);    % heave vel
    thetad(iT) = x(7);    % pitch vel
    tfad(iT)   = x(8);    % tower vel

end



%% Update floater states


%% Update floater states
T.floater.x(:,1)          = s;
T.floater.z(:,1)          = h;       
T.floater.vx(:,1)         = sd;
T.floater.vz(:,1)         = hd;   
T.floater.theta(:,1)      = theta;
T.floater.theta_dot(:,1)  = thetad;
T.floater.tfa(:,1)        = tfa;     
T.floater.tfa_dot(:,1)    = tfad;    
T.floater.y(:,1)           = y_disp;
T.floater.vy(:,1)          = y_dot;


%% Pitch-induced nacelle displacement

x_disp_pitch = z_hub.*sin(theta);

x_dot_pitch = (z_hub.*thetad).*cos(theta);


T.floater.pitch_ws = x_dot_pitch;



%% Update platform position

wd = T.States_WF(T.StartI,2);

alpha = mod(wd-270,360);



T.posBase(:,1) = T.pos0(:,1) ...
    +(T.floater.x(:,1)+x_disp_pitch).*cosd(alpha) ...
    -T.floater.y(:,1).*sind(alpha);



T.posBase(:,2) = T.pos0(:,2) ...
    +T.floater.y(:,1).*cosd(alpha) ...
    +(T.floater.x(:,1)+x_disp_pitch).*sind(alpha);



%% Outputs


x_disp = T.floater.x(:,1)+x_disp_pitch;

x_dot  = T.floater.vx(:,1)+x_dot_pitch;


y_disp = T.floater.y(:,1);

y_dot  = T.floater.vy(:,1);



%% Update observation point height

T.States_OP(T.StartI,3)=z_hub.*cos(theta);


end