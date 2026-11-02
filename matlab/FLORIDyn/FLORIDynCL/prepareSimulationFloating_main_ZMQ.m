function [T, Wind, Sim, Con, paramFLORIS, Vis] = prepareSimulationFloating_main_ZMQ(...
    Wind, Con, paramFLORIDyn, paramFLORIS, turbProp, Sim, Vis, EnKF, WS_init, WD_init)
% Versione semplificata per closed-loop ZMQ
% WS_init: wind speed iniziale da FAST.Farm [m/s]
% WD_init: wind direction iniziale da FAST.Farm [deg, convenzione FLORIDyn]

%% Wind — valori costanti da FAST.Farm
Wind.Dir = [0, WD_init; 1e9, WD_init]; % grande serie temporale per non arrivare al bound finale
Wind.Vel = WS_init;  
Wind.TI  = readmatrix('WindTIConstant.csv');
Wind.Shear.alpha = readmatrix('WindShearPowerLaw.csv');
Wind.Shear.z0    = 1;

%% Controller — yaw costante a 270 (neutro, FLORIDyn gestisce)
Con.YawData = WD_init;  % orientazione iniziale = direzione vento

%% Turbine struct
T.posBase = turbProp.Pos;
T.nT      = size(turbProp.Pos, 1);
t_data    = getTurbineData(turbProp.Type);
T.posNac  = t_data.NacPos;
T.D       = t_data.D;

states = States();

%% Floating features
if paramFLORIDyn.floating == true
    T.floater_use = true;
    T.floater.VCpCt           = turbProp.VCpCt;
    T.floater.MatFloat        = turbProp.MatFloat;
    T.NN_predictor_tower      = turbProp.NN_predictor_tower;
    T.NN_corrector_tower     = turbProp.NN_corrector_tower;
    T.NN_predictor_blade      = turbProp.NN_predictor_blade;
    T.NN_corrector_blade     = turbProp.NN_corrector_blade;
    T.floater.x               = zeros(T.nT, 1);
    T.floater.y               = zeros(T.nT, 1);
    T.floater.vx              = zeros(T.nT, 1);
    T.floater.vy              = zeros(T.nT, 1);
    T.floater.theta_dot       = zeros(T.nT, 1);
    T.floater.theta           = zeros(T.nT, 1);
    T.floater.theta_abs       = zeros(T.nT, 1);
    T.floater.DEL_tower       = zeros(T.nT, 1);
    T.floater.DEL_blade       = zeros(T.nT, 1);
    T.floater.States_names    = {'x_disp','y_disp','x_dot','y_dot','theta_disp','theta_dot','theta_abs'};
    T.floater.theta_undisturbed = ones(T.nT, 1) * turbProp.ptfm_initial;

    %% VOLTURN US MATRICES
    T.floater.VolturnUS.M = [2.01e7   -4.66e7;
                            -4.66e7    4.26e10];

    T.floater.VolturnUS.B = [1.6299e5   0;
                            0          1.2627e9];

    T.floater.VolturnUS.K = [7.192e4    1.145e6;
                            1.145e6    2.587e8];
    %% VOLTURN US MATRICES


    if paramFLORIDyn.CLC == true
        nStates = numel(T.floater.States_names);
        T.floater.T_EnKF_States = zeros(T.nT, nStates, EnKF.nE);
    end

    if ~isfield(T, 'pos0'); T.pos0 = T.posBase; end

    T.GeomYaw.active         = true;
    T.GeomYaw.ws_rated       = 10.59;
    T.GeomYaw.diameter       = T.D(1);
    T.GeomYaw.ws_hat         = [];
    T.GeomYaw.wd_hat         = [];
    T.GeomYaw.ws_std         = [];
    T.GeomYaw.wd_std         = [];
    T.GeomYaw.ws_hat_batch   = [];
    T.GeomYaw.wd_hat_batch   = [];
    T.GeomYaw.k_yaw          = 0.07;
    T.GeomYaw.yaw_max        = 20.;
    T.GeomYaw.yaw_m_x        = 0.2651;
    T.GeomYaw.yaw_m_y        = 15.00;
    T.GeomYaw.yaw_m_center_x = 1.0e-06;
    T.GeomYaw.yaw_c_center_y = 0.0;
    T.GeomYaw.k_farmEqFunc   = 0.1;
    T.GeomYaw.a_farmEqFunc   = 0.0;
    T.GeomYaw.c_farmEqFunc   = 1.0;

    T.floater.M_tot = 20.093e6;
    T.floater.M_add = 9.64e6;
    T.floater.M_eff = T.floater.M_tot + T.floater.M_add;
    T.floater.K_eff = 1.50e5;
    T.floater.C_eff = 9.225e5;
    T.floater.I_pitch = 1.695e10;
    T.floater.C_pitch = 1.676e10;
    T.floater.K_pitch = 2.190e9;
end

%% Stati OP
T.States_OP = zeros(paramFLORIDyn.n_op * T.nT, states.OP);
T.Names_OP  = states.OP_names;
T.States_T  = zeros(paramFLORIDyn.n_op * T.nT, states.Turbine);
T.Names_T   = states.T_names;
T.States_WF = zeros(paramFLORIDyn.n_op * T.nT, states.WF);
T.Names_WF  = states.WF_names;
T.StartI    = 1:paramFLORIDyn.n_op:paramFLORIDyn.n_op * T.nT;
T.nOP       = paramFLORIDyn.n_op;
T.red_arr   = ones(T.nT, T.nT);

if ~isfield(paramFLORIDyn, 'deltaUW')
    paramFLORIDyn.deltaUW = paramFLORIDyn.deltaDW;
end

%% Inizializza stati con WS e WD da FAST.Farm
[T.States_OP, T.States_T, T.States_WF] = InitStates(...
    T, Wind, turbProp.Init_States, paramFLORIS, Sim);

%% Sim
Sim.nSimSteps = length(Sim.StartTime:Sim.TimeStep:Sim.EndTime);
paramFLORIS.RotorPoints = Sim.RotorPoints;

%% Vis — disabilitata per closed-loop
Vis.Film.InProgress = false;

end