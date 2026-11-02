function state = CLC_init_floating_fun_ZMQ(WS_init, WD_init)

    % Assegna WS e WD nel base workspace prima di eseguire lo script
    assignin('base', 'WS_init', WS_init);
    assignin('base', 'WD_init', WD_init);
    
    % Esegui CLC_init_floating originale nel base workspace
    evalin('base', "cd('D:\zero_model\FLORIDyn_Matlab-main_3_floating\FLORIDyn_Matlab-main_3')");
    evalin('base', "pathToSimulation = ['2026_3T_DEL' filesep '3T_IEA15_floating_00'];");
    evalin('base', "CLC_init_floating_ZMQ");

    % Recupera variabili dal base workspace
    T             = evalin('base', 'T');
    EnKF          = evalin('base', 'EnKF');
    Sim           = evalin('base', 'Sim');
    Con           = evalin('base', 'Con');
    Wind          = evalin('base', 'Wind');
    Vis           = evalin('base', 'Vis');
    CLC           = evalin('base', 'CLC');
    paramFLORIS   = evalin('base', 'paramFLORIS');
    paramFLORIDyn = evalin('base', 'paramFLORIDyn');
    EnKF_Vis      = evalin('base', 'EnKF_Vis');
    turbProp      = evalin('base', 'turbProp');

    % Pacchetta state
    state.T             = T;
    state.EnKF          = EnKF;
    state.Sim           = Sim;
    state.Con           = Con;
    state.Wind          = Wind;
    state.Vis           = Vis;
    state.CLC           = CLC;
    state.paramFLORIS   = paramFLORIS;
    state.paramFLORIDyn = paramFLORIDyn;
    state.EnKF_Vis      = EnKF_Vis;
    state.NN_tower_const = parallel.pool.Constant(T.NN_corrector_tower);
    state.NN_blade_const = parallel.pool.Constant(T.NN_corrector_blade);

    state.time.it           = 0;
    state.time.current      = 0;
    state.time.current_nO   = 0;
    state.time.stepWorld    = 0.5;
    state.time.stepEnKF     = Sim.TimeStep * EnKF.nS;
    state.time.stepMPC      = Sim.TimeStep * CLC.Time.nS;

    state.pred.valid         = false;
    state.pred.t_applied     = NaN;
    state.pred.gamma_applied = zeros(1, state.T.nT);
    state.pred.P_pred        = NaN;
    state.pred.P_bl_pred     = NaN;
    state.pred.dP_pred       = NaN;
    state.pred.ws_hat        = NaN;
    state.pred.wd_hat        = NaN;
    state.pred.checked       = false;

    state.gamma_prev        = zeros(1, T.nT);
    state.yaw_hist          = [];
    state.ws_hist           = [];
    state.wd_hist           = [];
    state.enkf_log          = [];
    state.time_hist         = [];
    state.time_hist_control = [];
    state.DEL_NN_hist_tower = [];
    state.DEL_NN_hist_blade = [];
    state.feat_hist         = [];
    state.DEL_blade_yawed_opt = [];
    state.DEL_tower_yawed_opt = [];
    state.t_last_opt   = -inf;  
    state.rel_hist     = [];    
    state.wd_buffer      = [];    % rolling buffer of last 3 wd estimates
    state.wd_at_last_opt = NaN;   % wd at last optimization (for drift detection)

    state.opt_type     = 'CoL_DEL';
    state.lambda_blade = 1;
    state.lambda_tower = 1;
    state.lambda_yaw   = 1;
    state.lambda_P     = 5;
    state.lambda       = 1;
    state.N_batch      = 10;
    state.Hs           = 2;
    state.Tp           = 7;
    state.wvdir        = 0;
    state.TI_ambient   = Wind.TI;  % TI ambiente reale, da WindTIConstant.csv (non da turbProp.Init_States, che è la TI aggiunta di bootstrap di States_T)
    state.A            = pi*(T.D/2).^2;
    state.scale_factors = [0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0];
    state.controller_rate = 276; 

    disp('[FLORIDyn_init] Done.')
end