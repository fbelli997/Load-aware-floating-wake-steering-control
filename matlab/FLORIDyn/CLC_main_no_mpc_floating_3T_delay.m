% ======================================================================= %
% CLC_main_floating_v3.m
% Closed-loop control with FLORIS steady-state prediction + validation
%
% Key design choices:
%   - EnKF runs every 12s (wind estimation)
%   - Controller runs every controller_rate seconds (yaw optimization)
%   - Prediction uses runFLORIS with BASE positions (steady-state future)
%   - Validation uses RunMultipleFloris with CURRENT OP state (FLORIDyn
%     instantaneous power, including wake transitient effects)
%   - Reliability: at t+controller_rate, compares FLORIS SS prediction
%     made at t with FLORIDyn instantaneous power at t+controller_rate
%     --> if controller_rate < wake advection time (~400s), P_check will
%         be LOWER than baseline (wake not yet developed) and dP_check
%         will be negative, demonstrating the transitient effect.
% ======================================================================= %

pathToSimulation = ['2026_3T_DEL' filesep '3T_IEA15_floating_00'];

%% Initialize
CLC_init_floating;

%% Add UQ
addpath(genpath("UQ_CLC"))

%% Load Reference data
thrustFF = importSOWFAFile([Sim.PathToSim 'Data' filesep ...
    'FASTfarm_generatedThrust.csv']);
thrustFF(:,3) = thrustFF(:,3)./(10^6);

uSOWFA   = 0;
powSOWFA = 7e+6;

%% Run condition
run_controller = true;

%% Simulation parameters
Sim.StartTime = 0;
Sim.EndTime   = 4600;

%% Time keeper
time.it         = 0;
time.current    = Sim.StartTime;
time.current_nO = 0;
time.end        = Sim.EndTime;
time.stepWorld  = 0.5;
time.stepEnKF   = Sim.TimeStep * EnKF.nS;      % ~12 s

%% ===== CONTROLLER TIMING =====
% Set controller_rate here. Physics motivation:
%   Wake advection time ~ distance / U_inf
%   For 3T at 6D spacing, D=240m: total distance ~1440m
%   At 10 m/s: t_adv ~ 144s per gap, ~430s total
%
%   controller_rate = 48s  --> TRANSITIENT regime (wake not yet developed)
%                              P_check will be NEGATIVE vs baseline
%   controller_rate = 408s --> STEADY-STATE regime (wake fully propagated)
%                              P_check should match FLORIS prediction
%
controller_rate = 34 * time.stepEnKF;   % [s] -- change here to test different regimes

%% ===== OPTIMIZATION SETTINGS =====
opt_type       = 'CoL_DEL';
lambda_blade   = 1;
lambda_tower   = 1;
lambda_P       = 1;
lambda_yaw     = 1;
lambda         = 1;

% Wind uncertainty sampling
N_batch = 10;

% Yaw candidate grid: fractions of GeometricYaw
scale_factors = [0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5];
N_cand = numel(scale_factors);

%% Environmental inputs for DEL models (fixed per simulation)
Hs    = 4;
Tp    = 7;
wvdir = 0;

A = pi*(T.D/2).^2;

%% NN speedup
NN_tower = T.NN_corrector_tower;
NN_blade = T.NN_corrector_blade;
T.NN_corrector_tower = [];
T.NN_corrector_blade = [];
NN_tower_const = parallel.pool.Constant(NN_tower);
NN_blade_const = parallel.pool.Constant(NN_blade);
clear NN_tower NN_blade

%% State tracking
gamma_prev = zeros(1, T.nT);

% Initialize wind estimates
ws_hat = 10;  wd_hat = 270;
ws_std = 0.5; wd_std = 1;
wind_dir_training = 270 - wd_hat;
sin_wd = sin(deg2rad(wind_dir_training));
cos_wd = cos(deg2rad(wind_dir_training));

%% ===== HISTORY ARRAYS =====
ws_hist                = [];
wd_hist                = [];
ws_std_hist            = [];
wd_std_hist            = [];
time_hist              = [];
time_hist_control      = [];
yaw_hist               = [];
DEL_NN_hist_tower      = [];
DEL_NN_hist_blade      = [];
YD_hist                = [];
DD_hist                = [];
TI_hist                = [];
nds_hist               = [];
DEL_tower_bl_hist_full = [];
DEL_blade_bl_hist_full = [];

%% ===== PREDICTION TRACKING =====
pred.valid         = false;
pred.t_applied     = NaN;
pred.gamma_applied = zeros(1, T.nT);
pred.P_pred        = NaN;
pred.P_bl_pred     = NaN;
pred.dP_pred       = NaN;
pred.ws_hat        = NaN;
pred.wd_hat        = NaN;

%% ===== LOGS =====
log_folder = fullfile(Sim.PathToSim, 'Results', 'CLC_plots');
if ~exist(log_folder, 'dir'); mkdir(log_folder); end

log_path = fullfile(log_folder, 'controller_log_loads.csv');
fid_log  = fopen(log_path, 'w');
fprintf(fid_log, ['time,iter,cand_idx,is_selected,opt_type,' ...
    'scale_factor,gam_T1,gam_T2,gam_T3,gamma_prev_T1,gamma_prev_T2,' ...
    'ws_hat,wd_hat,ws_std,wd_std,' ...
    'E_i,std_i,E_bl,std_bl,dP,' ...
    'DEL_tower_sum,DEL_blade_sum,DEL_tower_bl_sum,DEL_blade_bl_sum,' ...
    'dTower,dBlade,DEL_pen,yaw_pen,' ...
    'J_avg,J_RA,J_CoL,J_avg_DEL,J_RA_DEL,J_CoL_DEL\n']);

rel_log_path = fullfile(log_folder, 'reliability_log.csv');
fid_rel      = fopen(rel_log_path, 'w');
fprintf(fid_rel, ['t_applied,t_check,controller_rate,' ...
    'gamma_T1,gamma_T2,gamma_T3,' ...
    'ws_pred,wd_pred,ws_check,wd_check,' ...
    'P_pred_FLORIS,P_bl_pred_FLORIS,dP_pred,' ...
    'P_check_FLORIDyn,P_bl_check_FLORIDyn,dP_check,' ...
    'rel_power,err_dP\n']);

%% ================================================================== %%
%%                         MAIN SIMULATION LOOP                        %%
%% ================================================================== %%
while run_controller
    t_loop = tic;
    time.it = time.it + 1;

    if time.current >= time.end; run_controller = false; end

    %% =================== ESTIMATOR (every 12s) ==================== %%
    if mod(time.current_nO, time.stepEnKF) == 0 && time.current > Sim.StartTime

        fprintf('[CLC_EnKF_Par START]\n');
        tt_enkf = tic;
        CLC_EnKF_Par_floating;
        fprintf('[CLC_EnKF_Par END] (%.3f s)\n', toc(tt_enkf));

        [T, M, stdM] = EnKF_getMeanState(EnKF, T, paramFLORIDyn);

        ws_en  = mean(EnKF.States.Vel(T.StartI,:), 1);
        ws_hat = mean(ws_en);
        ws_std = std(ws_en);
        wd_en  = mean(EnKF.States.Dir(T.StartI,:), 1);
        wd_hat = mean(wd_en);
        wd_std = std(wd_en);

        wind_dir_training = 270 - wd_hat;
        sin_wd = sin(deg2rad(wind_dir_training));
        cos_wd = cos(deg2rad(wind_dir_training));

        T.GeomYaw.ws_hat = ws_hat;
        T.GeomYaw.wd_hat = wd_hat;
        T.GeomYaw.ws_std = ws_std;
        T.GeomYaw.wd_std = wd_std;

        fprintf('[EnKF @ t=%.1f]  WS=%.2f+-%.2f  WD=%.2f+-%.2f\n', ...
            time.current, ws_hat, ws_std, wd_hat, wd_std);

        time_hist(end+1,1)   = time.current;
        ws_hist(end+1,1)     = ws_hat;
        wd_hist(end+1,1)     = wd_hat;
        ws_std_hist(end+1,1) = ws_std;
        wd_std_hist(end+1,1) = wd_std;

        if EnKF.StoreData
            folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName filesep 'EnKF_States'];
            if ~exist(folder_name, 'dir'); mkdir(folder_name); end
            writematrix(T.States_WF(:,1)', [folder_name filesep 'vel'],     'WriteMode','append')
            writematrix(T.States_WF(:,2)', [folder_name filesep 'dir'],     'WriteMode','append')
            writematrix(T.States_Var_Vel', [folder_name filesep 'vel_var'], 'WriteMode','append')
            writematrix(T.States_Var_Dir', [folder_name filesep 'dir_var'], 'WriteMode','append')
        end
    end

    %% ================ CONTROLLER (every controller_rate s) ========= %%
    if mod(time.current_nO, controller_rate) == 0 && time.current > controller_rate

        disp('////////// Run CONTROLLER //////////')
        turb = [1, 2, 3];

        %% ---- STEP 0: VALIDATE PREVIOUS PREDICTION ----
        % Reads FLORIDyn INSTANTANEOUS power using current OP state.
        % If controller_rate < wake advection time, scia has not yet
        % propagated --> P_check_FLORIDyn < P_bl_check --> dP_check < 0.
        % This demonstrates the wake transitient effect and justifies
        % using controller_rate >= advection time in practice.
        if pred.valid
            fprintf('\n--- VALIDATION from t=%.0f (controller_rate=%.0fs) ---\n', ...
                pred.t_applied, controller_rate);

            % --- FLORIDyn instantaneous power WITH applied yaw ---
            % Uses current OP state (scie in transitorio)
            % Must recompute dep and intOPs on current state
            T_inst = T;
            T_inst.States_T(T.StartI, 2) = pred.gamma_applied(:);
            T_inst.dep    = findTurbineGroups(T_inst.nT, T_inst.posBase, T_inst.States_OP, ...
                            T_inst.StartI, T_inst.States_WF, T_inst.D, T_inst.nOP, paramFLORIDyn);
            T_inst.intOPs = interpolateOPs(T_inst.nT, T_inst.StartI, T_inst.dep, ...
                            T_inst.States_OP, T_inst.posBase, T_inst.nOP);

            [~, Cp_inst, Ueff_inst, ~] = RunMultipleFloris(T_inst, paramFLORIS, Wind, turb);
            yaw_corr_inst = cosd(pred.gamma_applied(:)).^paramFLORIS.p_p;
            P_check = sum(0.5 .* paramFLORIS.airDen .* Cp_inst .* A ...
                .* Ueff_inst.^3 .* yaw_corr_inst);

            % --- FLORIDyn instantaneous counterfactual baseline (gamma=0) ---
            % Same OP state but yaw forced to 0 --> true counterfactual
            T_bl_inst = T;
            T_bl_inst.States_T(T.StartI, 2) = 0;
            T_bl_inst.dep    = findTurbineGroups(T_bl_inst.nT, T_bl_inst.posBase, T_bl_inst.States_OP, ...
                               T_bl_inst.StartI, T_bl_inst.States_WF, T_bl_inst.D, T_bl_inst.nOP, paramFLORIDyn);
            T_bl_inst.intOPs = interpolateOPs(T_bl_inst.nT, T_bl_inst.StartI, T_bl_inst.dep, ...
                               T_bl_inst.States_OP, T_bl_inst.posBase, T_bl_inst.nOP);

            [~, Cp_bl_inst, Ueff_bl_inst, ~] = RunMultipleFloris( ...
                T_bl_inst, paramFLORIS, Wind, turb);
            P_bl_check = sum(0.5 .* paramFLORIS.airDen .* Cp_bl_inst .* A ...
                .* Ueff_bl_inst.^3);

            dP_check  = (P_check - P_bl_check) / P_bl_check;
            rel_power = P_check / pred.P_pred;
            err_dP    = dP_check - pred.dP_pred;

            fprintf(fid_rel, ['%.2f,%.2f,%.1f,' ...
                '%.4f,%.4f,%.4f,' ...
                '%.4f,%.4f,%.4f,%.4f,' ...
                '%.6e,%.6e,%.4f,' ...
                '%.6e,%.6e,%.4f,' ...
                '%.6f,%.4f\n'], ...
                pred.t_applied, time.current, controller_rate, ...
                pred.gamma_applied(1), pred.gamma_applied(2), pred.gamma_applied(3), ...
                pred.ws_hat, pred.wd_hat, ws_hat, wd_hat, ...
                pred.P_pred, pred.P_bl_pred, pred.dP_pred*100, ...
                P_check, P_bl_check, dP_check*100, ...
                rel_power, err_dP*100);

            fprintf('  dP_pred  (FLORIS SS):    %.2f%%\n', pred.dP_pred*100);
            fprintf('  dP_check (FLORIDyn OP):  %.2f%%\n', dP_check*100);
            fprintf('  rel_power = %.4f  |  err_dP = %.2f%%\n', rel_power, err_dP*100);
            if dP_check < 0
                fprintf('  --> Wake not yet developed (expected at %.0fs < t_adv)\n', controller_rate);
            end
            fprintf('---------------------------------------\n\n');
        end


        %% ---- STEP 1: YAW CANDIDATES (grid) ----
        [geom_yaw_mean, ~, ~] = GeometricYaw( ...
            T.GeomYaw, T.posBase(:,1), T.posBase(:,2));

        gam = zeros(N_cand, T.nT);
        for k = 1:N_cand
            gam(k,:) = geom_yaw_mean(:)' * scale_factors(k);
        end

        fprintf('Gamma candidati:\n');
        for k = 1:N_cand
            fprintf('  [sf=%.2f] gamma = [%.2f, %.2f, %.2f] deg\n', ...
                scale_factors(k), gam(k,1), gam(k,2), gam(k,3));
        end


        %% ---- STEP 2: SAMPLE WIND SCENARIOS ----
        ws_batch = ws_hat + ws_std * randn(N_batch, 1);
        wd_batch = wd_hat + wd_std * randn(N_batch, 1);


        %% ---- STEP 3: POWER + DEL MATRIX (FLORIS steady-state, posBase) ----
        power_matrix     = zeros(N_cand, N_batch);
        T_Ueff_hist      = zeros(N_cand, N_batch, T.nT);
        TI_eff_hist      = zeros(N_cand, N_batch, T.nT);
        DEL_blade_matrix = zeros(N_cand, T.nT);
        DEL_tower_matrix = zeros(N_cand, T.nT);

        [DD_opt, YD_opt, turb_opt] = LayoutDisp(wd_hat, ...
            T.posBase(:,1), T.posBase(:,2), T.D(1));

        for i = 1:N_cand
            for j = 1:N_batch
                %% Steady-state: runFLORIS turbine by turbine with posBase
                nT_local  = T.nT;
                T_Ueff_ss = zeros(nT_local, 1);
                TI_add_ss = zeros(nT_local, 1);
                Cp_ss     = zeros(nT_local, 1);

                % Sort turbines upstream -> downstream
                phi_rot = angSOWFA2world(wd_batch(j));
                R2d     = [cos(phi_rot), sin(phi_rot); -sin(phi_rot), cos(phi_rot)];
                pos_rot = (R2d * T.posBase(:,1:2)')';
                [~, sort_idx] = sort(pos_rot(:,1), 'ascend');

                for k = 1:nT_local
                    iT       = sort_idx(k);
                    upstream = sort_idx(1:k-1);
                    n_up     = numel(upstream);
                    n_tmp    = n_up + 1;

                    if n_up == 0
                        % Most upstream: free-stream
                        [~, Cp_ss(iT)] = CalcCt_floating(ws_batch(j), T.floater.VCpCt.VCpCt);
                        T_Ueff_ss(iT)  = ws_batch(j);
                        TI_add_ss(iT)  = 0;
                    else
                        % Build arrays for runFLORIS
                        LocationT_tmp = [T.posBase(upstream,:); T.posBase(iT,:)];
                        States_WF_tmp = repmat([ws_batch(j), wd_batch(j), 0.06], n_tmp, 1);
                        States_T_tmp  = [zeros(n_up,1), gam(i,upstream)', zeros(n_up,1);
                                         0,             gam(i,iT),         0            ];
                        D_tmp         = [T.D(upstream); T.D(iT)];
                        ptilt_tmp     = T.floater.theta_abs([upstream(:); iT]);

                        % Ct of upstream turbines at their effective wind speed
                        Ct_up = zeros(n_up, 1);
                        for m = 1:n_up
                            [Ct_up(m), ~] = CalcCt_floating( ...
                                T_Ueff_ss(upstream(m)), T.floater.VCpCt.VCpCt);
                        end

                        [~, T_aTI_arr, T_Ueff_i, ~] = runFLORIS( ...
                            LocationT_tmp, States_WF_tmp, States_T_tmp, D_tmp, ...
                            paramFLORIS, Wind.Shear, Ct_up, ptilt_tmp, k);

                        T_Ueff_ss(iT) = T_Ueff_i;
                        TI_add_ss(iT) = sqrt(sum(T_aTI_arr.^2));
                        [~, Cp_ss(iT)] = CalcCt_floating(T_Ueff_i, T.floater.VCpCt.VCpCt);
                    end
                end

                T_Ueff_hist(i,j,:) = T_Ueff_ss;
                TI_eff_vec         = TI_add_ss + [0.06; 0.03; 0.03];
                TI_eff_hist(i,j,:) = TI_eff_vec;

                yaw_corr = cosd(gam(i,:)').^(paramFLORIS.p_p);
                power_matrix(i,j) = sum( ...
                    0.5 .* paramFLORIS.airDen .* Cp_ss .* A .* T_Ueff_ss.^3 .* yaw_corr);
            end

            %% DEL matrix
            T_Ueff_mean_i = squeeze(mean(T_Ueff_hist(i,:,:), 2));
            TI_eff_mean_i = squeeze(mean(TI_eff_hist(i,:,:), 2));

            for t = 1:T.nT
                if turb_opt(t) == 1
                    yaw0_deg_t = 0;  gamma_up_t = 0;
                else
                    up_t       = find(turb_opt == turb_opt(t)-1, 1);
                    yaw0_deg_t = abs(gam(i, up_t));
                    gamma_up_t = abs(gam(i, up_t));
                end

                [DEL_NN_b, ~, ~] = predictBladeNN( ...
                    T.NN_predictor_blade, gam(i,t), Hs, Tp, wvdir, ...
                    T_Ueff_mean_i(t), TI_eff_mean_i(t), wind_dir_training);

                [DEL_NN_tw, ~, ~] = predictTowerNN( ...
                    T.NN_predictor_tower, gam(i,t), Hs, Tp, wvdir, ...
                    T_Ueff_mean_i(t), TI_eff_mean_i(t), wind_dir_training);
                res_tw = residualTowerNN( ...
                    NN_tower_const.Value, DEL_NN_tw, DD_opt(t), YD_opt(t), ...
                    gamma_up_t, yaw0_deg_t, Hs, Tp, wvdir, ...
                    T_Ueff_mean_i(t), TI_eff_mean_i(t), turb_opt(t), sin_wd, cos_wd);

                DEL_blade_matrix(i,t) = DEL_NN_b;
                DEL_tower_matrix(i,t) = DEL_NN_tw - res_tw;
            end
        end


        %% ---- STEP 4: BASELINE (gamma=0, FLORIS SS) ----
        P_bl           = zeros(N_batch, 1);
        T_Ueff_bl_hist = zeros(N_batch, T.nT);
        DEL_blade_bl   = zeros(T.nT, 1);
        DEL_tower_bl   = zeros(T.nT, 1);

        for j = 1:N_batch
            % Baseline uses same posBase approach for consistency
            nT_local  = T.nT;
            T_Ueff_bl = zeros(nT_local, 1);
            TI_add_bl = zeros(nT_local, 1);
            Cp_bl_j   = zeros(nT_local, 1);

            phi_rot = angSOWFA2world(wd_batch(j));
            R2d     = [cos(phi_rot), sin(phi_rot); -sin(phi_rot), cos(phi_rot)];
            pos_rot = (R2d * T.posBase(:,1:2)')';
            [~, sort_idx_bl] = sort(pos_rot(:,1), 'ascend');

            for k = 1:nT_local
                iT       = sort_idx_bl(k);
                upstream = sort_idx_bl(1:k-1);
                n_up     = numel(upstream);
                n_tmp    = n_up + 1;

                if n_up == 0
                    [~, Cp_bl_j(iT)] = CalcCt_floating(ws_batch(j), T.floater.VCpCt.VCpCt);
                    T_Ueff_bl(iT)    = ws_batch(j);
                    TI_add_bl(iT)    = 0;
                else
                    LocationT_tmp = [T.posBase(upstream,:); T.posBase(iT,:)];
                    States_WF_tmp = repmat([ws_batch(j), wd_batch(j), 0.06], n_tmp, 1);
                    States_T_tmp  = [zeros(n_up, 3); 0, 0, 0];   % gamma=0 for all
                    D_tmp         = [T.D(upstream); T.D(iT)];
                    ptilt_tmp     = T.floater.theta_abs([upstream(:); iT]);

                    Ct_up = zeros(n_up, 1);
                    for m = 1:n_up
                        [Ct_up(m), ~] = CalcCt_floating( ...
                            T_Ueff_bl(upstream(m)), T.floater.VCpCt.VCpCt);
                    end

                    [~, T_aTI_bl, T_Ueff_bl_i, ~] = runFLORIS( ...
                        LocationT_tmp, States_WF_tmp, States_T_tmp, D_tmp, ...
                        paramFLORIS, Wind.Shear, Ct_up, ptilt_tmp, k);

                    T_Ueff_bl(iT) = T_Ueff_bl_i;
                    TI_add_bl(iT) = sqrt(sum(T_aTI_bl.^2));
                    [~, Cp_bl_j(iT)] = CalcCt_floating(T_Ueff_bl_i, T.floater.VCpCt.VCpCt);
                end
            end

            TI_eff_bl_j         = TI_add_bl + [0.06; 0.03; 0.03];
            T_Ueff_bl_hist(j,:) = T_Ueff_bl(:)';
            P_bl(j) = sum(0.5 .* paramFLORIS.airDen .* Cp_bl_j .* A .* T_Ueff_bl.^3);

            for t = 1:T.nT
                [DEL_NN_b_bl, ~, ~] = predictBladeNN( ...
                    T.NN_predictor_blade, 0, Hs, Tp, wvdir, ...
                    T_Ueff_bl(t), TI_eff_bl_j(t), wind_dir_training);
                [DEL_NN_tw_bl, ~, ~] = predictTowerNN( ...
                    T.NN_predictor_tower, 0, Hs, Tp, wvdir, ...
                    T_Ueff_bl(t), TI_eff_bl_j(t), wind_dir_training);
                res_tw_bl = residualTowerNN( ...
                    NN_tower_const.Value, DEL_NN_tw_bl, DD_opt(t), YD_opt(t), ...
                    0, 0, Hs, Tp, wvdir, ...
                    T_Ueff_bl(t), TI_eff_bl_j(t), turb_opt(t), sin_wd, cos_wd);

                DEL_blade_bl(t) = DEL_blade_bl(t) + DEL_NN_b_bl / N_batch;
                DEL_tower_bl(t) = DEL_tower_bl(t) + (DEL_NN_tw_bl - res_tw_bl) / N_batch;
            end
        end

        E_bl   = mean(P_bl);
        std_bl = std(P_bl);

        DEL_tower_bl_hist_full(end+1,:) = DEL_tower_bl.';
        DEL_blade_bl_hist_full(end+1,:) = DEL_blade_bl.';
        DEL_blade_bl_sum = sum(DEL_blade_bl);
        DEL_tower_bl_sum = sum(DEL_tower_bl);


        %% ---- STEP 5: COST FUNCTIONS ----
        J_avg       = zeros(N_cand, 1);
        J_RA        = zeros(N_cand, 1);
        J_CoL       = zeros(N_cand, 1);
        J_avg_DEL   = zeros(N_cand, 1);
        J_RA_DEL    = zeros(N_cand, 1);
        J_CoL_DEL   = zeros(N_cand, 1);
        log_rows    = struct([]);

        for i = 1:N_cand
            P_i   = power_matrix(i, :)';
            E_i   = mean(P_i);
            std_i = std(P_i);

            steered       = 1:2;
            delta_yaw_max = 20;
            yaw_pen = mean(((gam(i,steered) - gamma_prev(steered)) / delta_yaw_max).^2);

            dP = (E_i - E_bl) / E_bl;

            DEL_tower_yawed_sum = sum(DEL_tower_matrix(i,:));
            DEL_blade_yawed_sum = sum(DEL_blade_matrix(i,:));
            dTower  = (DEL_tower_yawed_sum - DEL_tower_bl_sum) / DEL_tower_bl_sum;
            dBlade  = (DEL_blade_yawed_sum - DEL_blade_bl_sum) / DEL_blade_bl_sum;
            DEL_pen = lambda_tower * dTower + lambda_blade * dBlade;

            J_avg(i)     = E_i;
            J_avg_DEL(i) = dP - DEL_pen;

            J_RA(i)      = E_i - lambda * std_i;
            J_RA_DEL(i)  = dP - (lambda * std_i / E_bl) - DEL_pen;

            cov_mat  = cov(P_i, P_bl);
            cov_i_bl = cov_mat(1,2);
            mu_diff  = E_i - E_bl;
            std_diff = sqrt(std_i^2 + std_bl^2 - 2*cov_i_bl);

            if std_diff < 1e-6
                J_CoL(i)     = 0.5;
                J_CoL_DEL(i) = 0.5 + DEL_pen;
            else
                J_CoL(i)     = normcdf(0, mu_diff, std_diff);
                J_CoL_DEL(i) = J_CoL(i) - lambda_P * dP + DEL_pen + lambda_yaw * yaw_pen;
            end

            log_rows(i).E_i           = E_i;
            log_rows(i).std_i         = std_i;
            log_rows(i).dP            = dP;
            log_rows(i).DEL_tower_sum = DEL_tower_yawed_sum;
            log_rows(i).DEL_blade_sum = DEL_blade_yawed_sum;
            log_rows(i).dTower        = dTower;
            log_rows(i).dBlade        = dBlade;
            log_rows(i).DEL_pen       = DEL_pen;
            log_rows(i).yaw_pen       = yaw_pen;
        end


        %% ---- STEP 6: SELECT OPTIMAL ----
        switch opt_type
            case 'avg',     [~, idx_opt] = max(J_avg);     cost_name = 'J_avg';
            case 'RA',      [~, idx_opt] = max(J_RA);      cost_name = 'J_RA';
            case 'CoL',     [~, idx_opt] = min(J_CoL);     cost_name = 'J_CoL';
            case 'avg_DEL', [~, idx_opt] = max(J_avg_DEL); cost_name = 'J_avg_DEL';
            case 'RA_DEL',  [~, idx_opt] = max(J_RA_DEL);  cost_name = 'J_RA_DEL';
            case 'CoL_DEL', [~, idx_opt] = min(J_CoL_DEL); cost_name = 'J_CoL_DEL';
            otherwise
                warning('opt_type non valido, uso avg');
                [~, idx_opt] = max(J_avg); cost_name = 'J_avg [default]';
        end

        gamma_opt = gam(idx_opt, :);

        fprintf('\n=== YAW OPTIMIZATION [%s] ===\n', cost_name);
        fprintf('  controller_rate:   %.0f s\n',   controller_rate);
        fprintf('  Scale factor:      %.2f\n',      scale_factors(idx_opt));
        fprintf('  Gamma opt:         [%.2f, %.2f, %.2f] deg\n', gamma_opt(1), gamma_opt(2), gamma_opt(3));
        fprintf('  BASELINE power:    %.4f MW\n',   E_bl/1e6);
        fprintf('  E[P] predicted:    %.4f MW\n',   mean(power_matrix(idx_opt,:))/1e6);
        fprintf('  dP predicted:      %.2f%%\n',    log_rows(idx_opt).dP*100);
        fprintf('  J_CoL:             %.2f%%\n',    J_CoL(idx_opt)*100);
        fprintf('  J_CoL_DEL:         %.6f\n',      J_CoL_DEL(idx_opt));
        fprintf('  DEL_pen:           %.4f\n',      log_rows(idx_opt).DEL_pen);
        fprintf('================================\n\n');


        %% ---- STEP 7: STORE PREDICTION ----
        pred.valid         = true;
        pred.t_applied     = time.current;
        pred.gamma_applied = gamma_opt;
        pred.P_pred        = mean(power_matrix(idx_opt, :));
        pred.P_bl_pred     = E_bl;
        pred.dP_pred       = (pred.P_pred - pred.P_bl_pred) / pred.P_bl_pred;
        pred.ws_hat        = ws_hat;
        pred.wd_hat        = wd_hat;


        %% ---- STEP 8: WRITE LOG ----
        for i = 1:N_cand
            fprintf(fid_log, ['%.2f,%d,%d,%d,%s,' ...
                '%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,' ...
                '%.4f,%.4f,%.4f,%.4f,' ...
                '%.6e,%.6e,%.6e,%.6e,%.6f,' ...
                '%.6e,%.6e,%.6e,%.6e,' ...
                '%.6f,%.6f,%.6f,%.6f,' ...
                '%.6e,%.6e,%.6f,%.6e,%.6f,%.6f\n'], ...
                time.current, time.it, i, (i==idx_opt), opt_type, ...
                scale_factors(i), gam(i,1), gam(i,2), gam(i,3), gamma_prev(1), gamma_prev(2), ...
                ws_hat, wd_hat, ws_std, wd_std, ...
                log_rows(i).E_i, log_rows(i).std_i, E_bl, std_bl, log_rows(i).dP, ...
                log_rows(i).DEL_tower_sum, log_rows(i).DEL_blade_sum, DEL_tower_bl_sum, DEL_blade_bl_sum, ...
                log_rows(i).dTower, log_rows(i).dBlade, log_rows(i).DEL_pen, log_rows(i).yaw_pen, ...
                J_avg(i), J_RA(i), J_CoL(i), J_avg_DEL(i), J_RA_DEL(i), J_CoL_DEL(i));
        end


        %% ---- STEP 9: APPLY GAMMA OPTIMAL ----
        T.States_T(T.StartI, 2) = gamma_opt(:);

        for iE = 1:EnKF.nE
            col_start = (iE-1)*EnKF.nStatesT + 1;
            col_end   = iE*EnKF.nStatesT;
            States_T_iE = EnKF.States_T(:, col_start:col_end);
            States_T_iE(T.StartI, 2) = gamma_opt(:);
            EnKF.States_T(:, col_start:col_end) = States_T_iE;
        end

        gamma_prev = gamma_opt;


        %% ---- STEP 10: STORE HISTORY ----
        yaw_hist(end+1,:) = gamma_opt;

        DEL_tower_yawed_opt        = DEL_tower_matrix(idx_opt,:);
        DEL_blade_yawed_opt        = DEL_blade_matrix(idx_opt,:);
        DEL_NN_hist_tower(end+1,:) = DEL_tower_yawed_opt;
        DEL_NN_hist_blade(end+1,:) = DEL_blade_yawed_opt;

        if T.GeomYaw.active == 1
            [DD, YD, turb_out] = LayoutDisp(T.States_WF(T.StartI(1), 2), ...
                T.posBase(:,1), T.posBase(:,2), T.D(1));
        else
            DD = 5*ones(T.nT,1); YD = zeros(T.nT,1); turb_out = ones(T.nT,1);
        end

        Mnum_4d   = reshape(EnKFM_out, T.nT, nSimSteps_enkf, 8, EnKF.nE);
        mean_wind = squeeze(mean(mean(Mnum_4d(:,:,4,:), 2), 4));
        TI_wind   = squeeze(mean(mean(Mnum_4d(:,:,3,:), 2), 4))/100 + [0.06, 0.03, 0.03].';

        YD_hist(end+1,:)           = YD.';
        DD_hist(end+1,:)           = DD.';
        TI_hist(end+1,:)           = TI_wind.';
        nds_hist(end+1,:)          = turb_out.';
        time_hist_control(end+1,1) = time.current;

        fprintf('  [controller wall time: %.3f s]\n', toc(t_loop));
    end

    %% TIME STEP
    time.current    = time.current    + time.stepWorld;
    time.current_nO = time.current_nO + time.stepWorld;
end

disp('Simulation done!')
fclose(fid_log);
fclose(fid_rel);
fprintf('Controller log:   %s\n', log_path);
fprintf('Reliability log:  %s\n', rel_log_path);


%% ================================================================== %%
%%                          POST-PROCESSING                            %%
%% ================================================================== %%

if EnKF.StoreData
    folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName filesep 'EnKF_Measurements'];
    if ~exist(folder_name, 'dir'); mkdir(folder_name); end
    writetable(cat(1,EnKF.M{:}), [folder_name filesep 'EnKF_measurements.csv'])
end

if and(EnKF.StoreData, CLC.Con.StoreData)
    postprocess_output([Sim.PathToSim 'Results' filesep CLC.CaseName], T.nT)
end

EnKF_Vis = plotK_CombinedMeasurement(EnKF_Vis, EnKF, T, Sim, paramFLORIS, 5, 1:3);
e = error_EnKF(EnKF, T.nT, thrustFF, uSOWFA, 1000);

plot_folder = fullfile(Sim.PathToSim, 'Results', 'CLC_plots');
if ~exist(plot_folder, 'dir'); mkdir(plot_folder); end
set(0, 'DefaultFigureVisible', 'on');

plot_CLC_results(time_hist, time_hist_control, ws_hist, wd_hist, ws_std_hist, wd_std_hist, ...
    yaw_hist, DEL_NN_hist_blade, DEL_NN_hist_tower, ...
    DD_hist, YD_hist, TI_hist, nds_hist, T, Sim, ...
    DEL_blade_bl_hist_full, DEL_tower_bl_hist_full);


%% ===== RELIABILITY PLOT =====
rel_data = readmatrix(rel_log_path);

if ~isempty(rel_data) && size(rel_data, 1) > 0

    t_applied  = rel_data(:,1);
    t_check    = rel_data(:,2);
    ws_pred    = rel_data(:,7);
    wd_pred    = rel_data(:,8);
    ws_check   = rel_data(:,9);
    wd_check   = rel_data(:,10);
    P_pred_F   = rel_data(:,11) / 1e6;
    P_bl_pred  = rel_data(:,12) / 1e6;
    dP_pred    = rel_data(:,13);
    P_check_F  = rel_data(:,14) / 1e6;
    P_bl_check = rel_data(:,15) / 1e6;
    dP_check   = rel_data(:,16);
    rel_power  = rel_data(:,17);
    err_dP     = rel_data(:,18);

    n_steps = numel(t_applied);
    x_steps = 1:n_steps;

    fig_rel = figure('Name','Reliability Analysis','Position',[50 50 1300 900]);

    %% Panel 1: Power absolute
    subplot(3,2,1)
    hold on
    plot(x_steps, P_pred_F,  'b-o', 'LineWidth',1.5, 'MarkerSize',6, 'DisplayName','P_{pred} FLORIS SS')
    plot(x_steps, P_check_F, 'r-o', 'LineWidth',1.5, 'MarkerSize',6, 'DisplayName','P_{check} FLORIDyn')
    plot(x_steps, P_bl_pred,  'b--','LineWidth',1,   'DisplayName','P_{bl,pred}')
    plot(x_steps, P_bl_check, 'r--','LineWidth',1,   'DisplayName','P_{bl,check}')
    hold off
    xlabel('Controller step [-]'); ylabel('Power [MW]')
    title(sprintf('Power pred vs realized  (rate=%.0fs)', controller_rate))
    legend('Location','best','FontSize',8); grid on

    %% Panel 2: dP bars
    subplot(3,2,2)
    hold on
    bar(x_steps-0.2, dP_pred,  0.35, 'b', 'DisplayName','dP_{pred} FLORIS SS [%]')
    bar(x_steps+0.2, dP_check, 0.35, 'r', 'DisplayName','dP_{check} FLORIDyn [%]')
    yline(0,'k--','LineWidth',1.2)
    hold off
    xlabel('Controller step [-]'); ylabel('Power gain [%]')
    title('Power gain: FLORIS SS pred vs FLORIDyn realized')
    legend('Location','best','FontSize',8); grid on

    %% Panel 3: rel_power
    subplot(3,2,3)
    hold on
    bar(x_steps, rel_power, 'FaceColor',[0.2 0.6 0.8])
    yline(1.0,  'k--','LineWidth',1.5,'DisplayName','Perfect')
    yline(1.05, 'r:' ,'LineWidth',1)
    yline(0.95, 'r:' ,'LineWidth',1)
    hold off
    xlabel('Controller step [-]'); ylabel('rel\_power = P_{check}/P_{pred} [-]')
    title('Power reliability  (1=perfect, <1=wake not yet developed)')
    grid on

    %% Panel 4: err_dP
    subplot(3,2,4)
    hold on
    bar(x_steps, err_dP, 'FaceColor',[0.8 0.4 0.2])
    yline(0,'k--','LineWidth',1.5)
    hold off
    xlabel('Controller step [-]'); ylabel('err\_dP [%]')
    title('Power gain error: dP_{check} - dP_{pred}  (negative = transitient)')
    grid on

    %% Panel 5: Wind speed shift
    subplot(3,2,5)
    hold on
    plot(x_steps, ws_pred,  'b-o','LineWidth',1.5,'DisplayName','WS @ t_{applied}')
    plot(x_steps, ws_check, 'r-o','LineWidth',1.5,'DisplayName','WS @ t_{check}')
    hold off
    xlabel('Controller step [-]'); ylabel('Wind speed [m/s]')
    title('Wind speed: at decision vs at validation')
    legend('Location','best','FontSize',8); grid on

    %% Panel 6: Wind direction shift
    subplot(3,2,6)
    hold on
    plot(x_steps, wd_pred,  'b-o','LineWidth',1.5,'DisplayName','WD @ t_{applied}')
    plot(x_steps, wd_check, 'r-o','LineWidth',1.5,'DisplayName','WD @ t_{check}')
    hold off
    xlabel('Controller step [-]'); ylabel('Wind direction [deg]')
    title('Wind direction: at decision vs at validation')
    legend('Location','best','FontSize',8); grid on

    sgtitle(sprintf(['Reliability Analysis  |  controller\_rate = %.0f s\n' ...
        'Mean rel\_power = %.3f  |  Mean |err\_dP| = %.2f%%'], ...
        controller_rate, mean(rel_power), mean(abs(err_dP))), 'FontSize',12)

    saveas(fig_rel, fullfile(plot_folder, 'reliability_analysis.png'))

    fprintf('\n===== RELIABILITY SUMMARY =====\n')
    fprintf('  controller_rate:   %.0f s\n',  controller_rate)
    fprintf('  N validation pts:  %d\n',       n_steps)
    fprintf('  Mean rel_power:    %.4f\n',     mean(rel_power))
    fprintf('  Std  rel_power:    %.4f\n',     std(rel_power))
    fprintf('  Mean err_dP:       %.2f%%\n',   mean(err_dP))
    fprintf('  Mean |err_dP|:     %.2f%%\n',   mean(abs(err_dP)))
    fprintf('  N steps dP<0:      %d / %d\n',  sum(dP_check<0), n_steps)
    fprintf('================================\n')
else
    fprintf('No reliability data (need >= 2 controller steps).\n')
end

%% Save all figures
figList = findall(0,'Type','figure');
for iFig = 1:numel(figList)
    saveas(figList(iFig), fullfile(plot_folder, ['fig_' num2str(figList(iFig).Number) '.png']));
end
disp('All plots saved.')