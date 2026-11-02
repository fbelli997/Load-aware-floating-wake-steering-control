% Copyright (C) <2024>, M Becker
%
% List of the contributors to the development of FLORIDyn: see LICENSE file.
% Description and complete License: see LICENSE file.
	
% This program (FLORIDyn) is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.

% You should have received a copy of the GNU Affero General Public License
% along with this program (see COPYING file).  If not, see <https://www.gnu.org/licenses/>.
% ======================================================================= %
% Updated: 16. Dez. 2024, M. Becker
% ======================================================================= %

%% Closed loop control main funtion
% Combines the EnKF with the MPC, both using FLORIDyn

pathToSimulation = ['2026_3T_DEL' filesep '3T_IEA15_floating_00'];
% pathToSimulation = ['2026_6T_Fede_floating' filesep '6T_IEA15_floating_00'];
% pathToSimulation = ['2024_CLC' filesep '10T_HKN_12_to_15_deg'];


%% Initialize
CLC_init_floating;

%% Add UQ
addpath(genpath("UQ_CLC"))

%% Load Reference data (temporary)

thrustFF = importSOWFAFile([Sim.PathToSim 'Data' filesep ...
    'FASTfarm_generatedThrust.csv']);
% thrustFF(:,3) = thrustFF(:,3)./(10^6*paramFLORIS.airDen);
thrustFF(:,3) = thrustFF(:,3)./(10^6);


uSOWFA = 0;
powSOWFA = 7e+6;

%% Run
% Run condition
run_controller  = true;

%% parametes changed sepcifically for this simulation   
Sim.StartTime = 0;
Sim.EndTime = 4600;

% "Global" time keeper
time.it         = 0;                % Current iteration
time.current    = Sim.StartTime;    % World time (LES + precursor time) 
time.current_nO = 0;                % World time reset to 0
time.end        = Sim.EndTime;      % World time to end controller
time.stepWorld  = 0.5;              % World / LES time step 
time.stepEnKF   = Sim.TimeStep * EnKF.nS;       % Time step of the EnKF
time.stepMPC    = Sim.TimeStep * CLC.Time.nS;   % Time step of the MPC

ws_hist             = [];
wd_hist             = [];
ws_std_hist         = [];
wd_std_hist         = [];
time_hist           = [];
time_hist_control   = [];
yaw_hist            = [];          
DEL_NN_hist_tower   = [];
DEL_NN_hist_blade   = [];
YD_hist             = [];
DD_hist             = [];
TI_hist             = [];
nds_hist            = [];
DEL_tower_bl_hist_full = [];
DEL_blade_bl_hist_full = [];
dP_hist                = [];

    
% Environmental / operating inputs for DEL models
Hs        = 4;
Tp        = 7;
wvdir      = 0;


A = pi*(T.D/2).^2;

%% ===== OPTIMIZATION SETTINGS =====
opt_type = 'CoL_DEL';  % 'avg', 'RA', 'CoL' 'avg_DEL', 'RA_DEL', 'CoL_DEL'
lambda_blade   = 1;
lambda_tower   = 1;
lambda_P       = 1;
lambda_yaw     = 1;

% For the average cost function
lambda         = 1;
%Number of wind conditions generated
N_batch        = 10;

%% ********************** SPEEDUP **********************
NN_tower = T.NN_corrector_tower;
NN_blade = T.NN_corrector_blade;
T.NN_corrector_tower = [];
T.NN_corrector_blade = [];

NN_tower_const = parallel.pool.Constant(NN_tower);
NN_blade_const = parallel.pool.Constant(NN_blade);
clear NN_tower NN_blade
%% ********************** SPEEDUP **********************

gamma_prev = zeros(1, T.nT);

%% ===== CONTROLLER LOG =====
log_folder = fullfile(Sim.PathToSim, 'Results', 'CLC_plots');
if ~exist(log_folder, 'dir'); mkdir(log_folder); end
log_path = fullfile(log_folder, 'controller_log_loads.csv');
fid_log = fopen(log_path, 'w');
fprintf(fid_log, ['time,iter,cand_idx,is_selected,opt_type,' ...
    'gam_T1,gam_T2,gam_T3,gamma_prev_T1,gamma_prev_T2,' ...
    'ws_hat,wd_hat,ws_std,wd_std,' ...
    'E_i,std_i,E_bl,std_bl,dP,' ...
    'DEL_tower_sum,DEL_blade_sum,DEL_tower_bl_sum,DEL_blade_bl_sum,' ...
    'dTower,dBlade,DEL_pen,yaw_pen,' ...
    'J_avg,J_RA,J_CoL,J_avg_DEL,J_RA_DEL,J_CoL_DEL\n']);
log_rows = struct([]);


while run_controller
    % Increase counter & check break condition
    t_loop = tic;
    % disp(['t = ' num2str(time.current,7)])
    time.it = time.it + 1;

    if time.current >= time.end; run_controller = false; end % Simulation end
    % fprintf('  [pre-EnKF check: %.4f s]\n', toc(t_loop));

    %% Estimator
    if mod(time.current_nO, time.stepEnKF) == 0 && time.current > Sim.StartTime
        disp('////////// Run EnKF //////////')
        % ======= Run Ensemble Kalman Filter ======= 
        
        fprintf('[CLC_EnKF_Par START]\n');
        tt_enkf_wall = tic;
        
        CLC_EnKF_Par_floating;
    
        t_EnKF = toc(tt_enkf_wall);
        fprintf('[CLC_EnKF_Par END]   (%.6f s)\n', t_EnKF);
        fprintf('  [post-EnKF: %.4f s]\n', toc(t_loop));
        ws_hat = mean(EnKF.States.Vel(T.StartI, :), 2);
        wd_hat = mean(EnKF.States.Dir(T.StartI, :), 2);

        % --- EnKF estimated wind at turbine locations (StartI) ---
        fprintf('[EnKF @ t=%.2f]\n', time.current);
        for k = 1:numel(ws_hat)
            fprintf('  T%02d: WS_hat = %.3f m/s, WD_hat = %.2f deg\n', k, ws_hat(k), wd_hat(k));
        end

        Uinfty = mean(T.States_WF(T.StartI,1));
        WDinfty = mean(T.States_WF(T.StartI,2));

        % ======= Collect "true state" from ensembles
        [T,M, stdM] = EnKF_getMeanState(EnKF, T, paramFLORIDyn);

        ws_en = mean(EnKF.States.Vel(T.StartI,:),1);
        ws_hat = mean(ws_en);
        ws_std = std(ws_en);
        wd_en = mean(EnKF.States.Dir(T.StartI,:),1);
        wd_hat = mean(wd_en);
        wd_std = std(wd_en);
        wind_dir_training = 270 - wd_hat;
        sin_wd = sin(deg2rad(wind_dir_training));
        cos_wd = cos(deg2rad(wind_dir_training));
        T.GeomYaw.ws_hat = ws_hat;
        T.GeomYaw.wd_hat = wd_hat;
        T.GeomYaw.ws_std = ws_std;
        T.GeomYaw.wd_std = wd_std;


        % ======= Store history of scalar wind estimates =======
        time_hist(end+1,1)   = time.current;
        ws_hist(end+1,1)     = ws_hat;
        wd_hist(end+1,1)     = wd_hat;
        ws_std_hist(end+1,1) = ws_std;
        wd_std_hist(end+1,1) = wd_std;

        fprintf('predicted wind speed:   = %.6f\n', ws_hat);
        fprintf('uncertainty on ws prediciton: %.6f\n', ws_std);
        fprintf('predicted wind direction:   = %.6f\n', wd_hat);
        fprintf('uncertainty on wd prediciton: %.6f\n', wd_std);

        if EnKF.StoreData
            folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName ...
                filesep 'EnKF_States'];
            if ~exist(folder_name, "dir")
                mkdir(folder_name)
            end
            writematrix(T.States_WF(:,1)',[folder_name filesep...
                'vel'],'WriteMode','append')
            writematrix(T.States_WF(:,2)',[folder_name filesep...
                'dir'],'WriteMode','append')
            writematrix(T.States_Var_Vel',[folder_name filesep...
                'vel_var'],'WriteMode','append')
            writematrix(T.States_Var_Dir',[folder_name filesep...
                'dir_var'],'WriteMode','append')
        end

    end

    
    controller_rate = 1 * time.stepEnKF;

    if mod(time.current_nO, time.stepEnKF) == 0 && time.current > Sim.StartTime
        %% BATCH
        ws_batch = T.GeomYaw.ws_hat + T.GeomYaw.ws_std * randn(N_batch, 1);
        wd_batch = T.GeomYaw.wd_hat + T.GeomYaw.wd_std * randn(N_batch, 1);
        wind_cond = [ws_batch, wd_batch];
        
        %% STEP 1: calcola gamma candidati (un vettore yaw per ogni scenario)
        gam = zeros(N_batch, T.nT);
        
        for i = 1:N_batch
            T.GeomYaw.wd_hat= wd_batch(i);
            T.GeomYaw.ws_hat = ws_batch(i);
        
            [geom_yaw_i, ~, ~] = GeometricYaw( ...
                T.GeomYaw, T.posBase(:,1), T.posBase(:,2));
        
            gam(i,:) = geom_yaw_i(:)';


        end

        T.GeomYaw.wd_hat = wd_hat;
        T.GeomYaw.ws_hat = ws_hat; 

        % % ===== Aggiungi candidato yaw=0 come extra =====
        % gam(end+1, :) = 0;
        % N_cand = size(gam, 1);   % = N_batch + 1

        % % ********************* TEST HALFED YAW ***********************
        % gam(:, 2) = gam(:, 2) * 0.5;  % dimezza yaw turbina 2
        
        % ================RIATTIVA PER VEDERE SET DI YAW==================
        fprintf('Gamma candidati:\n');
        for i = 1:N_batch
            fprintf('  Scenario %02d: gamma = %.2f deg\n', i, gam(i,1));
        end

    

        %% STEP 2: power matrix (N_batch x N_batch)
        
        power_matrix     = zeros(N_batch, N_batch);
        clos             = zeros(N_batch, N_batch, T.nT);
        TI_eff_hist      = zeros(N_batch, N_batch, T.nT);
        DEL_blade_matrix = zeros(N_batch, T.nT);
        DEL_tower_matrix = zeros(N_batch, T.nT);
        turb = [1, 2, 3];
        [DD_opt, YD_opt, turb_opt] = LayoutDisp(wd_hat, T.posBase(:,1), T.posBase(:,2), T.D(1));

        for i = 1:N_batch       % candidato yaw i
            for j = 1:N_batch   % scenario di vento j

                T_tmp = T;
                T_tmp.States_T(T.StartI, 2) = gam(i,:)';
                T_tmp.States_WF(:, 1) = ws_batch(j);  
                T_tmp.States_WF(:, 2) = wd_batch(j);  
                T_tmp.dep = findTurbineGroups(T_tmp.nT, T_tmp.posBase, T_tmp.States_OP , T_tmp.StartI, T_tmp.States_WF, T_tmp.D, T_tmp.nOP, paramFLORIDyn);
                T_tmp.intOPs = interpolateOPs(T_tmp.nT, T_tmp.StartI, T_tmp.dep, T_tmp.States_OP, T_tmp.posBase, T_tmp.nOP);
                
                [~, Cp,T_Ueff, TI_eff_vec] = RunMultipleFloris( ...
                    T_tmp, paramFLORIS, Wind, turb);
                
                T_Ueff_hist(i,j,:) = T_Ueff;
                TI_eff_vec = TI_eff_vec + [0.06, 0.03, 0.03]';
                TI_eff_hist(i,j,:) = TI_eff_vec;


                yaw_corr  = cosd(gam(i,:)').^(paramFLORIS.p_p);            % [nT x 1] - tutte le turbine del candidato i
                                
                power_matrix(i,j) = sum( ...
                    0.5 .* paramFLORIS.airDen .* Cp .* A .* T_Ueff.^3 ...
                    .* yaw_corr);

            end

            % DEL matrix calculation:
            T_Ueff_mean_i = squeeze(mean(T_Ueff_hist(i,:,:), 2));  % [nT x 1]
            TI_eff_mean_i = squeeze(mean(TI_eff_hist(i,:,:), 2));  % [nT x 1]
            
            for t = 1:T.nT
                T_Ueff_t = T_Ueff_mean_i(t);
            
                if turb_opt(t) == 1
                    yaw0_deg_t = 0;
                    gamma_up_t = 0;
                else
                    up_t = find(turb_opt == turb_opt(t)-1, 1);
                    yaw0_deg_t = abs(gam(i, up_t));
                    gamma_up_t = abs(gam(i, up_t));
                end
                  

                [DEL_NN_b, ~, ~]  = predictBladeNN(T.NN_predictor_blade, gam(i,t), Hs, Tp, wvdir, T_Ueff_t, TI_eff_mean_i(t), wind_dir_training);
                % res_b = residualBladeNN(NN_blade_const.Value, DEL_NN_b, DD_opt(t), YD_opt(t), gamma_up_t, yaw0_deg_t, Hs, Tp, wvdir, T_Ueff_t, TI_eff_mean_i(t), turb_opt(t), sin_wd, cos_wd);
            
                [DEL_NN_tw, ~, ~] = predictTowerNN(T.NN_predictor_tower, gam(i,t), Hs, Tp, wvdir, T_Ueff_t, TI_eff_mean_i(t), wind_dir_training);
                res_tw = residualTowerNN(NN_tower_const.Value, DEL_NN_tw, DD_opt(t), YD_opt(t), gamma_up_t, yaw0_deg_t, Hs, Tp, wvdir, T_Ueff_t, TI_eff_mean_i(t), turb_opt(t), sin_wd, cos_wd);
            
                DEL_blade_matrix(i,t) = DEL_NN_b; %- res_b;
                DEL_tower_matrix(i,t) = DEL_NN_tw - res_tw;
            end
        end 
    end
    
    %% CONTROLLER:
    % Default yaw command
     if mod(time.current_nO, controller_rate) == 0 && time.current > controller_rate

        disp('////////// Run CONTROLLER and DEL calculator //////////')   
        geom_yaw = zeros(T.nT,1);
        upstream_idx = zeros(T.nT,1);


        %% STEP 3: BASELINE (gamma = 0 per tutte le turbine)

        P_bl            = zeros(N_batch, 1);
        T_Ueff_bl_hist  = zeros(N_batch, T.nT);
        DEL_blade_bl    = zeros(T.nT, 1);
        DEL_tower_bl    = zeros(T.nT, 1);

        for j = 1:N_batch
            T_tmp_b = T;
            T_tmp_b.States_T(:, 2) = 0;  % gamma = 0
            T_tmp_b.States_WF(:, 1) = ws_batch(j);
            T_tmp_b.States_WF(:, 2) = wd_batch(j);
            T_tmp_b.dep = findTurbineGroups(T_tmp_b.nT, T_tmp_b.posBase, T_tmp_b.States_OP , T_tmp_b.StartI, T_tmp_b.States_WF, T_tmp_b.D, T_tmp_b.nOP, paramFLORIDyn);
            T_tmp_b.intOPs = interpolateOPs(T_tmp_b.nT, T_tmp_b.StartI, T_tmp_b.dep, T_tmp_b.States_OP, T_tmp_b.posBase, T_tmp_b.nOP);
            
            [~, Cp_bl, T_Ueff_bl, TI_eff_bl] = RunMultipleFloris( ...
                T_tmp_b, paramFLORIS, Wind, turb);

            TI_eff_bl = TI_eff_bl(:) + [0.06; 0.03; 0.03];
            
            for t = 1:T.nT
                [DEL_NN_b_bl, ~, ~] = predictBladeNN(T.NN_predictor_blade, 0, Hs, Tp, wvdir, T_Ueff_bl(t), TI_eff_bl(t), wind_dir_training);
                % res_b_bl = residualBladeNN(NN_blade_const.Value, DEL_NN_b_bl, DD_opt(t), YD_opt(t), 0, 0, Hs, Tp, wvdir, T_Ueff_bl(t), TI_eff_bl(t), turb_opt(t), sin_wd, cos_wd);
                
                [DEL_NN_tw_bl, ~, ~] = predictTowerNN(T.NN_predictor_tower, 0, Hs, Tp, wvdir, T_Ueff_bl(t), TI_eff_bl(t), wind_dir_training);
                res_tw_bl = residualTowerNN(NN_tower_const.Value, DEL_NN_tw_bl, DD_opt(t), YD_opt(t), 0, 0, Hs, Tp, wvdir, T_Ueff_bl(t), TI_eff_bl(t), turb_opt(t), sin_wd, cos_wd);
                
                DEL_blade_bl(t) = DEL_blade_bl(t) + DEL_NN_b_bl / N_batch;
                % DEL_blade_bl(t) = DEL_blade_bl(t) + (DEL_NN_b_bl - res_b_bl) / N_batch;
                DEL_tower_bl(t) = DEL_tower_bl(t) + (DEL_NN_tw_bl - res_tw_bl) / N_batch;
            end



            T_Ueff_bl_hist(j,:) = T_Ueff_bl(:)';

           
            P_bl(j) = sum(0.5 .* paramFLORIS.airDen .* Cp_bl .* A .* T_Ueff_bl.^3);
        end
        
        % BASELINE power with yaw = 0
        E_bl   = mean(P_bl);
        std_bl = std(P_bl);
        % optomized power
        P_opt = mean(power_matrix(1,:));

        DEL_tower_bl_hist_full(end+1,:) = DEL_tower_bl.';  
        DEL_blade_bl_hist_full(end+1,:) = DEL_blade_bl.';

        DEL_blade_bl_sum = sum(DEL_blade_bl);  % max sulle turbine per normalizzazione
        DEL_tower_bl_sum = sum(DEL_tower_bl);

        %% STEP 4: COST FUNCTIONS
        J_avg       = zeros(N_batch, 1);
        J_RA        = zeros(N_batch, 1);
        J_CoL       = zeros(N_batch, 1);
        J_avg_DEL   = zeros(N_batch, 1);
        J_RA_DEL    = zeros(N_batch, 1);
        J_CoL_DEL   = zeros(N_batch, 1);
        
        for i = 1:N_batch
            P_i   = power_matrix(i, :)';
            E_i   = mean(P_i);
            std_i = std(P_i);
            steered = 1:2; % indici dei turbine attivamente yawed
            delta_yaw_max = 20;
            yaw_pen = mean(((gam(i,steered) - gamma_prev(steered)) / delta_yaw_max).^2);
        
            % --- termini DEL normalizzati ---
            dP      = (E_i - E_bl) / E_bl;
            dP_hist(end+1,1) = dP;

            DEL_tower_yawed_sum = sum(DEL_tower_matrix(i,:));
            DEL_blade_yawed_sum = sum(DEL_blade_matrix(i,:));

            dTower  = (DEL_tower_yawed_sum - DEL_tower_bl_sum) / DEL_tower_bl_sum;
            dBlade  = (DEL_blade_yawed_sum - DEL_blade_bl_sum) / DEL_blade_bl_sum;
            DEL_pen = lambda_tower * dTower + lambda_blade * dBlade;
        
            % --- J_avlovely
            % g ---
            J_avg(i)      = E_i;
            J_avg_DEL(i)  = dP - DEL_pen;
        
            % --- J_RA ---
            J_RA(i)       = E_i - lambda * std_i;
            J_RA_DEL(i)   = dP - (lambda * std_i / E_bl) - DEL_pen;
        
            % --- J_CoL ---
            cov_mat  = cov(P_i, P_bl);
            cov_i_bl = cov_mat(1, 2);
            mu_diff  = E_i - E_bl;
            std_diff = sqrt(std_i^2 + std_bl^2 - 2*cov_i_bl);
        
            if std_diff < 1e-6
                J_CoL(i)     = 0.5;
                J_CoL_DEL(i) = 0.5 + DEL_pen;
            else
                J_CoL(i)     = normcdf(0, mu_diff, std_diff);
                % J_CoL_DEL(i) = J_CoL(i) + DEL_pen;
                J_CoL_DEL(i) = J_CoL(i) - lambda_P * dP + DEL_pen + lambda_yaw * yaw_pen;
                
            end
            % --- log ---
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
        
        %% STEP 5: SELEZIONE GAMMA OTTIMALE
        switch opt_type
            case 'avg'
                [~, idx_opt] = max(J_avg);
                cost_name    = 'J_avg (expected)';
            case 'RA'
                [~, idx_opt] = max(J_RA);
                cost_name    = 'J_RA (risk-averse)';
            case 'CoL'
                [~, idx_opt] = min(J_CoL);
                cost_name    = 'J_CoL (conservative)';
            case 'avg_DEL'
                [~, idx_opt] = max(J_avg_DEL);
                cost_name = 'J_avg_DEL';
            case 'RA_DEL'
                [~, idx_opt] = max(J_RA_DEL);
                cost_name = 'J_RA_DEL';
            case 'CoL_DEL'
                [~, idx_opt] = min(J_CoL_DEL);
                cost_name = 'J_CoL_DEL';
            otherwise
                warning('opt_type non valido, uso avg');
                [~, idx_opt] = max(J_avg);
                cost_name    = 'J_avg (expected) [default]';
        end

        gamma_opt = gam(idx_opt, :);  % vettore yaw ottimale [1 x nT]
        yaw_hist(end+1,:) = gamma_opt;

        % yaw_hist(end+1,:) = zeros(1, T.nT); % JUST TO CHECK

        % Per confronto timeseries Baseline/optimized
        DEL_tower_yawed_opt = DEL_tower_matrix(idx_opt,:);
        DEL_blade_yawed_opt = DEL_blade_matrix(idx_opt,:);



        fprintf('\n============POWER COMPARISON BASELINES===============\n')
        fprintf('P opt %.4f MW\n',  P_opt/1e6)
        fprintf('P base %.4f MW\n', E_bl/1e6)

        fprintf('\n=== OTTIMIZZAZIONE YAW [%s] ===\n', cost_name);
        fprintf('  BASELINE power    %.4f MW\n',    E_bl/1e6)
        fprintf('  Gamma ottimale T1: %.2f deg\n',  gamma_opt(1));
        fprintf('  Gamma ottimale T2: %.2f deg\n',  gamma_opt(2));
        fprintf('  E[P] attesa:       %.4f MW\n',   mean(power_matrix(idx_opt,:))/1e6);
        fprintf('  Std[P]:            %.4f MW\n',   std(power_matrix(idx_opt,:))/1e6);
        fprintf('  J_CoL:             %.2f%%\n',    J_CoL(idx_opt)*100);
        fprintf('  J_RA:              %.4f MW\n',   J_RA(idx_opt)/1e6);
        
        fprintf('================================\n\n');


        for i = 1:N_batch
            fprintf(fid_log, ['%.2f,%d,%d,%d,%s,' ...
                '%.4f,%.4f,%.4f,%.4f,%.4f,' ...
                '%.4f,%.4f,%.4f,%.4f,' ...
                '%.6e,%.6e,%.6e,%.6e,%.6f,' ...
                '%.6e,%.6e,%.6e,%.6e,' ...
                '%.6f,%.6f,%.6f,%.6f,' ...
                '%.6e,%.6e,%.6f,%.6e,%.6f,%.6f\n'], ...
                time.current, time.it, i, (i==idx_opt), opt_type, ...
                gam(i,1), gam(i,2), gam(i,3), gamma_prev(1), gamma_prev(2), ...
                ws_hat, wd_hat, ws_std, wd_std, ...
                log_rows(i).E_i, log_rows(i).std_i, E_bl, std_bl, log_rows(i).dP, ...
                log_rows(i).DEL_tower_sum, log_rows(i).DEL_blade_sum, DEL_tower_bl_sum, DEL_blade_bl_sum, ...
                log_rows(i).dTower, log_rows(i).dBlade, log_rows(i).DEL_pen, log_rows(i).yaw_pen, ...
                J_avg(i), J_RA(i), J_CoL(i), J_avg_DEL(i), J_RA_DEL(i), J_CoL_DEL(i));
        end

        gamma_prev = gamma_opt;

        %% STEP 6: APPLICA GAMMA OTTIMALE header
    
        T.States_T(T.StartI,2) = gamma_opt(:); 
         
        % T.States_T(T.StartI,2) = zeros(T.nT,1); % JUST FOR TEST
        
        % Propaga gamma_opt in tutti gli ensemble
        for iE = 1:EnKF.nE
            col_start = (iE-1)*EnKF.nStatesT + 1;
            col_end   = iE*EnKF.nStatesT;
            States_T_iE = EnKF.States_T(:, col_start:col_end);  
            States_T_iE(T.StartI, 2) = gamma_opt(:); % JUST FOR TEST
            % States_T_iE(T.StartI, 2) = zeros(T.nT,1);
            EnKF.States_T(:, col_start:col_end) = States_T_iE;
        end

        fprintf('  [loop time: %.4f s]\n', toc(t_loop));
 


        if T.GeomYaw.active == 1
    
            [DD, YD , turb] = LayoutDisp(T.States_WF(T.StartI(1), 2), T.posBase(:,1), T.posBase(:,2), T.D(1));
    
        else
            
            DD        = 5*ones(T.nT,1);  
            YD        = zeros(T.nT,1);
            turb      = ones(T.nT,1);
        end
   
    
        % Use geometric yaw as farm yaw command [deg]
        if exist('gamma_opt', 'var')
            yaw_deg = gamma_opt(:);
        else
            yaw_deg = geom_yaw(:);  
        end
        % yaw_deg = zeros(T.nT, 1);  % baseline: yaw = 0 % JUST TO CHECK

        % Reshape EnKFM_out in [nT x nSimSteps x 8 x nE]
        Mnum_4d = reshape(EnKFM_out, T.nT, nSimSteps_enkf, 8, EnKF.nE);
        
        % Media sui timestep (dim 2) e sugli ensemble (dim 4)
        mean_wind = squeeze(mean(mean(Mnum_4d(:,:,4,:), 2), 4));  
        TI_wind   = squeeze(mean(mean(Mnum_4d(:,:,3,:), 2), 4))/100 + [0.06, 0.03, 0.03].'; 
    
        % %% DEL EVALUATION:
        % 
        % for i = 1:T.nT
        % 
        %     % -------- Current turbine yaw --------
        %     yaw_i = yaw_deg(i);
        % 
        %     % -------- Upstream turbine info --------
        %     if upstream_idx(i) == 0
        %         gamma_upstream = 0;
        %         yaw0_deg       = 0;
        %     else
        %         gamma_upstream = abs(yaw_deg(upstream_idx(i)));  %da pensarci
        %         yaw0_deg       = abs(yaw_deg(upstream_idx(i)));  %da pensarci
        %     end
        % 
        % 
        %     nds = turb(i);
        % 
        %     %% Tower
        %     % -------- NN prediction --------
        %     [DEL_NN, sigma_NN, logvar_std] = predictTowerNN( ...
        %         T.NN_predictor_tower, yaw_i, Hs, Tp, wvdir, mean_wind(i), TI_wind(i), wind_dir_training);
        % 
        %     % -------- NN residual correction --------
        %     res_hat = residualTowerNN( ...
        %         NN_tower_const.Value, DEL_NN, DD(i), YD(i), gamma_upstream, yaw0_deg, ...
        %         Hs, Tp, wvdir, mean_wind(i), TI_wind(i), nds, sin_wd, cos_wd);
        % 
        %     T.floater.DEL_tower(i) = DEL_NN - res_hat;
        % 
        %     %% BLADE 
        %     % -------- NN prediction --------
        %     [DEL_NN, sigma_NN, logvar_std] = predictBladeNN( ...
        %         T.NN_predictor_blade, yaw_i, Hs, Tp, wvdir, mean_wind(i), TI_wind(i), wind_dir_training);
        % 
        %     % -------- NN residual correction --------
        %     res_hat = residualBladeNN( ...
        %         NN_blade_const.Value, DEL_NN, DD(i), YD(i), gamma_upstream, yaw0_deg, ...
        %         Hs, Tp, wvdir, mean_wind(i), TI_wind(i), nds, sin_wd, cos_wd);
        % 
        %     % -------- Final DEL --------
        % 
        %     T.floater.DEL_blade(i) = DEL_NN - res_hat;
        %     % T.floater.DEL_tower(i) = DEL_NN;   % use this only for NN-only testing
        % 
        % end
    
        
        % DEL_NN_hist_tower(end+1,:) = T.floater.DEL_tower.';
        % DEL_NN_hist_blade(end+1,:) = T.floater.DEL_blade.';

        %% Store DEL history
        DEL_NN_hist_tower(end+1,:) = DEL_tower_yawed_opt;
        DEL_NN_hist_blade(end+1,:) = DEL_blade_yawed_opt; 
        YD_hist(end+1,:) = YD.';
        DD_hist(end+1,:) = DD.';
        TI_hist(end+1,:) = TI_wind.';
        nds_hist(end+1,:) = turb.';
        time_hist_control(end+1,1) = time.current;
    
    
     end


    %% Time
    time.current    = time.current    + time.stepWorld;
    time.current_nO = time.current_nO + time.stepWorld;


end

disp('Simulation done!')

fclose(fid_log);
disp(['Controller log saved to: ' log_path]);

%% Store data
if EnKF.StoreData
    folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName ...
        filesep 'EnKF_Measurements'];
    if ~exist(folder_name, "dir")
        mkdir(folder_name)
    end
    writetable(cat(1,EnKF.M{:}),[folder_name filesep...
         'EnKF_measurements.csv'])
end

%% 
if and(EnKF.StoreData, CLC.Con.StoreData)
    % Generates plots from the outputs in .jpg and vector .pdf format
    postprocess_output([Sim.PathToSim 'Results' filesep CLC.CaseName], T.nT)
end




% Plot EnKF power
EnKF_Vis = plotK_CombinedMeasurement(EnKF_Vis, EnKF, T, Sim, paramFLORIS, 5, 1:3);

% % Plot Flow field
% PlotFlowField(T,[],Wind,Sim,Vis,paramFLORIDyn,paramFLORIS);
% colormap(flipud(mako(100))); clim([2,10])
% % generateWakeCenterVideo(EnKF, T, Sim);

% % Calculate outputs
e = error_EnKF(EnKF,T.nT,thrustFF,uSOWFA,1000);

% Crea cartella per i plot
plot_folder = fullfile(Sim.PathToSim, 'Results', 'CLC_plots');
if ~exist(plot_folder, 'dir')
    mkdir(plot_folder);
end
set(0, 'DefaultFigureVisible', 'on');  % on to visualize plot, off to not do

%% PLOT REAL TIME DATA: 
plot_CLC_results(time_hist, time_hist_control, ws_hist, wd_hist, ws_std_hist, wd_std_hist, ...
                 yaw_hist, DEL_NN_hist_blade, DEL_NN_hist_tower, ...
                 DD_hist, YD_hist, TI_hist, nds_hist, T, Sim, DEL_blade_bl_hist_full, DEL_tower_bl_hist_full);

% Salva tutti i grafici aperti
figList = findall(0, 'Type', 'figure');
for iFig = 1:numel(figList)
    saveas(figList(iFig), fullfile(plot_folder, ['fig_' num2str(figList(iFig).Number) '.png']));
end

figure;
plot(dP_hist,'LineWidth',1.5)
grid on
xlabel('Controller evaluation')
ylabel('\DeltaP = (E_i-E_{bl})/E_{bl}')
title('Relative power gain')

disp('Plots saved to CLC_plots folder');