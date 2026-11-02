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

A = pi*(T.D/2).^2;

%% ===== OPTIMIZATION SETTINGS =====
opt_type = 'avg';  % 'avg', 'RA', 'CoL'
lambda   = 1.0;             
N_batch  = 10;

%% ********************** SPEEDUP **********************
GPR_tower = T.GPR_corrector_tower;
GPR_blade = T.GPR_corrector_blade;
T.GPR_corrector_tower = [];
T.GPR_corrector_blade = [];

GPR_tower_const = parallel.pool.Constant(GPR_tower);
GPR_blade_const = parallel.pool.Constant(GPR_blade);
clear GPR_tower GPR_blade
%% ********************** SPEEDUP **********************

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

        % ********************* TEST HALFED YAW ***********************
        gam(:, 2) = gam(:, 2) * 0.5;  % dimezza yaw turbina 2
        
        % ================RIATTIVA PER VEDERE SET DI YAW==================
        fprintf('Gamma candidati:\n');
        for i = 1:N_batch
            fprintf('  Scenario %02d: gamma = %.2f deg\n', i, gam(i,1));
        end

        %% STEP 2: power matrix (N_batch x N_batch)
        power_matrix = zeros(N_batch, N_batch);
        T_Ueff_hist = zeros(N_batch, N_batch, T.nT);
        turb = [1, 2, 3];

        for i = 1:N_batch       % candidato yaw i
            for j = 1:N_batch   % scenario di vento j

                T_tmp = T;
                T_tmp.States_T(T.StartI, 2) = gam(i,:)';
                T_tmp.States_WF(:, 1) = ws_batch(j);  
                T_tmp.States_WF(:, 2) = wd_batch(j);  
                T_tmp.dep = findTurbineGroups(T_tmp.nT, T_tmp.posBase, T_tmp.States_OP , T_tmp.StartI, T_tmp.States_WF, T_tmp.D, T_tmp.nOP, paramFLORIDyn);
                T_tmp.intOPs = interpolateOPs(T_tmp.nT, T_tmp.StartI, T_tmp.dep, T_tmp.States_OP, T_tmp.posBase, T_tmp.nOP);
                
                [~, Cp,T_Ueff] = RunMultipleFloris( ...
                    T_tmp, paramFLORIS, Wind, turb);
                
                T_Ueff_hist(i,j,:) = T_Ueff;

                yaw_corr  = cosd(gam(i,:)').^(paramFLORIS.p_p);            % [nT x 1] - tutte le turbine del candidato i
                tilt_corr = cosd(T.floater.theta_abs).^(paramFLORIS.p_t);  % [nT x 1] - tilt per ogni turbina
                
                power_matrix(i,j) = sum( ...
                    0.5 .* paramFLORIS.airDen .* Cp .* A .* T_Ueff.^3 ...
                    .* yaw_corr .* tilt_corr);

            end
        end 
    end
    
    %% CONTROLLER:
    % Default yaw command
     if mod(time.current_nO, controller_rate) == 0 && time.current > controller_rate

        disp('////////// Run CONTROLLER and DEL calculator //////////')   
        geom_yaw = zeros(T.nT,1);
        upstream_idx = zeros(T.nT,1);


        %% STEP 3: BASELINE (gamma = 0 per tutti i turbine)
        P_bl = zeros(N_batch, 1);
        T_Ueff_bl_hist = zeros(N_batch, T.nT);
        for j = 1:N_batch
            T_tmp_b = T;
            T_tmp_b.States_T(:, 2) = 0;  % gamma = 0
            T_tmp_b.States_WF(:, 1) = ws_batch(j);
            T_tmp_b.States_WF(:, 2) = wd_batch(j);
            T_tmp_b.dep = findTurbineGroups(T_tmp_b.nT, T_tmp_b.posBase, T_tmp_b.States_OP , T_tmp_b.StartI, T_tmp_b.States_WF, T_tmp_b.D, T_tmp_b.nOP, paramFLORIDyn);
            T_tmp_b.intOPs = interpolateOPs(T_tmp_b.nT, T_tmp_b.StartI, T_tmp_b.dep, T_tmp_b.States_OP, T_tmp_b.posBase, T_tmp_b.nOP);
            
            [~, Cp_bl, T_Ueff_bl] = RunMultipleFloris( ...
                T_tmp_b, paramFLORIS, Wind, turb);

            T_Ueff_bl_hist(j,:) = T_Ueff_bl(:)';

            tilt_corr = cosd(T_tmp_b.floater.theta_abs).^paramFLORIS.p_t;  % ← cosd!
            P_bl(j) = sum(0.5 .* paramFLORIS.airDen .* Cp_bl .* A .* T_Ueff_bl.^3 .* tilt_corr);
        end
        
        % BASELINE power with yaw = 0
        E_bl   = mean(P_bl);
        std_bl = std(P_bl);
        % optomized power
        P_opt = mean(power_matrix(1,:));

        %% STEP 4: COST FUNCTIONS
        J_avg = zeros(N_batch, 1);
        J_RA  = zeros(N_batch, 1);
        J_CoL = zeros(N_batch, 1);

        for i = 1:N_batch
            P_i   = power_matrix(i, :)';
            E_i   = mean(P_i);
            std_i = std(P_i);

            % --- J_avg: avg --- 
            J_avg(i) = E_i;

            % --- J_RA: RA ---
            J_RA(i) = E_i - lambda * std_i;

            % --- J_CoL: CoL ---
            cov_mat  = cov(P_i, P_bl);
            cov_i_bl = cov_mat(1, 2);

            mu_diff  = E_i - E_bl;
            std_diff = sqrt(std_i^2 + std_bl^2 - 2*cov_i_bl);

            if std_diff < 1e-6
                J_CoL(i) = 0.5;
            else
                J_CoL(i) = normcdf(0, mu_diff, std_diff);
            end
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
            otherwise
                warning('opt_type non valido, uso avg');
                [~, idx_opt] = max(J_avg);
                cost_name    = 'J_avg (expected) [default]';
        end

        gamma_opt = gam(idx_opt, :);  % vettore yaw ottimale [1 x nT]
        yaw_hist(end+1,:) = gamma_opt;
        

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

        %% STEP 6: APPLICA GAMMA OTTIMALE
    
        T.States_T(T.StartI,2) = gamma_opt(:);

        % T.States_T(T.StartI,2) = zeros(T.nT,1);

        % Propaga gamma_opt in tutti gli ensemble
        for iE = 1:EnKF.nE
            col_start = (iE-1)*EnKF.nStatesT + 1;
            col_end   = iE*EnKF.nStatesT;
            States_T_iE = EnKF.States_T(:, col_start:col_end);  
            States_T_iE(T.StartI, 2) = gamma_opt(:);
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
    
        % Environmental / operating inputs for DEL models
        Hs        = 4;
        Tp        = 7;
        wvdir      = 0;
    
        % Reshape EnKFM_out in [nT x nSimSteps x 8 x nE]
        Mnum_4d = reshape(EnKFM_out, T.nT, nSimSteps_enkf, 8, EnKF.nE);
        
        % Media sui timestep (dim 2) e sugli ensemble (dim 4)
        mean_wind = squeeze(mean(mean(Mnum_4d(:,:,4,:), 2), 4));  
        TI_wind    = squeeze(mean(mean(Mnum_4d(:,:,3,:), 2), 4))/100 + [0.06, 0.03, 0.03].'; 
    
        %% DEL EVALUATION:
    
        for i = 1:T.nT
    
            % -------- Current turbine yaw --------
            yaw_i = yaw_deg(i);
    
            % -------- Upstream turbine info --------
            if upstream_idx(i) == 0
                gamma_upstream = 0;
                yaw0_deg       = 0;
            else
                gamma_upstream = abs(yaw_deg(upstream_idx(i)));  %da pensarci
                yaw0_deg       = abs(yaw_deg(upstream_idx(i)));  %da pensarci
            end
    
    
            nds = turb(i);
    
            %% Tower
            % -------- NN prediction --------
            [DEL_NN, sigma_NN, logvar_std] = predictTowerNN( ...
                T.NN_predictor_tower, yaw_i, Hs, Tp, wvdir, mean_wind(i), TI_wind(i));
    
            % -------- GPR residual correction --------
            res_hat = residualTowerGPR( ...
                GPR_tower_const.Value, DEL_NN, DD(i), YD(i), gamma_upstream, yaw0_deg, ...
                Hs, Tp, wvdir, mean_wind(i), TI_wind(i), nds);
    
            T.floater.DEL_tower(i) = DEL_NN - res_hat;
    
            %% BLADE 
            % -------- NN prediction --------
            [DEL_NN, sigma_NN, logvar_std] = predictBladeNN( ...
                T.NN_predictor_blade, yaw_i, Hs, Tp, wvdir, mean_wind(i), TI_wind(i));
    
            % -------- GPR residual correction --------
            res_hat = residualBladeGPR( ...
                GPR_tower_const.Value, DEL_NN, DD(i), YD(i), gamma_upstream, yaw0_deg, ...
                Hs, Tp, wvdir, mean_wind(i), TI_wind(i), nds);
    
            % -------- Final DEL --------
    
            T.floater.DEL_blade(i) = DEL_NN - res_hat;
            % T.floater.DEL_tower(i) = DEL_NN;   % use this only for NN-only testing
    
        end
    
        % Store DEL history
        DEL_NN_hist_tower(end+1,:) = T.floater.DEL_tower.';
        DEL_NN_hist_blade(end+1,:) = T.floater.DEL_blade.';
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
                 DD_hist, YD_hist, TI_hist, nds_hist, T, Sim);

% Salva tutti i grafici aperti
figList = findall(0, 'Type', 'figure');
for iFig = 1:numel(figList)
    saveas(figList(iFig), fullfile(plot_folder, ['fig_' num2str(figList(iFig).Number) '.png']));
end
disp('Plots saved to CLC_plots folder');