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

    %% CONTROLLER:
    % Default yaw command
     if mod(time.current_nO, controller_rate) == 0 && time.current > controller_rate

        disp('////////// Run CONTROLLER and DEL calculator //////////')   
        geom_yaw = zeros(T.nT,1);
        upstream_idx = zeros(T.nT,1);
    
        if T.GeomYaw.active == 1
    
    
            %% Get geometric yaw and upstream relationship
            [geom_yaw, ~, upstream_idx] = GeometricYaw( ...
                T.GeomYaw, T.posBase(:,1), T.posBase(:,2));
    
            [DD, YD , turb] = LayoutDisp(T.States_WF(T.StartI(1), 2), T.posBase(:,1), T.posBase(:,2), T.D(1));
    
        else
        DD        = 5*ones(T.nT,1);  
        YD        = zeros(T.nT,1);
        turb      = ones(T.nT,1);
        end
    
        % Store yaw history
        yaw_hist(end+1,:) = geom_yaw.';
    
        % Use geometric yaw as farm yaw command [deg]
        yaw_deg = geom_yaw(:);
        
        % Environmental / operating inputs for DEL models
        Hs        = 4;
        Tp        = 7;
        wvdir      = 0;
    
        mean_wind = mean(reshape(Mnum(end-8:end, 4), 3, 3), 2);
    
        %% when using 6T activate the first one
    
        % TI_wind   = tmpM(:,2) + [0.06, 0.0, 0.0, 0.06, 0.0, 0.0].';
        TI_wind   = mean(reshape(Mnum(end-8:end, 3), 3, 3), 2)/100 + [0.06, 0.03, 0.03].';
    
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
                T.GPR_corrector_tower, DEL_NN, DD(i), YD(i), gamma_upstream, yaw0_deg, ...
                Hs, Tp, wvdir, mean_wind(i), TI_wind(i), nds);
    
            T.floater.DEL_tower(i) = DEL_NN - res_hat;
    
            %% BLADE 
            % -------- NN prediction --------
            [DEL_NN, sigma_NN, logvar_std] = predictBladeNN( ...
                T.NN_predictor_blade, yaw_i, Hs, Tp, wvdir, mean_wind(i), TI_wind(i));
    
            % -------- GPR residual correction --------
            res_hat = residualBladeGPR( ...
                T.GPR_corrector_blade, DEL_NN, DD(i), YD(i), gamma_upstream, yaw0_deg, ...
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

        %% GEOM YAW END
    
        T.States_T(T.StartI,2) = yaw_deg;
        % T.States_T(T.StartI,2) = zeros(T.nT,1);

     end


    %% Time
    time.current    = time.current    + time.stepWorld;
    time.current_nO = time.current_nO + time.stepWorld;
    % fprintf('  [loop time: %.4f s]\n', toc(t_loop));   % ← aggiungi qui
% 
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

%% PLOT REAL TIME DATA: 
plot_CLC_results(time_hist, time_hist_control, ws_hist, wd_hist, ws_std_hist, wd_std_hist, ...
                 yaw_hist, DEL_NN_hist_blade, DEL_NN_hist_tower, ...
                 DD_hist, YD_hist, TI_hist, nds_hist, T, Sim);