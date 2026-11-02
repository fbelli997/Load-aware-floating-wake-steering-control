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

pathToSimulation = ['2026_6T_Fede_floating' filesep '6T_IEA15_floating_00'];
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

while run_controller
    % Increase counter & check break condition
    disp(['t = ' num2str(time.current,7)])
    time.it = time.it + 1;

    if time.current >= time.end; run_controller = false; end % Simulation end


    %% Estimator
    if mod(time.current_nO, time.stepEnKF) == 0 && time.current > Sim.StartTime
        disp('////////// Run EnKF //////////')
        % ======= Run Ensemble Kalman Filter ======= 
        
        fprintf('[CLC_EnKF_Par START]\n');
        t0 = tic;
        
        CLC_EnKF_Par_floating;
    
        t_EnKF = toc(t0);
        fprintf('[CLC_EnKF_Par END]   (%.6f s)\n', t_EnKF);

        ws_hat = mean(EnKF.States.Vel(T.StartI, :), 2);
        wd_hat = mean(EnKF.States.Dir(T.StartI, :), 2);
        
        % --- EnKF estimated wind at turbine locations (StartI) ---
        fprintf('[EnKF @ t=%.2f]\n', time.current);
        for k = 1:numel(ws_hat)
            fprintf('  T%02d: WS_hat = %.3f m/s, WD_hat = %.2f deg\n', k, ws_hat(k), wd_hat(k));
        end
        
        % ======= Collect "true state" from ensembles
        [T,M, stdM] = EnKF_getMeanState(EnKF, T, paramFLORIDyn);
    
        ws_en = mean(EnKF.States.Vel(T.StartI,:),1);
        ws_hat = mean(ws_en);
        ws_std = std(ws_en);
        wd_en = mean(EnKF.States.Dir(T.StartI,:),1);
        wd_hat = mean(wd_en);
        wd_std = std(wd_en);


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

      % %% Geometric Yaw
      % 
      % if mod(time.curren_nO, time.stepMPC) == 0 && time.current > Sim.StartTime
      % 
      %     disp('Collect ensemble states to run controller')
      % 
      %    if mod(time.current_nO, time.stepEnKF) ~= 0
      %        error("Control not calculated following an EnKF step!")
      %    end
      % 
      %   % ======= Collect "true state" from ensembles
      %   [T,M, stdM] = EnKF_getMeanState(EnKF, T, paramFLORIDyn);
      % 
      %   ws_en = mean(EnKF.States.Vel(T.StartI,:),1);
      %   ws_hat = mean(ws_en);
      %   ws_std = std(ws_en);
      %   wd_en = mean(EnKF.States.Dir(T.StartI,:),1);
      %   wd_hat = mean(wd_en);
      %   wd_std = std(wd_en);
      % 
      %   disp('////////// Run Controller ////////// ')
      %   Wind_Con = predict_wind(EnKF.Wind, time.current, CLC, Sim.TimeStep);
      %   [Con_Test, CLC_Test] = controller(T, Wind_Con, Sim, Con, Vis, paramFLORIDyn,...
      %         paramFLORIS,CLC, time.current);
      % 
      % 
      % end


%     %% MPC
%     if mod(time.current_nO, time.stepMPC) == 0 && time.current > Sim.StartTime
%         % ======= Advance ensembles to current time step
%         disp('Collect ensemble states to run controller')
%         if mod(time.current_nO, time.stepEnKF) ~= 0
%             error("Control not calculated following an EnKF step!")
%         end
% 
% 
%         % ======= Collect "true state" from ensembles
%         [T,M] = EnKF_getMeanState(EnKF, T, paramFLORIDyn);
% 
%         if EnKF.StoreData
%             folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName ...
%                 filesep 'EnKF_States'];
%             if ~exist(folder_name, "dir")
%                 mkdir(folder_name)
%             end
%             writematrix(T.States_WF(:,1)',[folder_name filesep...
%                 'vel'],'WriteMode','append')
%             writematrix(T.States_WF(:,2)',[folder_name filesep...
%                 'dir'],'WriteMode','append')
%             writematrix(T.States_Var_Vel',[folder_name filesep...
%                 'vel_var'],'WriteMode','append')
%             writematrix(T.States_Var_Dir',[folder_name filesep...
%                 'dir_var'],'WriteMode','append')
%         end
% 
% 
%         % DEBUGGING START
%         CLC.Time.StartTime = Sim.StartTime;
%         % DEBUGGING END
% 
%         % ======= MPC ======= 
%         disp('////////// Run Controller ////////// ')
%         % ===== Predict Wind conditions
%         Wind_Con = predict_wind(EnKF.Wind, time.current, CLC, Sim.TimeStep);
%         [Con_Test, CLC_Test] = controller(T, Wind_Con, Sim, Con, Vis, paramFLORIDyn,...
%             paramFLORIS,CLC, time.current);
% 
%         %  ======= Store data ======= 
%         if or(EnKF.StoreData, CLC.Con.StoreData)
%             folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName ...
%                 filesep 'FlowField'];
%             if ~exist(folder_name, "dir")
%                 mkdir(folder_name)
%             end
%             %   Generate flow field and mask
%             fieldLims = Vis.FlowField.Lims;
%             fieldRes  = Vis.FlowField.Res;
%             xAx = fieldLims(1,1):fieldRes:fieldLims(2,1);
%             yAx = fieldLims(1,2):fieldRes:fieldLims(2,2);
%             [X,Y] = meshgrid(xAx,yAx);
% 
%             % Get Mask
%             W = getWeightsT(X(:), Y(:), T,...
%                 Sim.Dyn.Vel, Sim.TimeStep);
%             W = sum(W,2);
%             W = reshape(W,size(X));
% 
%             % Get flow field data
%             U = GetFlowField(T,Wind,Vis,paramFLORIS);
% 
%             writematrix(U,[folder_name filesep...
%                 pad(num2str(time.current),7,"left",'0') '_U.csv'])
%             writematrix(W,[folder_name filesep...
%                 pad(num2str(time.current),7,"left",'0') '_W.csv'])
%             if ~exist([folder_name filesep 'X.csv'],"file")
%                 writematrix(X,[folder_name filesep 'X.csv'])
%                 writematrix(Y,[folder_name filesep 'Y.csv'])
%             end
%             writematrix([time.current T.States_OP(:,1)'],...
%                 [folder_name filesep 'OP_x.csv'],"WriteMode","append")
%             writematrix([time.current T.States_OP(:,2)'],...
%                 [folder_name filesep 'OP_y.csv'],"WriteMode","append")
% 
%             clear fieldLims fieldRes xAx yAx X Y W U
%         end
%     end
% 
% 
    %% Time
    time.current    = time.current    + time.stepWorld;
    time.current_nO = time.current_nO + time.stepWorld;

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
    postprocess_output([Sim.PathToSim 'Results' filesep CLC.CaseName])
end




% Plot EnKF power
EnKF_Vis = plotK_CombinedMeasurement(EnKF_Vis, EnKF, T, Sim, paramFLORIS, 5, 1:6);

% Plot Flow field
PlotFlowField(T,[],Wind,Sim,Vis,paramFLORIDyn,paramFLORIS);
colormap(flipud(mako(100))); clim([2,10])
% generateWakeCenterVideo(EnKF, T, Sim);

% % Calculate outputs
e = error_EnKF(EnKF,T.nT,thrustFF,uSOWFA,1000);

%% ===================== WAKE CENTERLINE PLOT ========================

M_all    = cat(1, EnKF.M{:});
VisTime  = EnKF.M{1}.("Time [s]");
timeVec  = VisTime(1:T.nT:end);
nSteps   = numel(timeVec);

% --- pick the correct column robustly (7 in your table is NOT guaranteed) ---
zVarName = "Wake Centerline Z [m]";
if any(strcmpi(M_all.Properties.VariableNames, zVarName))
    z_all = M_all.(zVarName);
else
    z_all = M_all{:, end};
end

% --- reshape correctly given M_all is [E blocks stacked] ---
nPerE = nSteps * T.nT;
wakeCenterZ = reshape(z_all(1:nPerE * EnKF.nE), nPerE, EnKF.nE);

wakeCenterZ_mean = mean(wakeCenterZ, 2, "omitnan");

colors = lines(T.nT);
nCols  = 2;
nRows  = ceil(T.nT / nCols);

figure('Position', [100 100 900 nRows*250]);
for iT = 1:T.nT
    subplot(nRows, nCols, iT);
    hold on; grid on; box on;

    rowsT = iT : T.nT : nSteps * T.nT;
    zMean = wakeCenterZ_mean(rowsT);
    zStd  = std(wakeCenterZ(rowsT,:), 0, 2, "omitnan");

    fill([timeVec; flipud(timeVec)], ...
         [zMean + zStd; flipud(zMean - zStd)], ...
         colors(iT,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    plot(timeVec, zMean, '-', 'Color', colors(iT,:), 'LineWidth', 1.5);

    xlabel('Time (s)');
    ylabel('Z (m)');
    title(sprintf('Turbine %d', iT-1), 'Color', colors(iT,:));
end

sgtitle('Wake Centerline Z');
folderName = [Sim.PathToSim filesep 'Results' filesep 'EnKF_WakeCenter'];
if ~isfolder(folderName); mkdir(folderName); end
saveas(gcf, [folderName filesep 'WakeCenterlineZ.png']);

%% =================== ROTOR WIND SPEEDS COMPARISON ======================
M_all    = cat(1, EnKF.M{:});
VisTime  = EnKF.M{1}.("Time [s]");
timeVec  = VisTime(1:T.nT:end);
nSteps   = numel(timeVec);

% --------- FLORIDyn effective wind speed (robust pick) ----------
uVarName = "Effective Wind Speed [ms^-1]";
if any(strcmpi(M_all.Properties.VariableNames, uVarName))
    u_all = M_all.(uVarName);
else
    u_all = M_all{:,4};
end

% reshape to [nSteps*T.nT , EnKF.nE] taking only the first nSteps*T.nT rows per ensemble
nPerE      = nSteps * T.nT;
u_mat_allE = reshape(u_all(1:nPerE * EnKF.nE), nPerE, EnKF.nE);

% mean/std across ensembles
u_mean = mean(u_mat_allE, 2, "omitnan");
u_std  = std(u_mat_allE,  0, 2, "omitnan");

% --------- FAST.Farm rotor wind speed from CSV ----------
ffFile = "C:\Users\fbell\OneDrive - Politecnico di Milano\FLORIDyn_Matlab-main_3_floating\FLORIDyn_Matlab-main_3\Simulations\2026_6T_Fede_floating\6T_IEA15_floating_00\Data\FFRotWindSpeeds.csv";

if ~isfile(ffFile)
    error("FFRotWindSpeeds.csv non trovato qui: %s", ffFile);
end

% header starts with # -> use readmatrix to ignore it safely
FF = readmatrix(ffFile);   % columns: Turbine, Time, dt, RotWindSpeed

tFF   = FF(:,2);
idFF  = FF(:,1);           % 0..(T.nT-1)
uFF   = FF(:,4);

% build FF matrix aligned on timeVec (interp per turbina)
uFF_on_timeVec = nan(nSteps, T.nT);
for iT = 1:T.nT
    mask = (idFF == (iT-1));
    if ~any(mask)
        warning("Nessun dato FAST.Farm per turbina %d nel CSV.", iT-1);
        continue
    end
    uFF_on_timeVec(:,iT) = interp1(tFF(mask), uFF(mask), timeVec, "linear", "extrap");
end

% --------- Plot per turbina: FLORIDyn mean +/- std vs FAST.Farm ----------
colors = lines(T.nT);
nCols  = 2;
nRows  = ceil(T.nT / nCols);

figure('Position', [100 100 1000 nRows*270]);
for iT = 1:T.nT
    subplot(nRows, nCols, iT);
    hold on; grid on; box on;

    rowsT   = iT : T.nT : nSteps * T.nT;
    uMeanT  = u_mean(rowsT);
    uStdT   = u_std(rowsT);

    % FLORIDyn band
    fill([timeVec; flipud(timeVec)], ...
         [uMeanT + uStdT; flipud(uMeanT - uStdT)], ...
         colors(iT,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');

    % FLORIDyn mean
    plot(timeVec, uMeanT, '-', 'Color', colors(iT,:), 'LineWidth', 1.5);

    % FAST.Farm rotor wind speed
    plot(timeVec, uFF_on_timeVec(:,iT), '--', 'Color', colors(iT,:), 'LineWidth', 1.2);

    xlabel('Time (s)');
    ylabel('Wind speed (m/s)');
    title(sprintf('Turbine %d', iT-1), 'Color', colors(iT,:));
    legend({'FLORIDyn \pm1\sigma','FLORIDyn mean','FAST.Farm U_{rot}'}, 'Location','best');
end

sgtitle('Effective wind speed comparison: FLORIDyn vs FAST.Farm');

% save
folderName = [Sim.PathToSim filesep 'Results' filesep 'EnKF_RotorWindSpeedCompare'];
if ~isfolder(folderName); mkdir(folderName); end
saveas(gcf, [folderName filesep 'Ueff_vs_FF_Urot.png']);