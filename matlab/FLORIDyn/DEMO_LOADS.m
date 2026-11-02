%% ===================== DEMO LOAD FLORIDYN ===========================
%% MainFLORIDyn Center-Line model
% Improved FLORIDyn approach over the gaussian FLORIDyn model
% Link folder with simulation data
% pathToSimulation = '2021_9T_Data_Fede';
% pathToSimulation = '2026_DEL_test';
pathToSimulation = ['2026_3T_DEL_DEMO' filesep '3T_IEA15_floating_00'];

%% Load data from the simulation
% Reset the Matlab Path and load essential paths & the simulation path
addPaths;
% Get the settings for the wind field, visuali2026_DEL_testzation, controller and Sim.
[Wind, Vis, Sim, Con] = setup();
% Add according functions to the search path
addFLORISPaths;
addFLORIDynPaths;
% Load linked dataclose 
turbProp        = turbineArrayPropertiesIEA15();
paramFLORIS     = parameterFLORIS();
paramFLORIDyn   = parameterFLORIDyn();
%to distinguish the main file which doesn't use the ensembles and using the
%same "prepareSimulationFloating_main.m" file.
paramFLORIDyn.CLC = false;

%% Preprocess loaded data
[T, Wind, Sim, Con, paramFLORIS, Vis] = ...
    prepareSimulationFloating_main(Wind,Con,paramFLORIDyn,paramFLORIS,turbProp,Sim,Vis);
clear turbProp
paramFLORIDyn.batch = false;

%% ====== Init simulation
% Run initial conditions until no more change happens
T = initSimulation(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);

%% ============ RUN SIMULATION ============
% T := Simulation state (OP states, Turbine states, wind field states(OPs))
% M := Measurements from the simulation (T_red, T_addedTI, T_Ueff, T_pow)

yaw_vec = -20:2:0;
nY      = length(yaw_vec);

transient_step = 250;  % step da cui iniziare la media (butta il transitorio)

% pre-allocazione (scalari a livello farm, uno per ciascun yaw)
DELs_TOWER = zeros(nY, 1);
DELs_BLADE = zeros(nY, 1);
POWERs     = zeros(nY, 1);
yaws       = zeros(nY, 1);

T_init = T;  % snapshot per ripartire dalle stesse condizioni

for i = 1:nY
    yaw   = [yaw_vec(i), yaw_vec(i), 0];
    T_run = T_init;
    
    tic
    [T_run, M, Vis, Mint, DEL_NN_hist_tower, DEL_NN_hist_blade, P] = ...
        FLORIDynDL_floating(T_run, Wind, Sim, Con, Vis, paramFLORIDyn, paramFLORIS, yaw);
    t = toc;
        
    % taglia il transitorio -> matrice (nSimSteps-250) x 3
    DEL_tower_cut = DEL_NN_hist_tower(transient_step+1:end, :);
    DEL_blade_cut = DEL_NN_hist_blade(transient_step+1:end, :);
    P_cut         = P(transient_step+1:end, :);
    
    % media nel tempo, per turbina -> vettore 1 x 3
    DEL_tower_mean_t = mean(DEL_tower_cut, 1);
    DEL_blade_mean_t = mean(DEL_blade_cut, 1);
    P_mean_t         = mean(P_cut, 1);
    
    % somma sulle 3 turbine -> scalare
    DELs_TOWER(i) = sum(DEL_tower_mean_t);
    DELs_BLADE(i) = sum(DEL_blade_mean_t);
    POWERs(i)     = sum(P_mean_t);
    yaws(i)       = yaw_vec(i);

    disp(['Sec. per sim. step: ' num2str(t/Sim.nSimSteps) ' with '...
        num2str(T.nT) ' turbine(s), total sim. time: ' num2str(t) ' s.'])
end

%% Plotting & visualization
%% ============ RIATIVARE SE VUOI VEDERE PLOT POTENZA ====================
% % Measurements
% if Vis.Msmnts.Output
%     PlotMeasurements(T,M,Wind,Sim,Vis,paramFLORIDyn,paramFLORIS);
%     pause(0.1)
% end



%% ================ PARETO PLOTS ====================
% Normalize DEL/Power over the baseline (yaw = 0) values
idx0 = find(yaw_vec == 0, 1);
if isempty(idx0)
    error('DEMO_LOADS:NoBaseline', 'yaw_vec must contain 0 deg to use as baseline.');
end

DELs_TOWER_norm = DELs_TOWER / DELs_TOWER(idx0);
DELs_BLADE_norm = DELs_BLADE / DELs_BLADE(idx0);
POWERs_norm     = POWERs     / POWERs(idx0);

figPareto = figure('Name', 'Pareto DEL vs Power (normalized to baseline)');

% ---- Left: Tower ----
ax1 = subplot(1, 2, 1);
plot(DELs_TOWER_norm, POWERs_norm, '-k', 'LineWidth', 1.5); hold on;
scatter(DELs_TOWER_norm, POWERs_norm, 100, yaws, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
xlabel('Normalized Farm DEL Tower', 'FontSize', 18);
ylabel('Normalized Farm power',     'FontSize', 18);
grid on;
set(ax1, 'FontSize', 16, 'LineWidth', 1.2);
title(ax1, '(a)', 'FontSize', 18, 'FontWeight', 'normal');

% ---- Right: Blade ----
ax2 = subplot(1, 2, 2);
plot(DELs_BLADE_norm, POWERs_norm, '-k', 'LineWidth', 1.5); hold on;
scatter(DELs_BLADE_norm, POWERs_norm, 100, yaws, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
xlabel('Normalized Farm DEL Blade', 'FontSize', 18);
ylabel('Normalized Farm power',     'FontSize', 18);
grid on;
set(ax2, 'FontSize', 16, 'LineWidth', 1.2);
title(ax2, '(b)', 'FontSize', 18, 'FontWeight', 'normal');

% ---- Shared colormap and single colorbar ----
colormap(parula);
cb = colorbar(ax2, 'Position', [0.92 0.15 0.02 0.75]);
cb.Label.String = 'Yaw T1=T2 [deg]';
cb.Label.FontSize = 18;
cb.FontSize = 16;

% Rendi lo spazio per la colorbar accorciando i subplot
% (ax2 e' piu' stretto per lasciare piu' margine rispetto alla colorbar)
ax1.Position(3) = 0.35;
ax2.Position(3) = 0.325;

% Opzionale: dimensiona la figura per un aspect ratio da paper
set(gcf, 'Position', [100 100 1400 550]);

% ---- Save figure as .fig and .png ----
resultsFolder = [Sim.PathToSim 'Results'];
if ~exist(resultsFolder, 'dir')
    mkdir(resultsFolder);
end
saveas(figPareto, fullfile(resultsFolder, 'Pareto_DEL_vs_Power_normalized.fig'));
saveas(figPareto, fullfile(resultsFolder, 'Pareto_DEL_vs_Power_normalized.png'));

% 3D Flow Field (.vtk for ParaView)
if Vis.FlowField3D.Generate
    gen3DFlowField(T,M,Wind,Sim,Vis,paramFLORIDyn,paramFLORIS);
end

% Final state to skip initialization
if Sim.SaveFinalState
    save([Sim.PathToSim 'T_final.mat'],'T')
end

% Save movie if online plotting was activated
if Vis.FlowField.Plot.Online
    Vis.Film.InProgress = false;
    % Write saved frames as a movie to the simulation folder
    v = VideoWriter(Vis.Film.MovFileEffU);
    v.FrameRate = 20;
    v.Quality   = 100;
    open(v)
    writeVideo(v,Vis.Film.FrmFileEffU)
    close(v)
end