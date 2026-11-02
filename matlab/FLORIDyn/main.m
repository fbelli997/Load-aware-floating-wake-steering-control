% Copyright (C) <2024>, M Becker
%
% List of the contributors to the development of FLORIDyn: see LICENSE file.
% Description and complete License: see LICENSE file.
	
% This program (FLORIDyn) is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITplatfoHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.

% You should have received a copy of the GNU Affero General Public License
% along with this program (see COPYING file).  If not, see <https://www.gnu.org/licenses/>.
% ======================================================================= %
% Updated: 16. Dez. 2024, M. Becker
% ======================================================================= %

%% MainFLORIDyn Center-Line model
% Improved FLORIDyn approach over the gaussian FLORIDyn model
% Link folder with simulation data
% pathToSimulation = '2021_9T_Data_Fede';
% pathToSimulation = '2026_DEL_test';
pathToSimulation = ['2026_3T_DEL' filesep '3T_IEA15_floating_00'];

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
tic
[T,M,Vis,Mint] = FLORIDynCL_floating(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
t = toc;

disp(['Sec. per sim. step: ' num2str(t/Sim.nSimSteps) ' with '...
    num2str(T.nT) ' turbine(s), total sim. time: ' num2str(t) ' s.'])
%% Plotting & visualization

% Measurements
if Vis.Msmnts.Output
    PlotMeasurements(T,M,Wind,Sim,Vis,paramFLORIDyn,paramFLORIS);
    pause(0.1)
end

% Flow Field
if Vis.FlowField.Plot.Post
    PlotFlowField(T,M,Wind,Sim,Vis,paramFLORIDyn,paramFLORIS)
    pause(0.1)
end

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
