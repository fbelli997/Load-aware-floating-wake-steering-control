% Copyright (C) <2024>, M BeckerFloridyn_
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

function [T,Mt,Vis,Mint] = FLORIDynCL_floating_ss(T,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS)
%FLORIDYNCL Summary of this function goes here
%   
% ======================================================================= %
% INPUTS
% T     : Turbine data and OP states
% Wind  : Data for the environment (wind speed, direction, ...)
% Sim   : Simulation parameters and settings
% Con   : Settings regarding the control design
% Vis   : Settings regarding the Visualization
% paramFLORIDyn : Parameter for FLORIDyn
% paramFLORIS   : Parameter for FLORIS
% ======================================================================= %
% OUTPUTS
% T := Simulation state (OP states, Turbine states, wind field states(OPs))
% M := Measurements from the simulation (Power, tbd)
% Vis := Visualization data
% Mint := Interaction matrices for the turbines
% ======================================================================= %
%% Run simulation
% T_red, T_addedTI, T_Ueff
M = zeros(Sim.nSimSteps*T.nT,8);
M(:,1) = 1;
M_int = cell(Sim.nSimSteps,1);
yaw_hist = zeros(Sim.nSimSteps, T.nT);
DEL_NN_hist_tower = zeros(Sim.nSimSteps, T.nT);
DEL_NN_hist_blade = zeros(Sim.nSimSteps, T.nT);
YD_hist = zeros(Sim.nSimSteps, T.nT);
DD_hist = zeros(Sim.nSimSteps, T.nT);
TI_hist = zeros(Sim.nSimSteps, T.nT);
nds_hist = zeros(Sim.nSimSteps, T.nT);
turb = ones(Sim.nSimSteps, T.nT);
VCpCt = T.floater.VCpCt.VCpCt;

Vis.StoreFlowField.Vel      = false;
Vis.StoreFlowField.Dir      = false;
Vis.StoreFlowField.Eff      = false;
Vis.StoreFlowField.OP       = false;
Vis.FlowField.Plot.parallel = true;

SimTime = Sim.StartTime;

if paramFLORIDyn.floating == true 
    T.x_float = zeros(Sim.nSimSteps, T.nT);
    T.v_float = zeros(Sim.nSimSteps, T.nT);
end

% per plot
if paramFLORIDyn.floating == true
    T.float_log.x  = zeros(Sim.nSimSteps, T.nT);
    T.float_log.y  = zeros(Sim.nSimSteps, T.nT);
    T.float_log.vx = zeros(Sim.nSimSteps, T.nT);
    T.float_log.vy = zeros(Sim.nSimSteps, T.nT);
    T.float_log.t  = (Sim.StartTime:Sim.TimeStep:Sim.StartTime+(Sim.nSimSteps-1)*Sim.TimeStep).';
end


for it = 1:Sim.nSimSteps



    Sim.SimStep = it;
    % ========== PREDICTION ==========
    % Iterate OPs and states

    % Update position for each starting OP in floating conditions
    T.States_OP(T.StartI,1) = T.posBase(:,1);
    T.States_OP(T.StartI,2) = T.posBase(:,2);


    % ====== MOVED AFTER TO USE Ct instead of axial induction factor ======
    % T = iterateOPs(T,Sim,paramFLORIS,paramFLORIDyn); % can I do that after seUpTmpWFAndRun?????
    % 
    % % ========== Wind Field Pertubation
    % T = pertubationOfTheWF(T,Wind);
    
    % ========== Get FLORIS reductions
    %   ======== Find groups %% for future serial refine (SR) opt important if
    %   SR will be needed (geom yaw + SR (just refine, not SR)

    T.dep = findTurbineGroups(T.nT, T.posBase, T.States_OP , T.StartI, T.States_WF, T.D, T.nOP, paramFLORIDyn);
    
    %   ======= Interpolate OPs
    T.intOPs = interpolateOPs(T.nT, T.StartI, T.dep, T.States_OP, T.posBase, T.nOP);
    
    %   ======= Set up temp Wind farms & run FLORIS
    [tmpM,T,Ct] = setUpTmpWFAndRun(T,paramFLORIS,Wind, turb);
    M((it-1)*T.nT+1:it*T.nT,2:4)    = tmpM;
    M((it-1)*T.nT+1:it*T.nT,1)      = SimTime;
    T.States_T(T.StartI,3)          = tmpM(:,2);
    M_int{it}                       = T.red_arr;

    if paramFLORIDyn.batch == false

        %% GEOMETRIC YAW IMPLEMENTATION, TO CHECK IF HERE OR BEFORE FLORIS CALCULATION:
    
        % ================= GEOMETRIC YAW + DEL PREDICTION =================
    
        % Default yaw command
        geom_yaw = zeros(T.nT,1);
        upstream_idx = zeros(T.nT,1);
    
    
        % ********** TO activate when running the main file and you want to activate yaw
        % angle control and DEL calculation ********************
        if T.GeomYaw.active && it > 1
    
            % Update GeomYaw inputs
            T.GeomYaw.ws_hat = tmpM(:,3);
            T.GeomYaw.wd_hat = Wind.Dir(it, 2); % When running the main file
            % T.GeomYaw.wd_hat = mean(T.States_WF(T.StartI, 2)); % When running the CLC_file
    
    
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
        yaw_hist(it,:) = geom_yaw.';
    
        % Use geometric yaw as farm yaw command [deg]
        yaw_deg = geom_yaw(:);
    
        % Environmental / operating inputs for DEL models
        Hs        = 4;
        Tp        = 7;
        wvdir      = 0;
        mean_wind = tmpM(:,3);
    
        %% when using 6T activate the first one
    
        % TI_wind   = tmpM(:,2) + [0.06, 0.0, 0.0, 0.06, 0.0, 0.0].';
        TI_wind   = tmpM(:,2) + [0.06, 0.03, 0.03].';
        wind_dir = T.States_T(T.StartI, 3);
        sin_wd = sin(deg2rad(wind_dir));
        cos_wd = cos(deg2rad(wind_dir));
    
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
                T.NN_predictor_tower, yaw_i, Hs, Tp, wvdir, mean_wind(i), TI_wind(i), wind_dir(i));
    
            % -------- GPR residual correction --------
            res_hat = residualTowerNN( ...
                T.NN_corrector_tower, DEL_NN, DD(i), YD(i), gamma_upstream, yaw0_deg, ...
                Hs, Tp, wvdir, mean_wind(i), TI_wind(i), nds, sin_wd(i), cos_wd(i));
    
            T.floater.DEL_tower(i) = DEL_NN - res_hat;
    
            %% BLADE 
            % -------- NN prediction --------
            [DEL_NN, sigma_NN, logvar_std] = predictBladeNN( ...
                T.NN_predictor_blade, yaw_i, Hs, Tp, wvdir, mean_wind(i), TI_wind(i), wind_dir(i));
    
            % % -------- GPR residual correction --------
            % res_hat = residualBladeNN( ...
            %     T.NN_corrector_corrector, DEL_NN, DD(i), YD(i), gamma_upstream, yaw0_deg, ...
            %     Hs, Tp, wvdir, mean_wind(i), TI_wind(i), nds, sin_wd(i), cos_wd(i));
    
            % -------- Final DEL --------
    
            % T.floater.DEL_blade(i) = DEL_NN - res_hat;
            T.floater.DEL_blade(i) = DEL_NN;   % use this only for NN-only testing
    
        end
    
        % Store DEL history
        DEL_NN_hist_tower(it,:) = T.floater.DEL_tower.';
        DEL_NN_hist_blade(it,:) = T.floater.DEL_blade.';
        YD_hist(it,:) = YD.';
        DD_hist(it,:) = DD.';
        TI_hist(it,:) = TI_wind.';
        nds_hist(it,:) = turb.';
    
        %% GEOM YAW END
    
        % T.States_T(T.StartI,2) = yaw_deg;
        % T.States_T(T.StartI,2) = zeros(T.nT,1);

    end
   
    T = iterateOPs(T, Sim,paramFLORIS,paramFLORIDyn, Ct); % can I do that after seUpTmpWFAndRun?????
    
    % ========== Wind Field Pertubation
    T = pertubationOfTheWF(T,Wind);

    % ========== Get Wind field variables & correct values
    [T,Wind] = correctVel(T,Wind,SimTime,paramFLORIS,tmpM);
    T = correctDir(T,Wind,SimTime);
    T = correctTi(T,Wind,SimTime);
    
    % Save free wind speed as measurement
    M((it-1)*T.nT+1:it*T.nT,5) = T.States_WF(T.StartI,1);
    
    % ========== Get Control settings
    % Yaw is read in orientation angles
    %T.States_T(T.StartI,1) = getAxInd(Con.YawData,(1:T.nT)',SimTime);     % PLACEHOLDER for axial induction control 
    %%%%%%%% SE RIATTIVI MPC; RIATTIVA GLI YAW QUI SOTTO %%%%%%%
    % T.States_T(T.StartI,2) = (T.States_WF(T.StartI,2) - ...
    %     getYaw(Con.YawData,(1:T.nT)',SimTime)');
    %%%%%%%% NO YAW MISALIGNMENT; SOLO PREDIZIONE WIND SPEED E WIND DIRECTION %%%%%%%

    
    %% FLOATER
    if paramFLORIDyn.floating == true
    %     % DispActiveFloat/Dispfloater
        % [T, x_disp, y_disp, x_dot, y_dot, theta_disp, theta_dot, Faero] = DispActiveFloat(T, tmpM, Sim, paramFLORIDyn, paramFLORIS, it, Ct);
        % tmpM_float = tmpM;
        % 
        % Ue   = tmpM(:,3);
        % Urel = Ue - x_dot;
        % Urel = max(Urel, 0.01);
        % 
        % tmpM_float      = tmpM;
        % tmpM_float(:,3) = Urel;
        % 
        % Ct_corr = CalcCt_floating(Urel, T.floater.VCpCt.VCpCt);
        % 
        % [T, x_disp, y_disp, x_dot, y_dot, theta_disp, theta_dot, Faero] = DispActiveFloat(T, tmpM_float, Sim, paramFLORIDyn, paramFLORIS, it, Ct_corr);

    % % 
    % %% FLOATER
    %     % 
    %     % DispActiveFloat/Dispfloater
    %     % --- Versione semplice (no predictor-corrector) ---
    %     [T, x_disp, y_disp, x_dot, y_dot, theta_disp, theta_dot, Faero] = ...
    %         DispActiveFloat(T, tmpM, Sim, paramFLORIDyn, paramFLORIS, it, Ct);


        % --- Predictor-Corrector ---
        T0 = T;
        Ue = tmpM(:,3);
        [~, ~, ~, x_dot_pred, ~, ~, ~, ~] = ...
            DispActiveFloat(T0, tmpM, Sim, paramFLORIDyn, paramFLORIS, it, Ct);
        
        Urel_pred       = max(Ue - x_dot_pred, 0.01);
        tmpM_float      = tmpM;
        tmpM_float(:,3) = Urel_pred;
        Ct_corr         = CalcCt_floating(Urel_pred, VCpCt);
        
        [T, x_disp, y_disp, x_dot, y_dot, theta_disp, theta_dot, Faero] = ...
            DispActiveFloat(T0, tmpM_float, Sim, paramFLORIDyn, paramFLORIS, it, Ct_corr);

    end
    
    T.floater.theta_abs = T.floater.theta_undisturbed + T.floater.theta(:,1);
    abs_ang_deg = rad2deg(T.floater.theta_abs);

    % per plot
    T.float_log.x(it,:)  = T.floater.x(:,1).';
    T.float_log.y(it,:)  = T.floater.y(:,1).';
    T.float_log.vx(it,:) = T.floater.vx(:,1).';
    T.float_log.vy(it,:) = T.floater.vy(:,1).';
    T.float_log.tilt(it,:) = rad2deg(T.floater.theta_undisturbed + T.floater.theta(:,1)).';
    T.float_log.theta_dot(it,:) = T.floater.theta_dot(:,1).';
    T.float_log.pitch_ws (it,:) = T.floater.pitch_ws;
    T.float_log.abs_disp(it,:) = x_disp.';
    T.float_log.abs_ws(it,:) = x_dot.';  % BOTH TILT AND SURGE CONTRIBUTES

    % ========== Calculate Power
    % P = getPower(T,tmpM,paramFLORIS,Con);
    Thrust = Faero(:,1); 

    % M((it-1)*T.nT+1:it*T.nT,6) = P;
    M((it-1)*T.nT+1:it*T.nT,6) = Thrust; % For using thrust instead of Power
    
    % for iT = 1:T.nT
    % 
    %     rowM = (it-1)*T.nT + iT;
    %     idx_half_way = T.StartI(iT)-40;
    %     if isempty(T.dep{iT})
    %         idxOP = max(idx_half_way, 1);
    %     else
    %         idxOP = max(idx_half_way - 1, 1);
    %     end
    %     M(rowM, 7) = T.States_OP(idxOP, 3);
    %     M(rowM, 8) = T.States_OP(idxOP, 5);
    % end
    for iT = 1:T.nT
    
        rowM = (it-1)*T.nT + iT;
    
        % --- Save incoming wake centerline deflection at this turbine ---
        % Convention:
        %   - For turbine 1: no upstream -> 0
        %   - For turbine iT>1: use upstream turbine (iT-1) wake centerline
        %     evaluated at the downwind distance to turbine iT
    
        if iT == 1
            % No upstream wake
            M(rowM, 7) = 0;   % Wake centerline Y [m]  (deltaY)
            M(rowM, 8) = 0;   % Wake centerline Z [m]  (deltaZ)
            continue
        end
    
        iUp = iT - 1;  % upstream turbine index
    
        % Wake-frame angle for upstream turbine (using wind direction at StartI)
        phi = angSOWFA2world(T.States_WF(T.StartI(iUp), 2));
    
        % Vector from upstream hub to current hub (world coords)
        ri = T.posBase(iUp,:) + T.posNac(iUp,:);
        rj = T.posBase(iT,:)  + T.posNac(iT,:);
        dr = (rj - ri).';
    
        % Rotate into upstream wake frame to get downwind distance
        R  = [cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
        dr_w = R * dr;
        xdw  = dr_w(1);  % downwind distance [m]
    
        % Pick OP from upstream turbine block closest to that xdw
        idxUp = T.StartI(iUp) + (0:T.nOP-1);
        [~, krel] = min(abs(T.States_OP(idxUp, 4) - xdw));
        idxOP = idxUp(krel);
    
        % Store deflection components (wake-frame centerline)
        M(rowM, 7) = T.States_OP(idxOP, 5);  % Wake centerline Y [m]  (deltaY)
        M(rowM, 8) = T.States_OP(idxOP, 6);  % Wake centerline Z [m]  (deltaZ
    end

    % ========== Live Plotting
    if Vis.FlowField.Plot.Online
        Vis = PlotFlowField(T,M,Wind,Sim,Vis,paramFLORIDyn,paramFLORIS,...
            SimTime);
    end
    if Vis.FlowField.Error.Online
        if sum(Vis.FlowField.Error.Steps==SimTime)==1
            errorVTK([Vis.FlowField.Error.ValidationPath ...
                num2str(SimTime) '\' Vis.FlowField.Error.SnapName],...
                [num2str(SimTime) 'FlowError'],...
                T,Wind,Sim,Vis,paramFLORIDyn,paramFLORIS)
        end
    end
    if Vis.StoreFlowField.Vel
        createAndStoreWeightedFlowField(T,Wind,Sim,Vis,paramFLORIS,...
            SimTime,'Vel')
    end
    if Vis.StoreFlowField.Dir
        createAndStoreWeightedFlowField(T,Wind,Sim,Vis,paramFLORIS,...
            SimTime,'Dir')
    end
    if Vis.StoreFlowField.Eff
        createAndStoreWeightedFlowField(T,Wind,Sim,Vis,paramFLORIS,...
            SimTime,'Eff')
    end
    if Vis.StoreFlowField.OP
        createAndStoreWeightedFlowField(T,Wind,Sim,Vis,paramFLORIS,...
            SimTime,'OP')
    end
    % ========== Increase time
    SimTime = SimTime + Sim.TimeStep;

    % ==== Store the centerline data ====
    if Vis.CL.Store
        storeCLInfo(T,Sim,it);
    end
end
    
% if paramFLORIDyn.floating == true
%     t  = T.float_log.t;
%     nT = T.nT;
%     nr = 3; nc = 2;
% 
%     %% Surge displacement (x)
%     figure('Name','Floater surge displacement (x)');
%     for iT = 1:nT
%         subplot(nr,nc,iT);
%         set(gca, 'FontSize', 20);
%         plot(t, T.float_log.x(:,iT)); grid on;
%         xlabel('Time [s]'); ylabel('x disp [m]');
%         title(sprintf('Turbine %d', iT));
%     end
%     sgtitle('Floater surge displacement (x)', 'FontSize', 24);
% 
%     % %% Sway displacement (y)
%     % figure('Name','Floater sway displacement (y)');
%     % for iT = 1:nT
%     %     subplot(nr,nc,iT);
%     %     set(gca, 'FontSize', 20);
%     %     plot(t, T.float_log.y(:,iT)); grid on;
%     %     xlabel('Time [s]'); ylabel('y disp [m]');
%     %     title(sprintf('Turbine %d', iT));
%     % end
%     % sgtitle('Floater sway displacement (y)', 'FontSize', 24);
% 
%     %% Surge velocity (vx)
%     figure('Name','Floater surge velocity (vx)');
%     for iT = 1:nT
%         subplot(nr,nc,iT);
%         set(gca, 'FontSize', 20);
%         plot(t, T.float_log.vx(:,iT)); grid on;
%         xlabel('Time [s]'); ylabel('vx [m/s]');
%         title(sprintf('Turbine %d', iT));
%     end
%     sgtitle('Floater surge velocity (vx)', 'FontSize', 24);
% 
%     % %% Sway velocity (vy)
%     % figure('Name','Floater sway velocity (vy)');
%     % for iT = 1:nT
%     %     subplot(nr,nc,iT);
%     %     set(gca, 'FontSize', 20);
%     %     plot(t, T.float_log.vy(:,iT)); grid on;
%     %     xlabel('Time [s]'); ylabel('vy [m/s]');
%     %     title(sprintf('Turbine %d', iT));
%     % end
%     % sgtitle('Floater sway velocity (vy)', 'FontSize', 24);
% 
%     %% Pitch (tilt)
%     figure('Name','Platform pitch (deg)');
%     for iT = 1:nT
%         subplot(nr,nc,iT);
%         set(gca, 'FontSize', 24);
%         plot(t, T.float_log.tilt(:,iT)); grid on;
%         xlabel('Time [s]'); ylabel('\theta [deg]');
%         title(sprintf('Turbine %d', iT));
%     end
%     sgtitle('Platform pitch (deg)', 'FontSize', 24);
% 
%     %% abs_disp
%     figure('Name','Platform x motion due to pitch and surge (m)');
%     for iT = 1:nT
%         subplot(nr,nc,iT);
%         set(gca, 'FontSize', 20);
%         plot(t, T.float_log.abs_disp(:,iT)); grid on;
%         xlabel('Time [s]'); ylabel('abs disp x [m]');
%         title(sprintf('Turbine %d', iT));
%     end
%     sgtitle('Platform x motion due to pitch and surge (m)', 'FontSize', 24);
% 
%     %% abs_ws
% 
%     figure('Name','Absolute wind speed seen by each turbine');
%     for iT = 1:nT
%         subplot(nr,nc,iT);
%         set(gca, 'FontSize', 20);
%         plot(t, T.float_log.abs_ws(:,iT)+10.59); grid on;
%         xlabel('Time [s]'); ylabel('Wind Speed [m/s]');
%         title(sprintf('Turbine %d', iT));
%     end
%     sgtitle('Absolute wind speed seen by each turbine', 'FontSize', 24);
% 
%     %% vx vs pitch_ws
%     figure('Name','vx and pitch-induced velocities');
%     for iT = 1:nT
%         subplot(nr,nc,iT);
%         set(gca, 'FontSize', 20);
%         hold on; grid on;
%         plot(t, T.float_log.vx(:,iT));
%         plot(t, T.float_log.pitch_ws(:,iT), 'LineWidth', 2);
%         xlabel('Time [s]'); ylabel('Velocity [m/s]');
%         title(sprintf('Turbine %d', iT));
%         legend({'vx','pitch\_ws'}, 'Location','best');
%     end
%     sgtitle('Floater surge velocity and Floater pitch velocities', 'FontSize', 24);
% 
%     % %% Trajectory (x-y)
%     % figure('Name','Floater trajectory (x-y)');
%     % for iT = 1:nT
%     %     subplot(nr,nc,iT);
%     %     set(gca, 'FontSize', 20);
%     %     plot(T.float_log.x(:,iT), T.float_log.y(:,iT)); grid on;
%     %     xlabel('x [m]'); ylabel('y [m]');
%     %     title(sprintf('Turbine %d', iT), 'FontSize', 24);
%     %     axis equal;
%     % end
%     % sgtitle('Floater trajectory (x-y)', 'FontSize', 24);
% end



if Vis.FlowField.Plot.Online
    Vis.Film.InProgress = false;
    % Write saved frames as a movie to the simulation folder
    v = VideoWriter(Vis.Film.MovFileEffU);
    v.FrameRate = 12;
    v.Quality   = 100;
    open(v)
    writeVideo(v,Vis.Film.FrmFileEffU)
    close(v)
end

%% Plots for DEL evaluation when using main.
if paramFLORIDyn.CLC == false

    t = (Sim.StartTime:Sim.TimeStep:Sim.StartTime+(Sim.nSimSteps-1)*Sim.TimeStep).';
    
    
    figure('Name','Yaw history');
    nCols = 2;
    nRows = ceil(T.nT/nCols);
    
    for iT = 1:T.nT
        subplot(nRows,nCols,iT)
        plot(t, yaw_hist(:,iT), 'LineWidth', 1.5)
        grid on
        xlabel('Time [s]')
        ylabel('Yaw [deg]')
        title(sprintf('Turbine %d', iT))
        set(gca, 'Fontsize', 20)
    end
    sgtitle('Geometric yaw history')
    
    figure('Name', 'DEL_ History blade')
    
    for iT = 1:T.nT
        subplot(nRows,nCols,iT)
        plot(t, DEL_NN_hist_blade(:,iT), 'LineWidth', 1.5)
        grid on
        xlabel('Time [s]')
        ylabel('DEL Blade [deg]')
        title(sprintf('Turbine %d', iT))
        set(gca, 'Fontsize', 20)
    end
    sgtitle('DEL blade history')
    
    figure('Name', 'DEL_Tower History')
    
    for iT = 1:T.nT
        subplot(nRows,nCols,iT)
        plot(t, DEL_NN_hist_tower(:,iT), 'LineWidth', 1.5)
        grid on
        xlabel('Time [s]')
        ylabel('DEL Tower [deg]')
        title(sprintf('Turbine %d', iT))
        set(gca, 'Fontsize', 20)
    end
    sgtitle('DEL tower history')
    
    figure('Name', 'DD History')
    
    for iT = 1:T.nT
        subplot(nRows,nCols,iT)
        plot(t, DD_hist(:,iT), 'LineWidth', 1.5)
        grid on
        xlabel('Time [s]')
        ylabel('DD [deg]')
        title(sprintf('Turbine %d', iT))
        set(gca, 'Fontsize', 20)
    end
    sgtitle('DD history')
    
    figure('Name', 'YD History')
    
    for iT = 1:T.nT
        subplot(nRows,nCols,iT)
        plot(t, YD_hist(:,iT), 'LineWidth', 1.5)
        grid on
        xlabel('Time [s]')
        ylabel('YD [deg]')
        title(sprintf('Turbine %d', iT))
        set(gca, 'Fontsize', 20)
    end
    sgtitle('YD history')
    
    figure('Name', 'TI History')
    
    for iT = 1:T.nT
        subplot(nRows,nCols,iT)
        plot(t, TI_hist(:,iT), 'LineWidth', 1.5)
        grid on
        xlabel('Time [s]')
        ylabel('TI [deg]')
        title(sprintf('Turbine %d', iT))
        set(gca, 'Fontsize', 20)
    end
    sgtitle('TI history')
    
    figure('Name', 'nds')
    
    for iT = 1:T.nT
        subplot(nRows,nCols,iT)
        plot(t, nds_hist(:,iT), 'LineWidth', 1.5)
        grid on
        xlabel('Time [s]')
        ylabel('nds [deg]')
        title(sprintf('Turbine %d', iT))
        set(gca, 'Fontsize', 20)
    end
    sgtitle('nds history')
end

% Convert M to table and scale measurements accordingly
Mt = array2table(M*diag([1;100;100;1;1;10^(-6);1;1]),...
    'VariableNames',{'Time [s]','Foreign Reduction [%]','Added turbulence [%]',...
    'Effective Wind Speed [ms^-1]','Free Wind Speed [ms^-1]', ...
    'Thrust generated [MN]', 'Wake centerline Y [m]', 'Wake centerline Z [m]'});
Mint = [Mt.("Time [s]") cell2mat(M_int)];
end
