% Copyright (C) <2024>, M Becker
%
% List of the contributors to the development of FLORIDyn: see LICENSE file.
% Description and complete License: see LICENSE file.
	T;
% This program (FLORIDyn) is free software: you can rFedistribute it and/or modify
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

% CLC EnKF script
tic
% Set inputs
Sim.StartTime   = time.current - time.stepEnKF;
Sim.EndTime     = time.current;
Sim.nSimSteps   = EnKF.nS;

% Prepare parallel computing
EnKFStatesVel = EnKF.States.Vel;
EnKFStatesDir = EnKF.States.Dir;
EnKFStates_OP = reshape(EnKF.States_OP,[],EnKF.nE);
EnKFStates_T  = reshape(EnKF.States_T,[],EnKF.nE);
EnKFnStatesOP = EnKF.nStatesOP;
EnKFnStatesT  = EnKF.nStatesT;
EnKFOutputThrust = EnKF.Output.Thrust;
EnKFM         = zeros(Sim.nSimSteps * T.nT * 8,EnKF.nE);



% ===================================================

% si parfor

% parfor iE = 1:EnKF.nE
%     % ===== Assign relevant ensemble states
%     T_tmp = T;
% 
%     % Ensemble States
%     T_tmp.States_WF(:,1) = EnKFStatesVel(:,iE);
%     T_tmp.States_WF(:,2) = EnKFStatesDir(:,iE);
% 
%     % OP States
%     T_tmp.States_OP = reshape(EnKFStates_OP(:,iE),[],EnKFnStatesOP);
%         %EnKF.States_OP(:,EnKF.nStatesOP*(iE-1)+1:EnKF.nStatesOP*iE);
% 
%     % T States
%     T_tmp.States_T = reshape(EnKFStates_T(:,iE),[],EnKFnStatesT);
%         %EnKF.States_T(:,EnKF.nStatesT*(iE-1)+1:EnKF.nStatesT*iE);
% 
%     %T_tmp = EnKF_AssignEnStates(EnKF,T,iE);
% 
%     % =========== Run simulation ================
%     [T_tmp,M,~] = ...
%         FLORIDynCL_floating(T_tmp,Wind,Sim,Con,Vis,paramFLORIDyn,paramFLORIS);
% 
%     % ===== Stores relevant measurements and states
%     % Store Ensemble states
%     EnKFStatesVel(:,iE) = T_tmp.States_WF(:,1);
%     EnKFStatesDir(:,iE) = T_tmp.States_WF(:,2);
% 
%     % Store OP states
%     EnKFStates_OP(:,iE) = T_tmp.States_OP(:);
% 
%     % Store T states
%     EnKFStates_T(:,iE)  = T_tmp.States_T(:);
% 
%     % M è una table, prendo solo le 6 colonne numeriche già in ordine:
%     Mnum = [ ...
%         M.("Time [s]"), ...
%         M.("Foreign Reduction [%]"), ...
%         M.("Added turbulence [%]"), ...
%         M.("Effective Wind Speed [ms^-1]"), ...
%         M.("Free Wind Speed [ms^-1]"), ...
%         M.("Thrust generated [MN]") ...
%     ];
% 
%     EnKFM(:,iE) = Mnum(:);
% 
%     % potenza ultimo step (ultime nT righe, colonna 6)
%     % EnKFOutputPow(:,iE) = Mnum(end-T.nT+1:end,6);
%     EnKFOutputThrust(:,iE) = Mnum(end-T.nT+1:end,6);
% 
% 
%     % % Store measurements
%     % EnKFM(:,iE) = reshape(table2array(M),[],1);
%     % %EnKF.M{iE} = [EnKF.M{iE}; M];
%     % EnKFOutputThrust(:,iE) = ...
%     %     table2array(M(end-T.nT+1:end,6));
%     % %EnKF = EnKF_StoreEnStates(EnKF,M, T_tmp, iE);
% end





t_FLORIDyn = zeros(EnKF.nE,1);

%% ********************** SPEEDUP **********************
nSimSteps_enkf = Sim.nSimSteps;

%% ********************** SPEEDUP **********************

% Preallocate outputs per parfor
EnKFStatesVel_out    = zeros(size(EnKFStatesVel));
EnKFStatesDir_out    = zeros(size(EnKFStatesDir));
EnKFStates_OP_out    = zeros(size(EnKFStates_OP));
EnKFStates_T_out     = zeros(size(EnKFStates_T));
EnKFOutputThrust_out = zeros(size(EnKFOutputThrust));
EnKFM_out            = zeros(size(EnKFM));
FloaterStates_out    = zeros(T.nT, 7, EnKF.nE);

parfor iE = 1:EnKF.nE
    T_tmp = T;
    T_tmp.States_WF(:,1) = EnKFStatesVel(:,iE);
    T_tmp.States_WF(:,2) = EnKFStatesDir(:,iE);
    T_tmp.NN_corrector_tower = NN_tower_const.Value;
    T_tmp.NN_corrector_blade = NN_blade_const.Value;
    T_tmp.States_OP = reshape(EnKFStates_OP(:,iE),[],EnKFnStatesOP);
    T_tmp.States_T  = reshape(EnKFStates_T(:,iE),[],EnKFnStatesT);

    Float_states = T.floater.T_EnKF_States(:,:,iE);
    T_tmp.floater.x         = Float_states(:,1);
    T_tmp.floater.y         = Float_states(:,2);
    T_tmp.floater.vx        = Float_states(:,3);
    T_tmp.floater.vy        = Float_states(:,4);
    T_tmp.floater.theta     = Float_states(:,5);
    T_tmp.floater.theta_dot = Float_states(:,6);
    T_tmp.floater.theta_abs = Float_states(:,7);

    Sim_tmp = Sim;
    Sim_tmp.nSimSteps = nSimSteps_enkf;

    [T_tmp, M, ~] = FLORIDynCL_floating(T_tmp, Wind, Sim_tmp, Con, Vis, paramFLORIDyn, paramFLORIS);

    % Floater states output
    Float_out_ie = zeros(T.nT, 7);
    Float_out_ie(:,1) = T_tmp.floater.x;
    Float_out_ie(:,2) = T_tmp.floater.y;
    Float_out_ie(:,3) = T_tmp.floater.vx;
    Float_out_ie(:,4) = T_tmp.floater.vy;
    Float_out_ie(:,5) = T_tmp.floater.theta;
    Float_out_ie(:,6) = T_tmp.floater.theta_dot;
    Float_out_ie(:,7) = T_tmp.floater.theta_abs;
    FloaterStates_out(:,:,iE) = Float_out_ie;

    % Ensemble states output
    EnKFStatesVel_out(:,iE) = T_tmp.States_WF(:,1);
    EnKFStatesDir_out(:,iE) = T_tmp.States_WF(:,2);
    EnKFStates_OP_out(:,iE) = T_tmp.States_OP(:);
    EnKFStates_T_out(:,iE)  = T_tmp.States_T(:);

    % Thrust
    idxTh = find(contains(M.Properties.VariableNames,"Thrust","IgnoreCase",true),1);
    if isempty(idxTh)
        error('No Thrust column found in M.');
    end

    Mnum = [ ...
        M.("Time [s]"), ...
        M.("Foreign Reduction [%]"), ...
        M.("Added turbulence [%]"), ...
        M.("Effective Wind Speed [ms^-1]"), ...
        M.("Free Wind Speed [ms^-1]"), ...
        M{:, idxTh}, ...
        M.("Wake centerline Y [m]"), ...
        M.("Wake centerline Z [m]") ...
    ];

    EnKFM_out(:,iE)            = Mnum(:);
    EnKFOutputThrust_out(:,iE) = M{end-T.nT+1:end, idxTh};
end

% Raccogli outputs
EnKFStatesVel    = EnKFStatesVel_out;
EnKFStatesDir    = EnKFStatesDir_out;
EnKFStates_OP    = EnKFStates_OP_out;
EnKFStates_T     = EnKFStates_T_out;
EnKFOutputThrust = EnKFOutputThrust_out;
EnKFM            = EnKFM_out;
T.floater.T_EnKF_States = FloaterStates_out;

Mnum_main      = EnKFM_out;
nSimSteps_main = nSimSteps_enkf;

fprintf('\nFLORIDynCL timing over ensembles:\n');
fprintf('  mean = %.6f s | median = %.6f s | max = %.6f s | min = %.6f s\n', ...
    mean(t_FLORIDyn), median(t_FLORIDyn), max(t_FLORIDyn), min(t_FLORIDyn))

% % % no parfor
% for iE = 1:EnKF.nE
%     T_tmp = T;
%     T_tmp.States_WF(:,1) = EnKFStatesVel(:,iE);
%     T_tmp.States_WF(:,2) = EnKFStatesDir(:,iE);
%     T_tmp.States_OP      = reshape(EnKFStates_OP(:,iE),[],EnKFnStatesOP);
%     T_tmp.States_T       = reshape(EnKFStates_T(:,iE),[],EnKFnStatesT);
% 
%     Float_states = T.floater.T_EnKF_States(:,:,iE);   % [nT x 6]
%     T_tmp.floater.x         = Float_states(:,1);
%     T_tmp.floater.y         = Float_states(:,2);
%     T_tmp.floater.vx        = Float_states(:,3);
%     T_tmp.floater.vy        = Float_states(:,4);
%     T_tmp.floater.theta     = Float_states(:,5);
%     T_tmp.floater.theta_dot = Float_states(:,6);
%     T_tmp.floater.theta_abs = Float_states(:,7);
% 
%     % fprintf('[iE=%d/%d] FLORIDynCL START\n', iE, EnKF.nE);
%     % t0 = tic;
%     [T_tmp, M, ~] = FLORIDynCL_floating(T_tmp, Wind, Sim, Con, Vis, paramFLORIDyn, paramFLORIS);
%     % t_FLORIDyn(iE) = toc(t0);
%     % fprintf('[iE=%d/%d] FLORIDynCL END   (%.6f s)\n', iE, EnKF.nE, t_FLORIDyn(iE));
% 
%     T.floater.T_EnKF_States(:,1,iE) = T_tmp.floater.x;
%     T.floater.T_EnKF_States(:,2,iE) = T_tmp.floater.y;
%     T.floater.T_EnKF_States(:,3,iE) = T_tmp.floater.vx;
%     T.floater.T_EnKF_States(:,4,iE) = T_tmp.floater.vy;
%     T.floater.T_EnKF_States(:,5,iE) = T_tmp.floater.theta;
%     T.floater.T_EnKF_States(:,6,iE) = T_tmp.floater.theta_dot;
%     T.floater.T_EnKF_States(:,7,iE) = T_tmp.floater.theta_abs;
% 
% 
%     % if iE == 1
%     %     disp(M.Properties.VariableNames');
%     % end
% 
%     EnKFStatesVel(:,iE) = T_tmp.States_WF(:,1);
%     EnKFStatesDir(:,iE) = T_tmp.States_WF(:,2);
%     EnKFStates_OP(:,iE) = T_tmp.States_OP(:);
%     EnKFStates_T(:,iE)  = T_tmp.States_T(:);
% 
%     % trova colonna thrust in modo robusto
%     idxTh = find(contains(M.Properties.VariableNames,"Thrust","IgnoreCase",true),1);
% 
%     % prendi anche le altre colonne (serve per EnKF.M{iE}, es. Foreign Reduction)
%     Mnum = [ ...
%         M.("Time [s]"), ...
%         M.("Foreign Reduction [%]"), ...
%         M.("Added turbulence [%]"), ...
%         M.("Effective Wind Speed [ms^-1]"), ...
%         M.("Free Wind Speed [ms^-1]"), ...
%         M{:, idxTh}, ...   % thrust colonna trovata
%         M.("Wake centerline Y [m]"), ...
%         M.("Wake centerline Z [m]"), ...
%     ];
% 
%     EnKFM(:,iE) = Mnum(:);
% 
%     % ultimo step thrust (ultime nT righe, colonna 6)
%     EnKFOutputThrust(:,iE) = Mnum(end-T.nT+1:end, 6);
% 
% 
%     if isempty(idxTh)
%         error('No Thrust column found in M.');
%     end
%     EnKFOutputThrust(:,iE) = M{end-T.nT+1:end, idxTh};
% end
% 
% fprintf('\nFLORIDynCL timing over ensembles:\n');
% fprintf('  mean = %.6f s | median = %.6f s | max = %.6f s | min = %.6f s\n', ...
%     mean(t_FLORIDyn), median(t_FLORIDyn), max(t_FLORIDyn), min(t_FLORIDyn))
% 
% 
%  % ===================================================


% Limit
if EnKF.Vel.Limits
    EnKFStatesVel(EnKFStatesVel<EnKF.Vel.Min) = EnKF.Vel.Min;
    EnKFStatesVel(EnKFStatesVel>EnKF.Vel.Max) = EnKF.Vel.Max;
end

% Collect parallel computing
EnKF.States.Vel = EnKFStatesVel;
EnKF.States.Dir = EnKFStatesDir;
EnKF.States_OP  = reshape(EnKFStates_OP,[], EnKF.nE * EnKF.nStatesOP);
EnKF.States_T   = reshape(EnKFStates_T, [], EnKF.nE * EnKF.nStatesT);
% EnKF.Output.Pow = EnKFOutputPow;
EnKF.Output.Thrust = EnKFOutputThrust;
for iE = 1:EnKF.nE
    % EnKF.M{iE} = [EnKF.M{iE}; ...
    %     array2table(reshape(EnKFM(:,iE),Sim.nSimSteps * T.nT,[]),...
    %     'VariableNames',...
    %     {'Time [s]','Foreign Reduction [%]','Added turbulence [%]',...
    %     'Effective Wind Speed [ms^-1]','Free Wind Speed [ms^-1]', ...
    %     'Thrust generated [MN]'})];

        % version with centerline postion: 
    EnKF.M{iE} = [EnKF.M{iE}; ...
        array2table(reshape(EnKFM(:,iE),Sim.nSimSteps * T.nT,[]),...
        'VariableNames',...
        {'Time [s]','Foreign Reduction [%]','Added turbulence [%]',...
        'Effective Wind Speed [ms^-1]','Free Wind Speed [ms^-1]', ...
        'Thrust generated [MN]', 'Wake Centerline Y [m]', 'Wake Centerline Z [m]'})];
end
time_simulate = toc;
clear EnKFStatesVel EnKFStatesDir EnKFStates_OP EnKFStates_T EnKFM
%% Combination
tic
% Calculate "true" state & state error covariance
if EnKF.Vel.Correct
    % ===== Projection on the mean "true" state
    %           Velocity
    [EnKF,~]        = EnKF_projectOntoTrueState(EnKF,Sim,T,1);
    %           Direction
    [EnKF, truePos] = EnKF_projectOntoTrueState(EnKF,Sim,T,2);

    % Calc C_phi and C_u based on true states
    [C_u, C_phi] = EnKF_calcC_uC_phi(truePos,...
        mean(EnKF.States.Dir,2), Sim.Dyn, T.posBase, T.nOP,...
        Sim.TimeStep);
    
    % ===== State Error covariance matrix
    C_xx_Vel = EnKF_CalcStateCov(EnKF.States.Vel,EnKF.nE,T.StartI);

    % % ===== Thrust based on projected states
    % Thrust = EnKF_calcThrust(C_u, EnKF, T.D, paramFLORIS, T.StartI, T.floater.theta, T.floater.VCpCt.VCpCt);  
    % 
    % ===== State-to-Output & Output error covariance matrices
    [C_xy_Vel, C_yy_Vel] = ...
        EnKF_CalcStateCov(EnKF.States.Vel, EnKF.nE,...
         'Method', 'InputOutput', 'Output',EnKF.Output.Thrust); %It was P;
end

if EnKF.Dir.Correct
    C_xx_Dir = EnKF_CalcStateCov(EnKF.States.Dir,EnKF.nE,T.StartI);
end

if EnKF.TI.Correct
    C_xx_TI = EnKF_CalcStateCov(EnKF.States.TI,EnKF.nE,T.StartI);
end
time_combine = toc;
%% Get measurements from validation source
%   Calls the functions to get wind speed, direction and amb.
%   turbulence as the normal simulation would and gets [U,phi,TI]
[d, EnKF] = EnKF_GetWFInputs(EnKF,Sim,T,paramFLORIS);
d_phi = d(:,2);
% Calculate reduction of free wind speed at the turbines & apply to meas.
% r = EnKF_GetAveForeignReduction(EnKF.M,T.nT); 
% d(:,1) = d(:,1)./(r*10^(-2)); 

if thrustFF == 0
    % Running closed loop
    %   Retrieve power generated from ZeroMQ message
    %   (Air density error correction not needed!)
    d_T   = generatorPowerArray'./(10^6);

    %   Retrieve wind direction from probe data
    if startsWith(EnKF.Dir.MeasurementSource, 'CLC Probes ')
        try
            d_phi = CLC_read_probe(EnKF.Dir.PathToProbeFile, ...
                EnKF.Dir.MeasurementAverageTime, EnKF.Dir.MeasurementDefault)';
        catch
            disp('ERROR READING WIND DIRECTION PROBE FILE!')
            d_phi = ones(T.nT,1)*EnKF.Dir.MeasurementDefault;
        end
        if endsWith(EnKF.Dir.MeasurementSource, 'Averaged')
            d_phi(:) = mean(d_phi);
        end
    end
else
    % Running offline
    d_T = interp1(thrustFF(1:T.nT:end,2),...
        reshape(thrustFF(:,3),T.nT,[])',Sim.EndTime)';
end


% ===== UNIT CHECK: measured thrust vs predicted thrust (ensemble) =====

if EnKF.Vel.Correct
    % predicted thrust of first ensemble (units must match d_T)
    T_pred_1 = EnKF.Output.Thrust(:,1);

    innov_1 = d_T - T_pred_1;

    fprintf('\n[UNIT CHECK] per-turbine values (measured vs predicted)\n');
    for i = 1:numel(d_T)
        fprintf('  T%02d: d_T=%g   T_pred(1)=%g   mismatch=%g\n', ...
            i, d_T(i), T_pred_1(i), innov_1(i));
    end

end


%% Correct Ensemble states
% [1] Eq.4.37
% x_a,j = x_f,j + C_f,xx * M^T (M C_f,xx M^T + C_ee)^-1 * (d_j-M*x_f,j)
% x_a,j = x_f,j + K * (d_j-M*x_f,j)
tic
for iE = 1:EnKF.nE

    if sum([EnKF.Vel.loc, EnKF.Dir.loc, EnKF.TI.loc])>0
        % Calculate the distance of the OPs to the turbines
        OPs_tmp = ...
            EnKF.States_OP(:,EnKF.nStatesOP*(iE-1)+1:EnKF.nStatesOP*(iE-1)+3);
        distOPs = EnKF_distOPs(OPs_tmp,T.StartI);
    end
    
    % ========= Correct Velocity =========
    if EnKF.Vel.Correct
        % Calculate locaization
        LocCov = GaspariAndCohn1999(EnKF.Vel.cutOffLength,distOPs);

        % Set measurements and pollute with C_ee_Vel
        % d_Vel_j = d(:,1) + ...
        %     (randn(1,T.nT)*EnKF.Vel.C_ee_Vel_Chol)';
        
        % Kalman gain based on the mismatch in Thrust generated
        K = EnKF_CalcKalmanGain('cxy_cyy_cee_loc', ...
            'Cxy',C_xy_Vel,'Cyy',C_yy_Vel,'Cee',EnKF.Vel.C_ee_Vel,...
            'Loc',LocCov);

        noise = (randn(1,T.nT)*EnKF.Output.C_ee_Thrust_Chol)';
        d_T_j = d_T + noise;

        % Correction of the state
        innovation = (d_T_j - EnKF.Output.Thrust(:,iE));
        EnKF.States.Vel(:,iE) = EnKF.States.Vel(:,iE) + ...
            K * innovation;% + ...
            %(randn(1,T.nT)*EnKF.Output.C_ee_Thrust_Chol)');

        % EnKF.Output.Thrust IS NOT SET!!! IS ZERO!!!

    end

    % ========= Correct Direction =========
    if EnKF.Dir.Correct
        % Set measurements and pollute with C_ee_Dir
        d_Dir_j = d_phi + ...
            (randn(1,T.nT)*EnKF.Dir.C_ee_Dir_Chol)';

        % create localization covariance and multiply with state
        % error covariance matrix
        LocCov = GaspariAndCohn1999(EnKF.Dir.cutOffLength,distOPs);
        K = EnKF_CalcKalmanGain('c_cxx_cee_loc', ...
            'Cxx',C_xx_Dir,'Cee',EnKF.Dir.C_ee_Dir,...
            'Loc',LocCov,'C',C_phi);

        EnKF.States.Dir(:,iE) = EnKF.States.Dir(:,iE) + ...
            K * (d_Dir_j - EnKF.States.Dir(T.StartI,iE));

    end

    % ========= Correct TI =========
    if EnKF.TI.Correct
        % Set measurements
        d_TI_j = d(:,3) + ...
            (randn(1,T.nT)*EnKF.TI.C_ee_TI_Chol)';

        if EnKF.TI.loc
            % create localization covariance and multiply with state
            % error covariance matrix
            LocCov = GaspariAndCohn1999(EnKF.TI.cutOffLength,distOPs);

            K = (LocCov.*C_xx_TI)/(LocCov(T.StartI,:).*...
                C_xx_TI(T.StartI,:) + EnKF.TI.C_ee_TI);

            % Plot Localization influence
            if and(iS == 1, iE==1)
                EnKF_Vis = plotK_Localization(EnKF_Vis, T, Vis, ...
                    EnKF.TI.cutOffLength, 2);
            end
        else
            % No localization
            K = C_xx_TI/(C_xx_TI(T.StartI,:) + EnKF.TI.C_ee_TI);
        end

        % Plotting of the combined correction
        EnKF_Vis = plotK_combined(EnKF_Vis, T.nT, EnKF, Vis, K,...
            (d_TI_j - EnKF.States.TI(T.StartI,iE)),iE,3);

        % Correction
        EnKF.States.TI(:,iE) = EnKF.States.TI(:,iE) + ...
            K * (d_TI_j - EnKF.States.TI(T.StartI,iE));
    end

end
time_correct = toc;

%% Store data
if EnKF.StoreData
    folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName ...
        filesep 'EnKF_Thrust'];
    if ~exist(folder_name, "dir")
        mkdir(folder_name)
    end

    % Store Thrust
    writematrix(d_T,[folder_name filesep...
        pad(num2str(time.current),7,"left",'0') '_thrust_measured.txt'])
    writematrix(EnKF.Output.Thrust,[folder_name filesep...
        pad(num2str(time.current),7,"left",'0') '_thrust_predicted.txt'])
    writematrix(EnKF_calcThrust(C_u, EnKF, T.D, paramFLORIS, T.StartI, T.floater.theta, T.floater.VCpCt.VCpCt),...
        [folder_name filesep...
        pad(num2str(time.current),7,"left",'0') '_thrust_corrected.txt'])
    
   % Store Wind direction
    writematrix(d_phi,[folder_name filesep...
        pad(num2str(time.current),7,"left",'0') '_dir_measured.txt'])

    % Store time
    folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName ...
        filesep 'EnKF_Time'];
    if ~exist(folder_name, "dir")
        mkdir(folder_name)
    end
    writematrix([time.current time_simulate time_combine time_correct],...
        [folder_name filesep 'time_sim_com_cor.csv'],"WriteMode","append")
    
end

% Clean up
clear K d_Dir_j d_TI_j d_Vel_j d C_xx_TI C_xx_Dir LocCov C_xy_Vel C_yy_Vel C_u C_phi d_P M time_simulate time_combine time_correct d_phi