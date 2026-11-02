function [T, EnKF, EnKFM_out, nSimSteps_enkf] = CLC_EnKF_Par_floating_fun(T, EnKF, Sim, Con, Wind, Vis, paramFLORIS, paramFLORIDyn, time, GPR_tower_const, GPR_blade_const, generatorPowerArray, thrustFF, CLC, WD_meas)

% CLC EnKF function - converted from script
tic

% Set inputs
Sim.StartTime   = time.current - time.stepEnKF;
Sim.EndTime     = time.current;
Sim.nSimSteps   = EnKF.nS;

% Prepare parallel computing
EnKFStatesVel    = EnKF.States.Vel;
EnKFStatesDir    = EnKF.States.Dir;
EnKFStates_OP    = reshape(EnKF.States_OP,[],EnKF.nE);
EnKFStates_T     = reshape(EnKF.States_T,[],EnKF.nE);
EnKFnStatesOP    = EnKF.nStatesOP;
EnKFnStatesT     = EnKF.nStatesT;
EnKFOutputThrust = EnKF.Output.Thrust;
EnKFM            = zeros(Sim.nSimSteps * T.nT * 8, EnKF.nE);

t_FLORIDyn = zeros(EnKF.nE,1);
nSimSteps_enkf = Sim.nSimSteps;

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
    T_tmp.GPR_corrector_tower = GPR_tower_const.Value;
    T_tmp.GPR_corrector_blade = GPR_blade_const.Value;
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

    Float_out_ie = zeros(T.nT, 7);
    Float_out_ie(:,1) = T_tmp.floater.x;
    Float_out_ie(:,2) = T_tmp.floater.y;
    Float_out_ie(:,3) = T_tmp.floater.vx;
    Float_out_ie(:,4) = T_tmp.floater.vy;
    Float_out_ie(:,5) = T_tmp.floater.theta;
    Float_out_ie(:,6) = T_tmp.floater.theta_dot;
    Float_out_ie(:,7) = T_tmp.floater.theta_abs;
    FloaterStates_out(:,:,iE) = Float_out_ie;

    EnKFStatesVel_out(:,iE) = T_tmp.States_WF(:,1);
    EnKFStatesDir_out(:,iE) = T_tmp.States_WF(:,2);
    EnKFStates_OP_out(:,iE) = T_tmp.States_OP(:);
    EnKFStates_T_out(:,iE)  = T_tmp.States_T(:);

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
T.floater.T_EnKF_States = FloaterStates_out;

fprintf('\nFLORIDynCL timing over ensembles:\n');
fprintf('  mean = %.6f s | median = %.6f s | max = %.6f s | min = %.6f s\n', ...
    mean(t_FLORIDyn), median(t_FLORIDyn), max(t_FLORIDyn), min(t_FLORIDyn))

% Limit
if EnKF.Vel.Limits
    EnKFStatesVel(EnKFStatesVel<EnKF.Vel.Min) = EnKF.Vel.Min;
    EnKFStatesVel(EnKFStatesVel>EnKF.Vel.Max) = EnKF.Vel.Max;
end

% Collect
EnKF.States.Vel    = EnKFStatesVel;
EnKF.States.Dir    = EnKFStatesDir;
EnKF.States_OP     = reshape(EnKFStates_OP, [], EnKF.nE * EnKF.nStatesOP);
EnKF.States_T      = reshape(EnKFStates_T,  [], EnKF.nE * EnKF.nStatesT);
EnKF.Output.Thrust = EnKFOutputThrust;

for iE = 1:EnKF.nE
    EnKF.M{iE} = [EnKF.M{iE}; ...
        array2table(reshape(EnKFM_out(:,iE), Sim.nSimSteps * T.nT, []), ...
        'VariableNames', ...
        {'Time [s]','Foreign Reduction [%]','Added turbulence [%]', ...
        'Effective Wind Speed [ms^-1]','Free Wind Speed [ms^-1]', ...
        'Thrust generated [MN]','Wake Centerline Y [m]','Wake Centerline Z [m]'})];
end
time_simulate = toc;

%% Combination
tic
if EnKF.Vel.Correct
    [EnKF,~]        = EnKF_projectOntoTrueState(EnKF,Sim,T,1);
    [EnKF, truePos] = EnKF_projectOntoTrueState(EnKF,Sim,T,2);
    [C_u, C_phi]    = EnKF_calcC_uC_phi(truePos, mean(EnKF.States.Dir,2), Sim.Dyn, T.posBase, T.nOP, Sim.TimeStep);
    C_xx_Vel        = EnKF_CalcStateCov(EnKF.States.Vel, EnKF.nE, T.StartI);
    [C_xy_Vel, C_yy_Vel] = EnKF_CalcStateCov(EnKF.States.Vel, EnKF.nE, 'Method', 'InputOutput', 'Output', EnKF.Output.Thrust);
end
if EnKF.Dir.Correct
    C_xx_Dir = EnKF_CalcStateCov(EnKF.States.Dir, EnKF.nE, T.StartI);
end
if EnKF.TI.Correct
    C_xx_TI = EnKF_CalcStateCov(EnKF.States.TI, EnKF.nE, T.StartI);
end
time_combine = toc;


%% Get measurements da FAST.Farm
d_T   = generatorPowerArray' ./ (10^6);  % thrust: N → MN
d_phi = WD_meas * ones(T.nT, 1);         % wind direction da FAST.Farm [deg]

if EnKF.Vel.Correct
    T_pred_1 = EnKF.Output.Thrust(:,1);
    innov_1  = d_T - T_pred_1;
    fprintf('\n[UNIT CHECK] per-turbine values (measured vs predicted)\n');
    for i = 1:numel(d_T)
        fprintf('  T%02d: d_T=%g   T_pred(1)=%g   mismatch=%g\n', i, d_T(i), T_pred_1(i), innov_1(i));
    end
end

%% Correct Ensemble states
tic
for iE = 1:EnKF.nE
    if sum([EnKF.Vel.loc, EnKF.Dir.loc, EnKF.TI.loc]) > 0
        OPs_tmp  = EnKF.States_OP(:, EnKF.nStatesOP*(iE-1)+1:EnKF.nStatesOP*(iE-1)+3);
        distOPs  = EnKF_distOPs(OPs_tmp, T.StartI);
    end

    if EnKF.Vel.Correct
        LocCov = GaspariAndCohn1999(EnKF.Vel.cutOffLength, distOPs);
        K = EnKF_CalcKalmanGain('cxy_cyy_cee_loc', ...
            'Cxy',C_xy_Vel,'Cyy',C_yy_Vel,'Cee',EnKF.Vel.C_ee_Vel,'Loc',LocCov);
        noise  = (randn(1,T.nT) * EnKF.Output.C_ee_Thrust_Chol)';
        d_T_j  = d_T + noise;
        innovation = d_T_j - EnKF.Output.Thrust(:,iE);
        EnKF.States.Vel(:,iE) = EnKF.States.Vel(:,iE) + K * innovation;
    end

    if EnKF.Dir.Correct
        LocCov  = GaspariAndCohn1999(EnKF.Dir.cutOffLength, distOPs);
        d_Dir_j = d_phi + (randn(1,T.nT) * EnKF.Dir.C_ee_Dir_Chol)';
        K = EnKF_CalcKalmanGain('c_cxx_cee_loc', ...
            'Cxx',C_xx_Dir,'Cee',EnKF.Dir.C_ee_Dir,'Loc',LocCov,'C',C_phi);
        EnKF.States.Dir(:,iE) = EnKF.States.Dir(:,iE) + K * (d_Dir_j - EnKF.States.Dir(T.StartI,iE));
    end

    if EnKF.TI.Correct
        d_TI_j = d(:,3) + (randn(1,T.nT) * EnKF.TI.C_ee_TI_Chol)';
        if EnKF.TI.loc
            LocCov = GaspariAndCohn1999(EnKF.TI.cutOffLength, distOPs);
            K = (LocCov.*C_xx_TI) / (LocCov(T.StartI,:).*C_xx_TI(T.StartI,:) + EnKF.TI.C_ee_TI);
        else
            K = C_xx_TI / (C_xx_TI(T.StartI,:) + EnKF.TI.C_ee_TI);
        end
        EnKF.States.TI(:,iE) = EnKF.States.TI(:,iE) + K * (d_TI_j - EnKF.States.TI(T.StartI,iE));
    end
end
time_correct = toc;

%% Store data
if EnKF.StoreData
    folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName filesep 'EnKF_Thrust'];
    if ~exist(folder_name, "dir"); mkdir(folder_name); end
    writematrix(d_T, [folder_name filesep pad(num2str(time.current),7,"left",'0') '_thrust_measured.txt'])
    writematrix(EnKF.Output.Thrust, [folder_name filesep pad(num2str(time.current),7,"left",'0') '_thrust_predicted.txt'])
    writematrix(d_phi, [folder_name filesep pad(num2str(time.current),7,"left",'0') '_dir_measured.txt'])

    folder_name = [Sim.PathToSim 'Results' filesep CLC.CaseName filesep 'EnKF_Time'];
    if ~exist(folder_name, "dir"); mkdir(folder_name); end
    writematrix([time.current time_simulate time_combine time_correct], ...
        [folder_name filesep 'time_sim_com_cor.csv'], "WriteMode", "append")
end

end