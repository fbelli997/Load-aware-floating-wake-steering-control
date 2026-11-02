% ======================================================================= %
% SENSITIVITY ANALYSIS - CoL_DEL controller
% Sweeps lambda_tower x lambda_blade (6x6 = 36 runs)
% ~3 ore totali (~5 min per run)
% ======================================================================= %

pathToSimulation = ['2026_3T_DEL' filesep '3T_IEA15_floating_00'];

%% Initialize ONCE
CLC_init_floating;
addpath(genpath("UQ_CLC"))

thrustFF = importSOWFAFile([Sim.PathToSim 'Data' filesep 'FASTfarm_generatedThrust.csv']);
thrustFF(:,3) = thrustFF(:,3)./(10^6);
uSOWFA = 0;

%% Simulation parameters
Sim.StartTime = 0;
Sim.EndTime   = 4600;

% Environmental / operating inputs for DEL models
Hs    = 4;
Tp    = 7;
wvdir = 0;
A     = pi*(T.D/2).^2;

%% Save initial states for restore
T_init    = T;
EnKF_init = EnKF;
Wind_init = Wind;

%% GPR Constants - created ONCE outside the loop
GPR_tower = T.GPR_corrector_tower;
GPR_blade = T.GPR_corrector_blade;
T_init.GPR_corrector_tower = [];
T_init.GPR_corrector_blade = [];
T.GPR_corrector_tower      = [];
T.GPR_corrector_blade      = [];
GPR_tower_const = parallel.pool.Constant(GPR_tower);
GPR_blade_const = parallel.pool.Constant(GPR_blade);
clear GPR_tower GPR_blade

%% Lambda sweep definition
lam_tower_values = linspace(0, 1, 6);  % [0, 0.2, 0.4, 0.6, 0.8, 1.0]
lam_blade_values = linspace(0, 1, 6);
nL = numel(lam_tower_values);
nB = numel(lam_blade_values);

%% Preallocate summary results
results_DEL_blade_cum = zeros(nL, nB);
results_DEL_tower_cum = zeros(nL, nB);
results_yaw_travel    = zeros(nL, nB, 3);  % 3 turbines
results_power_gain    = zeros(nL, nB);

%% Fixed controller settings
opt_type  = 'CoL_DEL';
lambda_RA = 1.0;
N_batch   = 10;
lambda_yaw = 0.15;

%% Results folder
res_folder = fullfile(Sim.PathToSim, 'Results', 'sensitivity');
if ~exist(res_folder, 'dir'); mkdir(res_folder); end

fprintf('Starting sensitivity analysis: %d x %d = %d runs\n', nL, nB, nL*nB);
fprintf('Estimated time: ~%.0f minutes\n', nL*nB*5);


%% Run baseline yaw=0 ONCE before the lambda sweep
fprintf('\n========================================\n');
fprintf('Run BASELINE: yaw = 0\n');
fprintf('========================================\n');

T    = T_init;
EnKF = EnKF_init;
Wind = Wind_init;

time.it=0; time.current=0; time.current_nO=0; time.end=4600;
time.stepWorld=0.5; time.stepEnKF=Sim.TimeStep*EnKF.nS; time.stepMPC=Sim.TimeStep*CLC.Time.nS;

ws_hist=[]; wd_hist=[]; ws_std_hist=[]; wd_std_hist=[];
time_hist=[]; time_hist_control=[]; yaw_hist=[];
DEL_NN_hist_tower=[]; DEL_NN_hist_blade=[];
YD_hist=[]; DD_hist=[]; TI_hist=[]; nds_hist=[];
power_bl_hist_baseline=[];

run_controller = true;
while run_controller
    t_loop = tic;
    time.it = time.it + 1;
    if time.current >= time.end; run_controller = false; end

    if mod(time.current_nO, time.stepEnKF)==0 && time.current > 0
        CLC_EnKF_Par_floating;
        [T,M,stdM] = EnKF_getMeanState(EnKF, T, paramFLORIDyn);
        ws_en=mean(EnKF.States.Vel(T.StartI,:),1); ws_hat=mean(ws_en); ws_std=std(ws_en);
        wd_en=mean(EnKF.States.Dir(T.StartI,:),1); wd_hat=mean(wd_en); wd_std=std(wd_en);
        T.GeomYaw.ws_hat=ws_hat; T.GeomYaw.wd_hat=wd_hat;
        T.GeomYaw.ws_std=ws_std; T.GeomYaw.wd_std=wd_std;
        time_hist(end+1,1)=time.current;
    end

    controller_rate = 1 * time.stepEnKF;
    if mod(time.current_nO,controller_rate)==0 && time.current > controller_rate
        % yaw=0 always
        gamma_opt = zeros(1, T.nT);
        yaw_hist(end+1,:) = gamma_opt;

        T.States_T(T.StartI,2) = gamma_opt(:);
        for iE=1:EnKF.nE
            cs=(iE-1)*EnKF.nStatesT+1; ce=iE*EnKF.nStatesT;
            St=EnKF.States_T(:,cs:ce); St(T.StartI,2)=gamma_opt(:);
            EnKF.States_T(:,cs:ce)=St;
        end

        if T.GeomYaw.active==1
            [DD,YD,turb]=LayoutDisp(T.States_WF(T.StartI(1),2),T.posBase(:,1),T.posBase(:,2),T.D(1));
        else
            DD=5*ones(T.nT,1); YD=zeros(T.nT,1); turb=ones(T.nT,1);
        end
        yaw_deg = gamma_opt(:);

        Mnum_4d   = reshape(Mnum_main,T.nT,nSimSteps_main,8,EnKF.nE);
        mean_wind = squeeze(mean(mean(Mnum_4d(:,:,4,:),2),4));
        TI_wind   = squeeze(mean(mean(Mnum_4d(:,:,3,:),2),4))/100 + [0.06,0.03,0.03].';

        for i=1:T.nT
            [DEL_NN,~,~]=predictTowerNN(T.NN_predictor_tower,0,Hs,Tp,wvdir,mean_wind(i),TI_wind(i));
            res_hat=residualTowerGPR(GPR_tower_const.Value,DEL_NN,DD(i),YD(i),0,0,Hs,Tp,wvdir,mean_wind(i),TI_wind(i),turb(i));
            T.floater.DEL_tower(i)=DEL_NN-res_hat;
            [DEL_NN,~,~]=predictBladeNN(T.NN_predictor_blade,0,Hs,Tp,wvdir,mean_wind(i),TI_wind(i));
            res_hat=residualBladeGPR(GPR_blade_const.Value,DEL_NN,DD(i),YD(i),0,0,Hs,Tp,wvdir,mean_wind(i),TI_wind(i),turb(i));
            T.floater.DEL_blade(i)=DEL_NN-res_hat;
        end
        DEL_NN_hist_tower(end+1,:)=T.floater.DEL_tower.';
        DEL_NN_hist_blade(end+1,:)=T.floater.DEL_blade.';
        time_hist_control(end+1,1)=time.current;
    end
    time.current=time.current+time.stepWorld;
    time.current_nO=time.current_nO+time.stepWorld;
end

%% Baseline metrics
DEL_blade_baseline = sum(DEL_NN_hist_blade(:));
DEL_tower_baseline = sum(DEL_NN_hist_tower(:));
yaw_travel_baseline = zeros(1,T.nT);  % sempre zero per yaw=0

save(fullfile(res_folder,'baseline_yaw0.mat'), ...
    'DEL_blade_baseline','DEL_tower_baseline','yaw_travel_baseline', ...
    'DEL_NN_hist_blade','DEL_NN_hist_tower','time_hist_control');
fprintf('Baseline saved.\n');

%% ===== MAIN LOOP =====
for f = 1:nL
    for h = 1:nB

        lambda_tower = lam_tower_values(f);
        lambda_blade = lam_blade_values(h);

        fprintf('\n========================================\n');
        fprintf('Run %d/%d: lambda_tower=%.2f | lambda_blade=%.2f\n', ...
            (f-1)*nB+h, nL*nB, lambda_tower, lambda_blade);
        fprintf('========================================\n');

        %% Restore initial states - no reinit needed, use saved post-init states
        T    = T_init;
        EnKF = EnKF_init;
        Wind = Wind_init;

        %% Reset time
        time.it         = 0;
        time.current    = 0;
        time.current_nO = 0;
        time.end        = 4600;
        time.stepWorld  = 0.5;
        time.stepEnKF   = Sim.TimeStep * EnKF.nS;
        time.stepMPC    = Sim.TimeStep * CLC.Time.nS;

        %% Reset history arrays
        ws_hist=[]; wd_hist=[]; ws_std_hist=[]; wd_std_hist=[];
        time_hist=[]; time_hist_control=[]; yaw_hist=[];
        DEL_NN_hist_tower=[]; DEL_NN_hist_blade=[];
        YD_hist=[]; DD_hist=[]; TI_hist=[]; nds_hist=[];
        power_opt_hist=[]; power_bl_hist=[];
        gamma_prev = zeros(1, T.nT);  
        run_controller = true;

        %% ===== SIMULATION WHILE LOOP =====
        while run_controller

            t_loop  = tic;
            time.it = time.it + 1;
            if time.current >= time.end; run_controller = false; end

            %% EnKF
            if mod(time.current_nO, time.stepEnKF) == 0 && time.current > 0

                CLC_EnKF_Par_floating;

                [T, M, stdM] = EnKF_getMeanState(EnKF, T, paramFLORIDyn);

                ws_en  = mean(EnKF.States.Vel(T.StartI,:), 1);
                ws_hat = mean(ws_en); ws_std = std(ws_en);
                wd_en  = mean(EnKF.States.Dir(T.StartI,:), 1);
                wd_hat = mean(wd_en); wd_std = std(wd_en);

                T.GeomYaw.ws_hat = ws_hat; T.GeomYaw.wd_hat = wd_hat;
                T.GeomYaw.ws_std = ws_std; T.GeomYaw.wd_std = wd_std;

                time_hist(end+1,1)   = time.current;
                ws_hist(end+1,1)     = ws_hat;
                wd_hist(end+1,1)     = wd_hat;
                ws_std_hist(end+1,1) = ws_std;
                wd_std_hist(end+1,1) = wd_std;
            end

            controller_rate = 1 * time.stepEnKF;

            %% Power and DEL matrix
            if mod(time.current_nO, time.stepEnKF) == 0 && time.current > 0

                ws_batch = T.GeomYaw.ws_hat + T.GeomYaw.ws_std * randn(N_batch, 1);
                wd_batch = T.GeomYaw.wd_hat + T.GeomYaw.wd_std * randn(N_batch, 1);

                %% STEP 1: gamma candidates
                gam = zeros(N_batch, T.nT);
                for i = 1:N_batch
                    T.GeomYaw.wd_hat = wd_batch(i);
                    T.GeomYaw.ws_hat = ws_batch(i);
                    [geom_yaw_i, ~, ~] = GeometricYaw(T.GeomYaw, T.posBase(:,1), T.posBase(:,2));
                    gam(i,:) = geom_yaw_i(:)';
                end
                T.GeomYaw.wd_hat = wd_hat;
                T.GeomYaw.ws_hat = ws_hat;
                gam(:, 2) = gam(:, 2) * 0.5;

                %% STEP 2: power matrix + DEL matrix
                power_matrix     = zeros(N_batch, N_batch);
                TI_eff_hist      = zeros(N_batch, N_batch, T.nT);
                T_Ueff_hist      = zeros(N_batch, N_batch, T.nT);
                DEL_blade_matrix = zeros(N_batch, T.nT);
                DEL_tower_matrix = zeros(N_batch, T.nT);
                turb = [1, 2, 3];
                [DD_opt, YD_opt, turb_opt] = LayoutDisp(wd_hat, T.posBase(:,1), T.posBase(:,2), T.D(1));

                for i = 1:N_batch
                    for j = 1:N_batch
                        T_tmp = T;
                        T_tmp.States_T(T.StartI, 2) = gam(i,:)';
                        T_tmp.States_WF(:, 1) = ws_batch(j);
                        T_tmp.States_WF(:, 2) = wd_batch(j);
                        T_tmp.dep    = findTurbineGroups(T_tmp.nT, T_tmp.posBase, T_tmp.States_OP, T_tmp.StartI, T_tmp.States_WF, T_tmp.D, T_tmp.nOP, paramFLORIDyn);
                        T_tmp.intOPs = interpolateOPs(T_tmp.nT, T_tmp.StartI, T_tmp.dep, T_tmp.States_OP, T_tmp.posBase, T_tmp.nOP);
                        [~, Cp, T_Ueff, TI_eff_vec] = RunMultipleFloris(T_tmp, paramFLORIS, Wind, turb);
                        T_Ueff_hist(i,j,:) = T_Ueff;
                        TI_eff_vec = TI_eff_vec + [0.06, 0.03, 0.03]';
                        TI_eff_hist(i,j,:) = TI_eff_vec;
                        yaw_corr  = cosd(gam(i,:)').^(paramFLORIS.p_p);
                        tilt_corr = cosd(T.floater.theta_abs).^(paramFLORIS.p_t);
                        power_matrix(i,j) = sum(0.5.*paramFLORIS.airDen.*Cp.*A.*T_Ueff.^3.*yaw_corr.*tilt_corr);
                    end

                    T_Ueff_mean_i = squeeze(mean(T_Ueff_hist(i,:,:), 2));
                    TI_eff_mean_i = squeeze(mean(TI_eff_hist(i,:,:), 2));

                    for t = 1:T.nT
                        T_Ueff_t = T_Ueff_mean_i(t);
                        if turb_opt(t) == 1
                            yaw0_deg_t = 0; gamma_up_t = 0;
                        else
                            up_t = find(turb_opt == turb_opt(t)-1, 1);
                            yaw0_deg_t = abs(gam(i, up_t));
                            gamma_up_t = abs(gam(i, up_t));
                        end
                        [DEL_NN_b, ~, ~]  = predictBladeNN(T.NN_predictor_blade, gam(i,t), Hs, Tp, wvdir, T_Ueff_t, TI_eff_mean_i(t));
                        res_b             = residualBladeGPR(GPR_blade_const.Value, DEL_NN_b, DD_opt(t), YD_opt(t), gamma_up_t, yaw0_deg_t, Hs, Tp, wvdir, T_Ueff_t, TI_eff_mean_i(t), turb_opt(t));
                        [DEL_NN_tw, ~, ~] = predictTowerNN(T.NN_predictor_tower, gam(i,t), Hs, Tp, wvdir, T_Ueff_t, TI_eff_mean_i(t));
                        res_tw            = residualTowerGPR(GPR_tower_const.Value, DEL_NN_tw, DD_opt(t), YD_opt(t), gamma_up_t, yaw0_deg_t, Hs, Tp, wvdir, T_Ueff_t, TI_eff_mean_i(t), turb_opt(t));
                        DEL_blade_matrix(i,t) = DEL_NN_b  - res_b;
                        DEL_tower_matrix(i,t) = DEL_NN_tw - res_tw;
                    end
                end
            end

            %% CONTROLLER
            if mod(time.current_nO, controller_rate) == 0 && time.current > controller_rate

                upstream_idx = zeros(T.nT, 1);

                %% STEP 3: Baseline at yaw=0
                P_bl         = zeros(N_batch, 1);
                DEL_blade_bl = zeros(T.nT, 1);
                DEL_tower_bl = zeros(T.nT, 1);

                for j = 1:N_batch
                    T_tmp_b = T;
                    T_tmp_b.States_T(:, 2) = 0;
                    T_tmp_b.States_WF(:, 1) = ws_batch(j);
                    T_tmp_b.States_WF(:, 2) = wd_batch(j);
                    T_tmp_b.dep    = findTurbineGroups(T_tmp_b.nT, T_tmp_b.posBase, T_tmp_b.States_OP, T_tmp_b.StartI, T_tmp_b.States_WF, T_tmp_b.D, T_tmp_b.nOP, paramFLORIDyn);
                    T_tmp_b.intOPs = interpolateOPs(T_tmp_b.nT, T_tmp_b.StartI, T_tmp_b.dep, T_tmp_b.States_OP, T_tmp_b.posBase, T_tmp_b.nOP);
                    [~, Cp_bl, T_Ueff_bl, TI_eff_bl] = RunMultipleFloris(T_tmp_b, paramFLORIS, Wind, turb);
                    TI_eff_bl = TI_eff_bl + [0.06, 0.03, 0.03];
                    for t = 1:T.nT
                        [DEL_NN_b_bl, ~, ~]  = predictBladeNN(T.NN_predictor_blade, 0, Hs, Tp, wvdir, T_Ueff_bl(t), TI_eff_bl(t));
                        [DEL_NN_tw_bl, ~, ~] = predictTowerNN(T.NN_predictor_tower, 0, Hs, Tp, wvdir, T_Ueff_bl(t), TI_eff_bl(t));
                        DEL_blade_bl(t) = DEL_blade_bl(t) + DEL_NN_b_bl  / N_batch;
                        DEL_tower_bl(t) = DEL_tower_bl(t) + DEL_NN_tw_bl / N_batch;
                    end
                    tilt_corr = cosd(T_tmp_b.floater.theta_abs).^paramFLORIS.p_t;
                    P_bl(j) = sum(0.5.*paramFLORIS.airDen.*Cp_bl.*A.*T_Ueff_bl.^3.*tilt_corr);
                end

                E_bl              = mean(P_bl);
                std_bl            = std(P_bl);
                DEL_blade_bl_mean = max(DEL_blade_bl);
                DEL_tower_bl_mean = max(DEL_tower_bl);

                %% STEP 4: Cost functions
                J_avg     = zeros(N_batch, 1);
                J_CoL     = zeros(N_batch, 1);
                J_avg_DEL = zeros(N_batch, 1);
                J_CoL_DEL = zeros(N_batch, 1);

                for i = 1:N_batch
                    P_i   = power_matrix(i,:)';
                    E_i   = mean(P_i);
                    std_i = std(P_i);

                    dP      = (E_i - E_bl) / E_bl;
                    dTower  = (mean(DEL_tower_matrix(i,:)) - DEL_tower_bl_mean) / DEL_tower_bl_mean;
                    dBlade  = (mean(DEL_blade_matrix(i,:)) - DEL_blade_bl_mean) / DEL_blade_bl_mean;
                    DEL_pen = lambda_tower * dTower + lambda_blade * dBlade; % - lambda_blade * dBlade here we have the minus cause dBlade is always negative
                    yaw_pen = lambda_yaw * sum((deg2rad(gam(i,:) - gamma_prev)).^2);

                    J_avg(i)     = E_i;
                    J_avg_DEL(i) = dP - DEL_pen;

                    cov_mat  = cov(P_i, P_bl);
                    cov_i_bl = cov_mat(1, 2);
                    mu_diff  = E_i - E_bl;
                    std_diff = sqrt(std_i^2 + std_bl^2 - 2*cov_i_bl);

                    if std_diff < 1e-6
                        J_CoL(i)     = 0.5;
                        J_CoL_DEL(i) = 0.5 + DEL_pen + yaw_pen;
                    else
                        J_CoL(i)     = normcdf(0, mu_diff, std_diff);
                        J_CoL_DEL(i) = J_CoL(i) + DEL_pen + yaw_pen;
                    end
                end

                %% STEP 5: Select optimal gamma
                switch opt_type
                    case 'avg',     [~, idx_opt] = max(J_avg);
                    case 'CoL',     [~, idx_opt] = min(J_CoL);
                    case 'avg_DEL', [~, idx_opt] = max(J_avg_DEL);
                    case 'CoL_DEL', [~, idx_opt] = min(J_CoL_DEL);
                end

                gamma_opt = gam(idx_opt, :);
                gamma_prev = gamma_opt;
                yaw_hist(end+1,:)       = gamma_opt;
                power_opt_hist(end+1,1) = mean(power_matrix(idx_opt,:));
                power_bl_hist(end+1,1)  = E_bl;

                fprintf('  lt=%.2f lb=%.2f | T1=%.2f deg | P_opt=%.3f MW | P_bl=%.3f MW | [%.4f s]\n', ...
                    lambda_tower, lambda_blade, gamma_opt(1), ...
                    mean(power_matrix(idx_opt,:))/1e6, E_bl/1e6, toc(t_loop));

                %% STEP 6: Apply gamma
                T.States_T(T.StartI, 2) = gamma_opt(:);
                for iE = 1:EnKF.nE
                    cs = (iE-1)*EnKF.nStatesT+1; ce = iE*EnKF.nStatesT;
                    St = EnKF.States_T(:, cs:ce);
                    St(T.StartI, 2) = gamma_opt(:);
                    EnKF.States_T(:, cs:ce) = St;
                end

                if T.GeomYaw.active == 1
                    [DD, YD, turb] = LayoutDisp(T.States_WF(T.StartI(1),2), T.posBase(:,1), T.posBase(:,2), T.D(1));
                else
                    DD = 5*ones(T.nT,1); YD = zeros(T.nT,1); turb = ones(T.nT,1);
                end
                yaw_deg = gamma_opt(:);

                Mnum_4d   = reshape(Mnum_main, T.nT, nSimSteps_main, 8, EnKF.nE);
                mean_wind = squeeze(mean(mean(Mnum_4d(:,:,4,:), 2), 4));
                TI_wind   = squeeze(mean(mean(Mnum_4d(:,:,3,:), 2), 4))/100 + [0.06, 0.03, 0.03].';

                %% DEL Evaluation (applied yaw)
                for i = 1:T.nT
                    yaw_i = yaw_deg(i);
                    if upstream_idx(i) == 0
                        gamma_upstream = 0; yaw0_deg = 0;
                    else
                        gamma_upstream = abs(yaw_deg(upstream_idx(i)));
                        yaw0_deg       = abs(yaw_deg(upstream_idx(i)));
                    end
                    nds = turb(i);
                    [DEL_NN, ~, ~] = predictTowerNN(T.NN_predictor_tower, yaw_i, Hs, Tp, wvdir, mean_wind(i), TI_wind(i));
                    res_hat = residualTowerGPR(GPR_tower_const.Value, DEL_NN, DD(i), YD(i), gamma_upstream, yaw0_deg, Hs, Tp, wvdir, mean_wind(i), TI_wind(i), nds);
                    T.floater.DEL_tower(i) = DEL_NN - res_hat;
                    [DEL_NN, ~, ~] = predictBladeNN(T.NN_predictor_blade, yaw_i, Hs, Tp, wvdir, mean_wind(i), TI_wind(i));
                    res_hat = residualBladeGPR(GPR_blade_const.Value, DEL_NN, DD(i), YD(i), gamma_upstream, yaw0_deg, Hs, Tp, wvdir, mean_wind(i), TI_wind(i), nds);
                    T.floater.DEL_blade(i) = DEL_NN - res_hat;
                end

                DEL_NN_hist_tower(end+1,:) = T.floater.DEL_tower.';
                DEL_NN_hist_blade(end+1,:) = T.floater.DEL_blade.';
                YD_hist(end+1,:)           = YD.';
                DD_hist(end+1,:)           = DD.';
                TI_hist(end+1,:)           = TI_wind.';
                nds_hist(end+1,:)          = turb.';
                time_hist_control(end+1,1) = time.current;
            end

            time.current    = time.current    + time.stepWorld;
            time.current_nO = time.current_nO + time.stepWorld;
        end

        fprintf('Run %d/%d done.\n', (f-1)*nB+h, nL*nB);

        %% Cumulative metrics
        DEL_blade_cum = sum(DEL_NN_hist_blade(:));
        DEL_tower_cum = sum(DEL_NN_hist_tower(:));
        if size(yaw_hist, 1) > 1
            yaw_travel = sum(abs(diff(yaw_hist, 1, 1)), 1);
        else
            yaw_travel = zeros(1, T.nT);
        end
        power_gain = sum(power_opt_hist - power_bl_hist);

        %% Store in summary matrices
        results_DEL_blade_cum(f,h) = DEL_blade_cum;
        results_DEL_tower_cum(f,h) = DEL_tower_cum;
        results_yaw_travel(f,h,:)  = yaw_travel;
        results_power_gain(f,h)    = power_gain;

        %% Save individual .mat
        save(fullfile(res_folder, sprintf('lt%.2f_lb%.2f.mat', lambda_tower, lambda_blade)), ...
            'lambda_tower', 'lambda_blade', 'opt_type', ...
            'DEL_blade_cum', 'DEL_tower_cum', 'yaw_travel', 'power_gain', ...
            'yaw_hist', 'DEL_NN_hist_blade', 'DEL_NN_hist_tower', 'time_hist_control', ...
            'power_opt_hist', 'power_bl_hist');
        fprintf('Saved: lt=%.2f lb=%.2f\n', lambda_tower, lambda_blade);

    end
end

%% Save summary
save(fullfile(res_folder, 'all_results.mat'), ...
    'results_DEL_blade_cum', 'results_DEL_tower_cum', ...
    'results_yaw_travel', 'results_power_gain', ...
    'lam_tower_values', 'lam_blade_values', 'opt_type');
fprintf('\nAll %d runs done. Summary saved.\n', nL*nB);

%% ===== SENSITIVITY PLOTS =====
set(0, 'DefaultFigureVisible', 'off');
sp_folder = fullfile(Sim.PathToSim, 'Results', 'sensitivity', 'plots');
if ~exist(sp_folder, 'dir'); mkdir(sp_folder); end

%% Heatmap DEL blade cumulative
figure('Position',[100 100 700 500]);
imagesc(lam_blade_values, lam_tower_values, results_DEL_blade_cum);
colorbar; axis xy;
xlabel('\lambda_{blade}'); ylabel('\lambda_{tower}');
title(['DEL blade cumulativo - ' opt_type]); set(gca,'FontSize',14);
saveas(gcf, fullfile(sp_folder, 'heatmap_DEL_blade.png'));

%% Heatmap DEL tower cumulative
figure('Position',[100 100 700 500]);
imagesc(lam_blade_values, lam_tower_values, results_DEL_tower_cum);
colorbar; axis xy;
xlabel('\lambda_{blade}'); ylabel('\lambda_{tower}');
title(['DEL tower cumulativo - ' opt_type]); set(gca,'FontSize',14);
saveas(gcf, fullfile(sp_folder, 'heatmap_DEL_tower.png'));

%% Heatmap yaw travel T1
figure('Position',[100 100 700 500]);
imagesc(lam_blade_values, lam_tower_values, results_yaw_travel(:,:,1));
colorbar; axis xy;
xlabel('\lambda_{blade}'); ylabel('\lambda_{tower}');
title(['Yaw travel T1 [deg] - ' opt_type]); set(gca,'FontSize',14);
saveas(gcf, fullfile(sp_folder, 'heatmap_yaw_travel_T1.png'));

%% Heatmap yaw travel T2
figure('Position',[100 100 700 500]);
imagesc(lam_blade_values, lam_tower_values, results_yaw_travel(:,:,2));
colorbar; axis xy;
xlabel('\lambda_{blade}'); ylabel('\lambda_{tower}');
title(['Yaw travel T2 [deg] - ' opt_type]); set(gca,'FontSize',14);
saveas(gcf, fullfile(sp_folder, 'heatmap_yaw_travel_T2.png'));

%% Heatmap power gain
figure('Position',[100 100 700 500]);
imagesc(lam_blade_values, lam_tower_values, results_power_gain/1e6);
colorbar; axis xy;
xlabel('\lambda_{blade}'); ylabel('\lambda_{tower}');
title(['Power gain cumulativo [MWh] - ' opt_type]); set(gca,'FontSize',14);
saveas(gcf, fullfile(sp_folder, 'heatmap_power_gain.png'));

%% DEL blade vs lambda_blade for each lambda_tower
figure('Position',[100 100 900 500]);
hold on; grid on; box on;
cmap = lines(nL);
for f = 1:nL
    plot(lam_blade_values, results_DEL_blade_cum(f,:), 'o-', ...
        'Color', cmap(f,:), 'LineWidth', 2, 'MarkerSize', 8, ...
        'DisplayName', sprintf('\\lambda_T=%.2f', lam_tower_values(f)));
end
xlabel('\lambda_{blade}'); ylabel('DEL blade cumulativo');
title(['DEL blade vs \lambda_{blade} - ' opt_type]);
yline(DEL_blade_baseline, 'k--', 'LineWidth', 2, 'DisplayName', 'baseline yaw=0');
legend('Location','best'); set(gca,'FontSize',14);
saveas(gcf, fullfile(sp_folder, 'DEL_blade_vs_lambda_blade.png'));

%% DEL tower vs lambda_tower for each lambda_blade
figure('Position',[100 100 900 500]);
hold on; grid on; box on;
cmap = lines(nB);
for h = 1:nB
    plot(lam_tower_values, results_DEL_tower_cum(:,h), 's-', ...
        'Color', cmap(h,:), 'LineWidth', 2, 'MarkerSize', 8, ...
        'DisplayName', sprintf('\\lambda_B=%.2f', lam_blade_values(h)));
end
xlabel('\lambda_{tower}'); ylabel('DEL tower cumulativo');
title(['DEL tower vs \lambda_{tower} - ' opt_type]);
yline(DEL_tower_baseline, 'k--', 'LineWidth', 2, 'DisplayName', 'baseline yaw=0');
legend('Location','best'); set(gca,'FontSize',14);
saveas(gcf, fullfile(sp_folder, 'DEL_tower_vs_lambda_tower.png'));

%% Yaw travel T1 vs lambda_blade for each lambda_tower
figure('Position',[100 100 900 500]);
hold on; grid on; box on;
cmap = lines(nL);
for f = 1:nL
    plot(lam_blade_values, results_yaw_travel(f,:,1), 'o-', ...
        'Color', cmap(f,:), 'LineWidth', 2, 'MarkerSize', 8, ...
        'DisplayName', sprintf('\\lambda_T=%.2f', lam_tower_values(f)));
end
xlabel('\lambda_{blade}'); ylabel('Yaw travel T1 [deg]');
title(['Yaw travel T1 vs \lambda_{blade} - ' opt_type]);
yline(0, 'k--', 'LineWidth', 2, 'DisplayName', 'baseline yaw=0');
legend('Location','best'); set(gca,'FontSize',14);
saveas(gcf, fullfile(sp_folder, 'yaw_travel_T1_vs_lambda_blade.png'));

%% Power gain vs lambda_blade for each lambda_tower
figure('Position',[100 100 900 500]);
hold on; grid on; box on;
cmap = lines(nL);
for f = 1:nL
    plot(lam_blade_values, results_power_gain(f,:)/1e6, 'o-', ...
        'Color', cmap(f,:), 'LineWidth', 2, 'MarkerSize', 8, ...
        'DisplayName', sprintf('\\lambda_T=%.2f', lam_tower_values(f)));
end
xlabel('\lambda_{blade}'); ylabel('Power gain cumulativo [MWh]');
title(['Power gain vs \lambda_{blade} - ' opt_type]);
legend('Location','best'); set(gca,'FontSize',14);
saveas(gcf, fullfile(sp_folder, 'power_gain_vs_lambda_blade.png'));

%% Pareto: power gain vs DEL blade (scatter, colore = lambda_blade)
figure('Position',[100 100 900 600]);
hold on; grid on; box on;
cmap = parula(nB);
for h = 1:nB
    scatter(results_DEL_blade_cum(:,h), results_power_gain(:,h)/1e6, ...
        80, cmap(h,:), 'filled', ...
        'DisplayName', sprintf('\\lambda_B=%.2f', lam_blade_values(h)));
end
xlabel('DEL blade cumulativo'); ylabel('Power gain cumulativo [MWh]');
title(['Pareto: Power gain vs DEL blade - ' opt_type]);
xline(DEL_blade_baseline, 'k--', 'LineWidth', 2, 'DisplayName', 'baseline yaw=0');
legend('Location','best'); set(gca,'FontSize',14);
colorbar; clim([lam_blade_values(1), lam_blade_values(end)]);
saveas(gcf, fullfile(sp_folder, 'pareto_power_vs_DEL_blade.png'));

%% Pareto: power gain vs DEL tower
figure('Position',[100 100 900 600]);
hold on; grid on; box on;
cmap = parula(nL);
for f = 1:nL
    scatter(results_DEL_tower_cum(f,:), results_power_gain(f,:)/1e6, ...
        80, cmap(f,:), 'filled', ...
        'DisplayName', sprintf('\\lambda_T=%.2f', lam_tower_values(f)));
end
xlabel('DEL tower cumulativo'); ylabel('Power gain cumulativo [MWh]');
title(['Pareto: Power gain vs DEL tower - ' opt_type]);
xline(DEL_tower_baseline, 'k--', 'LineWidth', 2, 'DisplayName', 'baseline yaw=0');  % <-- aggiungi qui
legend('Location','best'); set(gca,'FontSize',14);
saveas(gcf, fullfile(sp_folder, 'pareto_power_vs_DEL_tower.png'));