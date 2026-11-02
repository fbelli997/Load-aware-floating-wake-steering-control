function [gamma_opt, state] = FLORIDyn_step(state, WS_meas, WD_meas, Thrust_meas, Power_meas, TI_meas)
% FLORIDyn_step v5 — closed-loop controller (ZMQ version)
%
% Optimization triggered by:
%   1. Wind direction drifts > 3deg from last opt AND stabilizes (range < 1deg in 36s)
%   2. Reliability check completed (276s after last applied yaw)
%
% Reliability check: once per optimization, 276s after yaw applied.
% Checks power (FAST.Farm vs FLORIS SS) and DEL (surrogate re-evaluated at current conditions).

    %% ===== UNPACK STATE =====
    T               = state.T;
    EnKF            = state.EnKF;
    Sim             = state.Sim;
    Con             = state.Con;
    Wind            = state.Wind;
    Vis             = state.Vis;
    CLC             = state.CLC;
    paramFLORIS     = state.paramFLORIS;
    paramFLORIDyn   = state.paramFLORIDyn;
    time            = state.time;
    N_batch         = state.N_batch;
    Hs              = state.Hs;
    Tp              = state.Tp;
    wvdir           = state.wvdir;
    TI_ambient      = state.TI_ambient;
    A               = state.A;
    gamma_prev      = state.gamma_prev;
    NN_tower_const  = state.NN_tower_const;
    NN_blade_const  = state.NN_blade_const;
    lambda_blade    = state.lambda_blade;
    lambda_tower    = state.lambda_tower;
    lambda_yaw      = state.lambda_yaw;
    lambda          = state.lambda;
    lambda_P        = state.lambda_P;
    opt_type        = state.opt_type;
    scale_factors   = state.scale_factors;
    N_cand          = numel(scale_factors);
    controller_rate = state.controller_rate;   % 276s wake advection time
    t_last_opt      = state.t_last_opt;
    pred            = state.pred;
    wd_buffer       = state.wd_buffer;         % rolling buffer of last 3 wd estimates
    wd_at_last_opt  = state.wd_at_last_opt;    % wd when last optimization was applied

    %% ===== TIME =====
    time.it      = time.it + 1;
    time.current = time.current + time.stepEnKF;

    turb = (1:T.nT)';

    %% ===== EnKF (every 12s) =====
    disp('////////// Run EnKF //////////')

    generatorPowerArray = Thrust_meas(:)' * 1000;   % kN → N
    thrustFF = 0;

    [T, EnKF, EnKFM_out, nSimSteps_enkf] = CLC_EnKF_Par_floating_fun_ZMQ( ...
        T, EnKF, Sim, Con, Wind, Vis, paramFLORIS, paramFLORIDyn, ...
        time, NN_tower_const, NN_blade_const, generatorPowerArray, thrustFF, CLC, WD_meas);

    [T, ~, ~] = EnKF_getMeanState(EnKF, T, paramFLORIDyn);

    ws_en  = mean(EnKF.States.Vel(T.StartI,:), 1);
    ws_hat = mean(ws_en);
    ws_std = std(ws_en);
    wd_en  = mean(EnKF.States.Dir(T.StartI,:), 1);
    wd_hat = mean(wd_en);
    wd_std = std(wd_en);

    wind_dir_training = 270 - wd_hat;
    sin_wd = sin(deg2rad(wind_dir_training));
    cos_wd = cos(deg2rad(wind_dir_training));

    T.GeomYaw.ws_hat = ws_hat;
    T.GeomYaw.wd_hat = wd_hat;
    T.GeomYaw.ws_std = ws_std;
    T.GeomYaw.wd_std = wd_std;

    fprintf('[EnKF @ t=%.1f]  WS=%.2f+-%.2f  WD=%.2f+-%.2f\n', ...
        time.current, ws_hat, ws_std, wd_hat, wd_std);

    %% ===== TI COMPUTATION =====
    TI_meas = TI_meas(:);
    % TI_synth = [TI_ambient; TI_ambient; TI_ambient];
    % 
    % TI_use = sqrt(TI_meas.^2 + TI_synth.^2);


    %% ===== UPDATE WIND DIRECTION BUFFER =====
    % Rolling buffer of last 3 wd estimates (36s window)
    wd_buffer(end+1) = wd_hat;
    if numel(wd_buffer) > 3
        wd_buffer = wd_buffer(end-2:end);
    end

    % Wind drift: current wd vs wd at last optimization
    if isnan(wd_at_last_opt)
        wd_drift  = 0;
    else
        wd_drift  = abs(wd_hat - wd_at_last_opt);
    end

    % Wind stability: range of wd in last 36s < 1deg
    if numel(wd_buffer) >= 3
        wd_stable = (max(wd_buffer) - min(wd_buffer)) < 1.0;
    else
        wd_stable = false;   % not enough data yet
    end

    fprintf('[WD buffer: %.2f %.2f %.2f]  drift=%.2f deg  stable=%d\n', ...
        wd_buffer(max(1,end-2):end), wd_drift, wd_stable);

    % Update history
    state.ws_hist(end+1,1)   = ws_hat;
    state.wd_hist(end+1,1)   = wd_hat;
    state.time_hist(end+1,1) = time.current;
    state.wd_buffer          = wd_buffer;

    %% ===== LOG EnKF ESTIMATES =====
    state.enkf_log(end+1,:) = [time.current, ws_hat, ws_std, wd_hat, wd_std, ...
                                wd_drift, double(wd_stable), double(pred.checked), ...
                                TI_meas(:)']; %TI must be averaged on a bigger window, so, don't take them into account

    % Save EnKF log to Excel every step
    enkf_varnames = {'time_s', 'ws_hat', 'ws_std', 'wd_hat', 'wd_std', ...
                     'wd_drift_deg', 'wd_stable', 'pred_checked'};
    for t = 1:T.nT
        enkf_varnames{end+1} = sprintf('TI_use_T%d', t);
    end
    enkf_table = array2table(state.enkf_log, 'VariableNames', enkf_varnames);
    writetable(enkf_table, fullfile(Sim.PathToSim, 'Results', 'enkf_log_15501.xlsx'));


    %% ===== RELIABILITY CHECK =====
    % Triggered once per optimization, controller_rate seconds after yaw applied.
    do_check = pred.valid && ~pred.checked && ...
               (time.current - pred.t_applied) >= controller_rate;

    if do_check
        fprintf('\n--- RELIABILITY CHECK @ t=%.0f (%.0fs after opt @ t=%.0f) ---\n', ...
            time.current, time.current - pred.t_applied, pred.t_applied);

        %% Power
        P_FF_farm = sum(Power_meas(:)) * 1000;   % kW → W

        nT_local    = T.nT;
        T_Ueff_bl_v = zeros(nT_local, 1);
        TI_bl_v     = zeros(nT_local, 1);
        Cp_bl_v     = zeros(nT_local, 1);

        phi_rot_v = angSOWFA2world(wd_hat);
        R2d_v     = [cos(phi_rot_v), sin(phi_rot_v); -sin(phi_rot_v), cos(phi_rot_v)];
        pos_rot_v = (R2d_v * T.posBase(:,1:2)')';
        [~, sort_v] = sort(pos_rot_v(:,1), 'ascend');

        for k = 1:nT_local
            iT_v   = sort_v(k);
            up_v   = sort_v(1:k-1);
            n_up_v = numel(up_v);
            if n_up_v == 0
                [~, Cp_bl_v(iT_v)] = CalcCt_floating(ws_hat, T.floater.VCpCt.VCpCt);
                T_Ueff_bl_v(iT_v)  = ws_hat;
                TI_bl_v(iT_v)      = 0;   % no upstream -> no wake-added TI
            else
                LocationT_v = [T.posBase(up_v,:); T.posBase(iT_v,:)];
                States_WF_v = repmat([ws_hat, wd_hat, TI_ambient], n_up_v+1, 1); %TI must be averaged on a bigger window, so, don't take them into account

                % TI locale delle turbine upstream = TI gia' accumulata dalla loro scia
                % (come main.m: T.States_T(T.StartI,3) = tmpM(:,2), puro added TI,
                % senza l'ambiente che e' gia' TI0 in States_WF), non zero.
                States_T_v  = [zeros(n_up_v,1), zeros(n_up_v,1), TI_bl_v(up_v); 0, 0, 0];
                D_v         = [T.D(up_v); T.D(iT_v)];
                ptilt_v     = T.floater.theta_abs([up_v(:); iT_v]);
                Ct_v        = zeros(n_up_v, 1);
                for m = 1:n_up_v
                    [Ct_v(m), ~] = CalcCt_floating(T_Ueff_bl_v(up_v(m)), T.floater.VCpCt.VCpCt);
                end
                [~, T_aTI_v, T_Ueff_bl_v(iT_v), ~] = runFLORIS( ...
                    LocationT_v+T.posNac(iT_v,:), States_WF_v, States_T_v, D_v, ...
                    paramFLORIS, Wind.Shear, Ct_v, ptilt_v, k);
                TI_bl_v(iT_v) = sqrt(sum(T_aTI_v.^2));   % pure added TI, no ambient (avoids double-counting TI0)
                [~, Cp_bl_v(iT_v)] = CalcCt_floating(T_Ueff_bl_v(iT_v), T.floater.VCpCt.VCpCt);
            end
        end
        P_bl_check = sum(0.5 .* paramFLORIS.airDen .* Cp_bl_v .* A .* T_Ueff_bl_v.^3);

        dP_check  = (P_FF_farm - P_bl_check) / P_bl_check;
        rel_power = P_FF_farm  / pred.P_pred;
        err_dP    = dP_check   - pred.dP_pred;

        %% DEL re-evaluated at current conditions with applied yaw
        [DD_chk, YD_chk, turb_chk] = LayoutDisp(wd_hat, T.posBase(:,1), T.posBase(:,2), T.D(1));

        T_Ueff_chk = zeros(nT_local, 1);
        TI_add_chk = zeros(nT_local, 1);   % total (ambient+added) TI, used for the DEL-model TI feature
        TIadd_c    = zeros(nT_local, 1);   % pure wake-added TI only, used for States_T (avoids double-counting TI0)

        phi_rot_c = angSOWFA2world(wd_hat);
        R2d_c     = [cos(phi_rot_c), sin(phi_rot_c); -sin(phi_rot_c), cos(phi_rot_c)];
        pos_rot_c = (R2d_c * T.posBase(:,1:2)')';
        [~, sort_c] = sort(pos_rot_c(:,1), 'ascend');

        for k = 1:nT_local
            iT_c   = sort_c(k);
            up_c   = sort_c(1:k-1);
            n_up_c = numel(up_c);
            if n_up_c == 0
                T_Ueff_chk(iT_c) = ws_hat;
                TI_add_chk(iT_c) = 0;
            else
                LocationT_c = [T.posBase(up_c,:); T.posBase(iT_c,:)];
                States_WF_c = repmat([ws_hat, wd_hat, TI_ambient], n_up_c+1, 1); %TI must be averaged on a bigger window, so, don't take them into account
                % TI locale delle turbine upstream = TI gia' accumulata dalla loro scia
                % (come main.m: T.States_T(T.StartI,3) = tmpM(:,2), puro added TI), non zero.
                States_T_c  = [zeros(n_up_c,1), pred.gamma_applied(up_c)', TIadd_c(up_c);
                                0, pred.gamma_applied(iT_c), 0];
                D_c         = [T.D(up_c); T.D(iT_c)];
                ptilt_c     = T.floater.theta_abs([up_c(:); iT_c]);
                Ct_c        = zeros(n_up_c, 1);
                for m = 1:n_up_c
                    [Ct_c(m), ~] = CalcCt_floating(T_Ueff_chk(up_c(m)), T.floater.VCpCt.VCpCt);
                end
                [~, T_aTI_c, T_Ueff_chk(iT_c), ~] = runFLORIS( ...
                    LocationT_c+T.posNac(iT_c,:), States_WF_c, States_T_c, D_c, ...
                    paramFLORIS, Wind.Shear, Ct_c, ptilt_c, k);
                TIadd_c(iT_c)    = sqrt(sum(T_aTI_c.^2));            % pure added, for States_T
                TI_add_chk(iT_c) = sqrt(TI_ambient^2 + sum(T_aTI_c.^2));  % total, for the DEL feature
            end
        end
        TI_eff_chk = TI_add_chk;

        DEL_tower_check = zeros(T.nT, 1);
        DEL_blade_check = zeros(T.nT, 1);

        for t = 1:T.nT
            if turb_chk(t) == 1
                gamma_up_chk = 0;  yaw0_chk = 0;
            else
                up_t = find(turb_chk == turb_chk(t)-1, 1);
                gamma_up_chk = abs(pred.gamma_applied(up_t));
                yaw0_chk     = abs(pred.gamma_applied(up_t));
            end

            [DEL_NN_b_c, ~, ~] = predictBladeNN( ...
                T.NN_predictor_blade, pred.gamma_applied(t), Hs, Tp, wvdir, ...
                T_Ueff_chk(t), TI_eff_chk(t), wind_dir_training);
            % Blade: solo 1o stadio, correttore NN disattivato (era un errore)

            [DEL_NN_tw_c, ~, ~] = predictTowerNN( ...
                T.NN_predictor_tower, pred.gamma_applied(t), Hs, Tp, wvdir, ...
                T_Ueff_chk(t), TI_eff_chk(t), wind_dir_training);
            res_tw_c = residualTowerNN( ...
                NN_tower_const.Value, DEL_NN_tw_c, DD_chk(t), YD_chk(t), ...
                gamma_up_chk, yaw0_chk, Hs, Tp, wvdir, ...
                T_Ueff_chk(t), TI_eff_chk(t), turb_chk(t), sin_wd, cos_wd);

            DEL_blade_check(t) = DEL_NN_b_c;
            DEL_tower_check(t) = DEL_NN_tw_c - res_tw_c;
        end

        rel_DEL_tower = DEL_tower_check ./ pred.DEL_tower_pred(:);
        rel_DEL_blade = DEL_blade_check ./ pred.DEL_blade_pred(:);

        fprintf('  P_FASTfarm:         %.4f MW\n', P_FF_farm/1e6);
        fprintf('  P_pred (FLORIS SS): %.4f MW\n', pred.P_pred/1e6);
        fprintf('  dP_pred:            %.2f%%\n',  pred.dP_pred*100);
        fprintf('  dP_check:           %.2f%%\n',  dP_check*100);
        fprintf('  rel_power:          %.4f\n',    rel_power);
        fprintf('  err_dP:             %.2f%%\n',  err_dP*100);
        for t = 1:T.nT
            fprintf('  T%d: rel_DEL_tower=%.4f  rel_DEL_blade=%.4f\n', ...
                t, rel_DEL_tower(t), rel_DEL_blade(t));
        end
        fprintf('---------------------------------------\n\n');

        state.rel_hist(end+1,:) = [ ...
            pred.t_applied, time.current, controller_rate, ...
            pred.gamma_applied, ...
            pred.ws_hat, pred.wd_hat, ws_hat, wd_hat, ...
            pred.P_pred, pred.P_bl_pred, pred.dP_pred*100, ...
            P_FF_farm, P_bl_check, dP_check*100, ...
            rel_power, err_dP*100, ...
            rel_DEL_tower(:)', rel_DEL_blade(:)'];

        pred.checked = true;
    end


    %% ===== OPTIMIZATION TRIGGER =====
    % Trigger 0: first optimization ever (t_last_opt = -inf)
    trigger_first = isinf(t_last_opt);

    % Trigger 1: wind direction drifted > 3deg AND is now stable (range < 1deg in 36s)
    trigger_drift = (wd_drift > 3.0) && wd_stable;

    % Trigger 2: reliability check just completed (276s passed, check done)
    trigger_post_check = pred.checked && ...
                         (time.current - pred.t_applied) >= controller_rate;

    do_optimize = trigger_first || trigger_drift || trigger_post_check;

    if trigger_first;      fprintf('[TRIGGER] First optimization\n'); end
    if trigger_drift;      fprintf('[TRIGGER] WD drift %.2f deg > 3deg AND stable\n', wd_drift); end
    if trigger_post_check; fprintf('[TRIGGER] Post-check re-optimization\n'); end

    if ~do_optimize
        gamma_opt = gamma_prev;
        fprintf('[Controller] No trigger. WD drift=%.2f deg, stable=%d, checked=%d\n', ...
            wd_drift, wd_stable, pred.checked);
        state.T             = T;
        state.EnKF          = EnKF;
        state.time          = time;
        state.gamma_prev    = gamma_prev;
        state.pred          = pred;
        state.t_last_opt    = t_last_opt;
        state.wd_at_last_opt = wd_at_last_opt;
        return;
    end

    %% ===== CONTROLLER =====
    fprintf('\n========== CONTROLLER @ t=%.0f s ==========\n', time.current);
    fprintf('  Trigger: drift=%.2f deg, post_check=%d\n', wd_drift, trigger_post_check);

    [DD_opt, YD_opt, turb_opt] = LayoutDisp(wd_hat, T.posBase(:,1), T.posBase(:,2), T.D(1));

    %% STEP 1: YAW CANDIDATES
    [geom_yaw_mean, ~, ~] = GeometricYaw(T.GeomYaw, T.posBase(:,1), T.posBase(:,2));

    gam = zeros(N_cand, T.nT);
    for k = 1:N_cand
        gam(k,:) = geom_yaw_mean(:)' * scale_factors(k);
    end

    fprintf('Yaw candidates:\n');
    for k = 1:N_cand
        fprintf('  [sf=%.2f] gamma = [%.2f, %.2f, %.2f] deg\n', ...
            scale_factors(k), gam(k,1), gam(k,2), gam(k,3));
    end

    %% STEP 2: WIND SCENARIOS
    ws_batch = ws_hat + ws_std * randn(N_batch, 1);
    wd_batch = wd_hat + wd_std * randn(N_batch, 1);

    %% STEP 3: POWER + DEL MATRIX (FLORIS SS, posBase)
    power_matrix     = zeros(N_cand, N_batch);
    T_Ueff_hist      = zeros(N_cand, N_batch, T.nT);
    TI_eff_hist      = zeros(N_cand, N_batch, T.nT);
    DEL_blade_matrix = zeros(N_cand, T.nT);
    DEL_tower_matrix = zeros(N_cand, T.nT);

    for i = 1:N_cand
        for j = 1:N_batch
            nT_local  = T.nT;
            T_Ueff_ss = zeros(nT_local, 1);
            TI_FLORIS = zeros(nT_local, 1);
            TIadd_ss  = zeros(nT_local, 1);   % pure wake-added TI only, used for States_T
            Cp_ss     = zeros(nT_local, 1);

            phi_rot = angSOWFA2world(wd_batch(j));
            R2d     = [cos(phi_rot), sin(phi_rot); -sin(phi_rot), cos(phi_rot)];
            pos_rot = (R2d * T.posBase(:,1:2)')';
            [~, sort_idx] = sort(pos_rot(:,1), 'ascend');

            for k = 1:nT_local
                iT       = sort_idx(k);
                upstream = sort_idx(1:k-1);
                n_up     = numel(upstream);

                if n_up == 0
                    [~, Cp_ss(iT)] = CalcCt_floating(ws_batch(j), T.floater.VCpCt.VCpCt);
                    T_Ueff_ss(iT)  = ws_batch(j);
                    % TI_FLORIS(iT)  = TI_meas(1);
                    TI_FLORIS(iT)  = TI_ambient;
                    TIadd_ss(iT)   = 0;   % no upstream -> no wake-added TI
                else
                    LocationT_tmp = [T.posBase(upstream,:); T.posBase(iT,:)];
                    States_WF_tmp = repmat([ws_batch(j), wd_batch(j), TI_FLORIS(1)], n_up+1, 1);
                    % TI locale delle turbine upstream = TI gia' accumulata dalla loro scia
                    % (come main.m: T.States_T(T.StartI,3) = tmpM(:,2), puro added TI), non zero.
                    States_T_tmp  = [zeros(n_up,1), gam(i,upstream)', TIadd_ss(upstream);
                                     0,             gam(i,iT),         0            ];
                    D_tmp         = [T.D(upstream); T.D(iT)];
                    ptilt_tmp     = T.floater.theta_abs([upstream(:); iT]);
                    Ct_up         = zeros(n_up, 1);
                    for m = 1:n_up
                        [Ct_up(m), ~] = CalcCt_floating( ...
                            T_Ueff_ss(upstream(m)), T.floater.VCpCt.VCpCt);
                    end
                    [~, T_aTI_arr, T_Ueff_i, ~] = runFLORIS( ...
                        LocationT_tmp+T.posNac(iT,:), States_WF_tmp, States_T_tmp, D_tmp, ...
                        paramFLORIS, Wind.Shear, Ct_up, ptilt_tmp, k);
                    T_Ueff_ss(iT) = T_Ueff_i;
                    TIadd_ss(iT)  = sqrt(sum(T_aTI_arr.^2));            % pure added, for States_T
                    TI_FLORIS(iT) = sqrt(TI_ambient^2 + sum(T_aTI_arr.^2));
                    [~, Cp_ss(iT)] = CalcCt_floating(T_Ueff_i, T.floater.VCpCt.VCpCt);
                end
            end

            T_Ueff_hist(i,j,:) = T_Ueff_ss;
            TI_eff_vec         = TI_FLORIS;
            TI_eff_hist(i,j,:) = TI_eff_vec;

            yaw_corr = cosd(gam(i,:)').^(paramFLORIS.p_p);
            power_matrix(i,j) = sum( ...
                0.5 .* paramFLORIS.airDen .* Cp_ss .* A .* T_Ueff_ss.^3 .* yaw_corr);
        end

        T_Ueff_mean_i = squeeze(mean(T_Ueff_hist(i,:,:), 2));
        TI_eff_mean_i = squeeze(mean(TI_eff_hist(i,:,:), 2));

        for t = 1:T.nT
            if turb_opt(t) == 1
                yaw0_deg_t = 0;  gamma_up_t = 0;
            else
                up_t       = find(turb_opt == turb_opt(t)-1, 1);
                yaw0_deg_t = abs(gam(i, up_t));
                gamma_up_t = abs(gam(i, up_t));
            end

            [DEL_NN_b, ~, ~] = predictBladeNN( ...
                T.NN_predictor_blade, gam(i,t), Hs, Tp, wvdir, ...
                T_Ueff_mean_i(t), TI_eff_mean_i(t), wind_dir_training);
            % Blade: solo 1o stadio, correttore NN disattivato (era un errore)

            [DEL_NN_tw, ~, ~] = predictTowerNN( ...
                T.NN_predictor_tower, gam(i,t), Hs, Tp, wvdir, ...
                T_Ueff_mean_i(t), TI_eff_mean_i(t), wind_dir_training);
            res_tw = residualTowerNN( ...
                NN_tower_const.Value, DEL_NN_tw, DD_opt(t), YD_opt(t), ...
                gamma_up_t, yaw0_deg_t, Hs, Tp, wvdir, ...
                T_Ueff_mean_i(t), TI_eff_mean_i(t), turb_opt(t), sin_wd, cos_wd);

            DEL_blade_matrix(i,t) = DEL_NN_b;
            DEL_tower_matrix(i,t) = DEL_NN_tw - res_tw;
        end
    end

    %% STEP 4: BASELINE (gamma=0)
    P_bl         = zeros(N_batch, 1);
    DEL_blade_bl = zeros(T.nT, 1);
    DEL_tower_bl = zeros(T.nT, 1);

    for j = 1:N_batch
        nT_local  = T.nT;
        T_Ueff_bl = zeros(nT_local, 1);
        TI_FLORIS_bl = zeros(nT_local, 1);
        TIadd_bl  = zeros(nT_local, 1);   % pure wake-added TI only, used for States_T
        Cp_bl_j   = zeros(nT_local, 1);

        phi_rot = angSOWFA2world(wd_batch(j));
        R2d     = [cos(phi_rot), sin(phi_rot); -sin(phi_rot), cos(phi_rot)];
        pos_rot = (R2d * T.posBase(:,1:2)')';
        [~, sort_bl] = sort(pos_rot(:,1), 'ascend');

        for k = 1:nT_local
            iT       = sort_bl(k);
            upstream = sort_bl(1:k-1);
            n_up     = numel(upstream);

            if n_up == 0
                [~, Cp_bl_j(iT)] = CalcCt_floating(ws_batch(j), T.floater.VCpCt.VCpCt);
                T_Ueff_bl(iT)    = ws_batch(j);
                % TI_FLORIS_bl(iT)    = TI_meas(1);
                TI_FLORIS_bl(iT)    = TI_ambient;
                TIadd_bl(iT)        = 0;   % no upstream -> no wake-added TI
            else
                LocationT_tmp = [T.posBase(upstream,:); T.posBase(iT,:)];
                States_WF_tmp = repmat([ws_batch(j), wd_batch(j), TI_FLORIS_bl(1)], n_up+1, 1);
                % TI locale delle turbine upstream = TI gia' accumulata dalla loro scia
                % (come main.m: T.States_T(T.StartI,3) = tmpM(:,2), puro added TI), non zero.
                States_T_tmp  = [zeros(n_up,1), zeros(n_up,1), TIadd_bl(upstream); 0, 0, 0];
                D_tmp         = [T.D(upstream); T.D(iT)];
                ptilt_tmp     = T.floater.theta_abs([upstream(:); iT]);
                Ct_up         = zeros(n_up, 1);
                for m = 1:n_up
                    [Ct_up(m), ~] = CalcCt_floating( ...
                        T_Ueff_bl(upstream(m)), T.floater.VCpCt.VCpCt);
                end
                [~, T_aTI_bl, T_Ueff_bl_i, ~] = runFLORIS( ...
                    LocationT_tmp+T.posNac(iT,:), States_WF_tmp, States_T_tmp, D_tmp, ...
                    paramFLORIS, Wind.Shear, Ct_up, ptilt_tmp, k);
                T_Ueff_bl(iT) = T_Ueff_bl_i;
                TIadd_bl(iT)     = sqrt(sum(T_aTI_bl.^2));               % pure added, for States_T
                TI_FLORIS_bl(iT) = sqrt(TI_ambient^2 + sum(T_aTI_bl.^2));
                [~, Cp_bl_j(iT)] = CalcCt_floating(T_Ueff_bl_i, T.floater.VCpCt.VCpCt);
            end
        end

        TI_bl_j = TI_FLORIS_bl;
        P_bl(j) = sum(0.5 .* paramFLORIS.airDen .* Cp_bl_j .* A .* T_Ueff_bl.^3);

        for t = 1:T.nT
            [DEL_NN_b_bl, ~, ~] = predictBladeNN( ...
                T.NN_predictor_blade, 0, Hs, Tp, wvdir, ...
                T_Ueff_bl(t), TI_bl_j(t), wind_dir_training);
            % Blade: solo 1o stadio, correttore NN disattivato (era un errore)

            [DEL_NN_tw_bl, ~, ~] = predictTowerNN( ...
                T.NN_predictor_tower, 0, Hs, Tp, wvdir, ...
                T_Ueff_bl(t), TI_bl_j(t), wind_dir_training);
            res_tw_bl = residualTowerNN( ...
                NN_tower_const.Value, DEL_NN_tw_bl, DD_opt(t), YD_opt(t), ...
                0, 0, Hs, Tp, wvdir, ...
                T_Ueff_bl(t), TI_bl_j(t), turb_opt(t), sin_wd, cos_wd);

            DEL_blade_bl(t) = DEL_blade_bl(t) + DEL_NN_b_bl / N_batch;
            DEL_tower_bl(t) = DEL_tower_bl(t) + (DEL_NN_tw_bl - res_tw_bl) / N_batch;
        end
    end

    E_bl             = mean(P_bl);
    std_bl           = std(P_bl);
    DEL_blade_bl_sum = sum(DEL_blade_bl);
    DEL_tower_bl_sum = sum(DEL_tower_bl);

    %% STEP 5: COST FUNCTIONS
    J_avg     = zeros(N_cand, 1);
    J_RA      = zeros(N_cand, 1);
    J_CoL     = zeros(N_cand, 1);
    J_avg_DEL = zeros(N_cand, 1);
    J_RA_DEL  = zeros(N_cand, 1);
    J_CoL_DEL = zeros(N_cand, 1);

    for i = 1:N_cand
        P_i     = power_matrix(i,:)';
        E_i     = mean(P_i);
        std_i   = std(P_i);
        yaw_pen = mean(((gam(i,1:2) - gamma_prev(1:2)) / 20).^2);
        dP      = (E_i - E_bl) / E_bl;

        DEL_tower_yawed_sum = sum(DEL_tower_matrix(i,:));
        DEL_blade_yawed_sum = sum(DEL_blade_matrix(i,:));
        dTower  = (DEL_tower_yawed_sum - DEL_tower_bl_sum) / DEL_tower_bl_sum;
        dBlade  = (DEL_blade_yawed_sum - DEL_blade_bl_sum) / DEL_blade_bl_sum;
        DEL_pen = lambda_tower * dTower + lambda_blade * dBlade;

        J_avg(i)     = E_i;
        J_avg_DEL(i) = dP - DEL_pen;
        J_RA(i)      = E_i - lambda * std_i;
        J_RA_DEL(i)  = dP - (lambda * std_i / E_bl) - DEL_pen;

        cov_mat  = cov(P_i, P_bl);
        cov_i_bl = cov_mat(1,2);
        mu_diff  = E_i - E_bl;
        std_diff = sqrt(std_i^2 + std_bl^2 - 2*cov_i_bl);

        if std_diff < 1e-6
            J_CoL(i)     = 0.5;
            J_CoL_DEL(i) = 0.5 + DEL_pen;
        else
            J_CoL(i)     = normcdf(0, mu_diff, std_diff);
            J_CoL_DEL(i) = J_CoL(i) - lambda_P * dP + DEL_pen; %+ lambda_yaw * yaw_pen;
        end
    end

    %% STEP 6: SELECT OPTIMAL
    switch opt_type
        case 'avg',     [~, idx_opt] = max(J_avg);     cost_name = 'J_avg';
        case 'RA',      [~, idx_opt] = max(J_RA);      cost_name = 'J_RA';
        case 'CoL',     [~, idx_opt] = min(J_CoL);     cost_name = 'J_CoL';
        case 'avg_DEL', [~, idx_opt] = max(J_avg_DEL); cost_name = 'J_avg_DEL';
        case 'RA_DEL',  [~, idx_opt] = max(J_RA_DEL);  cost_name = 'J_RA_DEL';
        case 'CoL_DEL', [~, idx_opt] = min(J_CoL_DEL); cost_name = 'J_CoL_DEL';
        otherwise
            warning('Unknown opt_type, using avg');
            [~, idx_opt] = max(J_avg); cost_name = 'J_avg [default]';
    end

    gamma_opt = gam(idx_opt, :);

    %% Features actually used to compute the DEL at the optimum (for CLC_history logging)
    T_Ueff_opt = squeeze(mean(T_Ueff_hist(idx_opt,:,:), 2))';
    TI_eff_opt = squeeze(mean(TI_eff_hist(idx_opt,:,:), 2))';
    gamma_up_opt = zeros(1, T.nT);
    for t = 1:T.nT
        if turb_opt(t) == 1
            gamma_up_opt(t) = 0;
        else
            up_t = find(turb_opt == turb_opt(t)-1, 1);
            gamma_up_opt(t) = abs(gamma_opt(up_t));
        end
    end

    fprintf('\n=== YAW OPTIMIZATION [%s] ===\n', cost_name);
    fprintf('  Scale factor:   %.2f\n',    scale_factors(idx_opt));
    fprintf('  Gamma opt:      [%.2f, %.2f, %.2f] deg\n', gamma_opt(1), gamma_opt(2), gamma_opt(3));
    fprintf('  BASELINE power: %.4f MW\n', E_bl/1e6);
    fprintf('  E[P] predicted: %.4f MW\n', mean(power_matrix(idx_opt,:))/1e6);
    fprintf('  dP predicted:   %.2f%%\n',  ((mean(power_matrix(idx_opt,:))-E_bl)/E_bl)*100);
    fprintf('  J_CoL:          %.2f%%\n',  J_CoL(idx_opt)*100);
    fprintf('  J_CoL_DEL:      %.6f\n',    J_CoL_DEL(idx_opt));
    fprintf('================================\n\n');

    %% STEP 7: STORE PREDICTION
    pred.valid          = true;
    pred.t_applied      = time.current;
    pred.gamma_applied  = gamma_opt;
    pred.P_pred         = mean(power_matrix(idx_opt,:));
    pred.P_bl_pred      = E_bl;
    pred.dP_pred        = (pred.P_pred - pred.P_bl_pred) / pred.P_bl_pred;
    pred.ws_hat         = ws_hat;
    pred.wd_hat         = wd_hat;
    pred.DEL_tower_pred = DEL_tower_matrix(idx_opt,:)';
    pred.DEL_blade_pred = DEL_blade_matrix(idx_opt,:)';

    %% STEP 8: APPLY GAMMA
    T.States_T(T.StartI, 2) = gamma_opt(:);

    for iE = 1:EnKF.nE
        col_start = (iE-1)*EnKF.nStatesT + 1;
        col_end   = iE*EnKF.nStatesT;
        S = EnKF.States_T(:, col_start:col_end);
        S(T.StartI, 2) = gamma_opt(:);
        EnKF.States_T(:, col_start:col_end) = S;
    end

    gamma_prev      = gamma_opt;
    t_last_opt      = time.current;
    wd_at_last_opt  = wd_hat;   % update reference WD for drift detection
    pred.checked    = false;    % enable next reliability check in 276s

    %% STEP 9: STORE HISTORY
    % Tempo FAST.Farm dell'ottimizzazione: il .py chiama MATLAB da t=60s (1a chiamata = init),
    % quindi t_FF = time.current + 60
    if ~isfield(state,'time_opt_hist'); state.time_opt_hist = []; end
    state.time_opt_hist(end+1,1) = time.current + 60;
    state.yaw_hist(end+1,:)          = gamma_opt;
    state.DEL_NN_hist_tower(end+1,:) = DEL_tower_matrix(idx_opt,:);
    state.DEL_NN_hist_blade(end+1,:) = DEL_blade_matrix(idx_opt,:);
    state.feat_hist(end+1,:)         = [T_Ueff_opt, TI_eff_opt, DD_opt(:)', YD_opt(:)', ...
                                         turb_opt(:)', gamma_up_opt, ...
                                         ws_hat, wd_hat, Hs, Tp, wvdir, ...
                                         wind_dir_training, sin_wd, cos_wd];

    %% SAVE TO EXCEL
    T_yaw   = array2table(state.yaw_hist, ...
        'VariableNames', arrayfun(@(i) sprintf('T%d_yaw_deg', i), 1:T.nT, 'UniformOutput', false));
    T_tower = array2table(state.DEL_NN_hist_tower, ...
        'VariableNames', arrayfun(@(i) sprintf('T%d_DEL_tower', i), 1:T.nT, 'UniformOutput', false));
    T_blade = array2table(state.DEL_NN_hist_blade, ...
        'VariableNames', arrayfun(@(i) sprintf('T%d_DEL_blade', i), 1:T.nT, 'UniformOutput', false));
    feat_names = [ ...
        arrayfun(@(i) sprintf('T%d_WS_eff', i), 1:T.nT, 'UniformOutput', false), ...
        arrayfun(@(i) sprintf('T%d_TI_eff', i), 1:T.nT, 'UniformOutput', false), ...
        arrayfun(@(i) sprintf('T%d_DD', i), 1:T.nT, 'UniformOutput', false), ...
        arrayfun(@(i) sprintf('T%d_YD', i), 1:T.nT, 'UniformOutput', false), ...
        arrayfun(@(i) sprintf('T%d_turb_row', i), 1:T.nT, 'UniformOutput', false), ...
        arrayfun(@(i) sprintf('T%d_yaw_upstream_deg', i), 1:T.nT, 'UniformOutput', false), ...
        {'WS_hat','WD_hat','Hs','Tp','WaveDir','WindDirTraining','sin_wd','cos_wd'}];
    T_feat = array2table(state.feat_hist, 'VariableNames', feat_names);
    writetable([array2table(state.time_opt_hist, 'VariableNames', {'time'}), T_yaw, T_tower, T_blade, T_feat], ...
        fullfile(Sim.PathToSim, 'Results', 'CLC_history_15501.xlsx'));

    if ~isempty(state.rel_hist)
        n_rel_cols = size(state.rel_hist, 2);
        base_names = {'t_applied','t_check','controller_rate', ...
            'gam_T1','gam_T2','gam_T3', ...
            'ws_pred','wd_pred','ws_check','wd_check', ...
            'P_pred_W','P_bl_pred_W','dP_pred_pct', ...
            'P_FASTfarm_W','P_bl_check_W','dP_check_pct', ...
            'rel_power','err_dP_pct'};
        del_names = {};
        for t = 1:T.nT; del_names{end+1} = sprintf('rel_DEL_tower_T%d',t); end
        for t = 1:T.nT; del_names{end+1} = sprintf('rel_DEL_blade_T%d',t); end
        all_names = [base_names, del_names];
        rel_table = array2table(state.rel_hist, 'VariableNames', all_names(1:n_rel_cols));
        writetable(rel_table, fullfile(Sim.PathToSim, 'Results', 'reliability_log_15501.xlsx'));
    end

    %% REPACK STATE
    state.T              = T;
    state.EnKF           = EnKF;
    state.time           = time;
    state.gamma_prev     = gamma_prev;
    state.pred           = pred;
    state.t_last_opt     = t_last_opt;
    state.wd_at_last_opt = wd_at_last_opt;

end