function [res_hat, res_std] = residualTowerGPR_wd( ...
    gprS, DEL_NN, DD, YD, gamma_upstream, yaw0_deg, ...
    Hs, Tp, wdir, mean_wind, TI_wind, turb, sin_wind_dir, cos_wind_dir)
% GPR residual corrector for the tower DEL, "wind_dir" version (13 features).
% Drop-in replacement of residualTowerNN: same arguments, same first output
% (res_hat is subtracted from the NN DEL by the caller: DEL = DEL_NN - res_hat).
%
% Optional 2nd output res_std = predictive standard deviation of the residual [same units as res_hat].
% It needs the Cholesky factor gprS.L (n_train x n_train) and is computed ONLY when nargout > 1,
% so the control path (mean only) stays cheap.
%
% Model: sklearn Pipeline(StandardScaler + GaussianProcessRegressor) exported by
%   Simulations/2026_3T_DEL/pts_joblibs/joblib_converter_for_GPR.py
%   kernel = ConstantKernel * RBF(ARD length scales) + WhiteKernel, normalize_y = True
% gprS = load('.../pts_joblibs/gpr_residual_tower_wind_dir.mat', fields used below).
%
% std definition = sklearn predict(return_std=True): k(x,x) includes the WhiteKernel noise,
%   var_norm = const + noise - ||L \ k(x,X)||^2 ,  std = sqrt(var_norm) * y_train_std
%
% Feature order must match training:
% ["DEL_NN", "DD", "YD", "gamma_upstream", "yaw0_deg",
%  "Hs", "Tp", "wvdir", "mean_wind", "TI_wind", "turb", "sin_wind_dir", "cos_wind_dir"]

    x = [DEL_NN, DD, YD, gamma_upstream, yaw0_deg, ...
         Hs, Tp, wdir, mean_wind, TI_wind, turb, sin_wind_dir, cos_wind_dir];

    % 1) Standardize input with the exported scaler
    x_scaled = (x - gprS.x_mean(:).') ./ gprS.x_scale(:).';   % 1 x 13

    % 2) Cross-kernel k(x, X_train) = c * exp(-0.5 * sum(((x - xi)/l).^2))
    %    (WhiteKernel only acts on the training diagonal, not here)
    DX      = (gprS.X_train_scaled - x_scaled) ./ gprS.length_scale(:).';   % n_train x 13
    k_trans = gprS.const_value(1) * exp(-0.5 * sum(DX.^2, 2));              % n_train x 1

    % 3) Posterior mean in normalized y-space, then undo normalize_y
    y_mean_norm = k_trans.' * gprS.alpha(:);
    if gprS.normalize_y(1) == 1
        res_hat = y_mean_norm * gprS.y_train_std(1) + gprS.y_train_mean(1);
        y_scale = gprS.y_train_std(1);
    else
        res_hat = y_mean_norm;
        y_scale = 1;
    end

    % 4) Predictive std (only if requested)
    if nargout > 1
        v          = linsolve(gprS.L, k_trans, struct('LT', true));   % L is lower triangular
        y_var_norm = gprS.const_value(1) + gprS.noise_level(1) - (v.' * v);
        res_std    = sqrt(max(y_var_norm, 0)) * y_scale;
    end
end
