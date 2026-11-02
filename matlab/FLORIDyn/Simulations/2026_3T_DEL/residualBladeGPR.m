% function [res_hat, res_std] = residualTowerGPR_withStd( ...
%     gprS, DEL_NN, DD, YD, gamma_upstream, yaw0_deg, ...
%     Hs, Tp, wdir, mean_wind, TI_wind, turb)
% 
%     % Feature order must match training:
%     % ["DEL_NN", "DD", "YD", "gamma_upstream", "yaw0_deg",
%     %  "Hs", "Tp", "wvdir", "mean_wind", "TI_wind", "turb"]
% 
%     x = [DEL_NN, DD, YD, gamma_upstream, yaw0_deg, ...
%          Hs, Tp, wdir, mean_wind, TI_wind, turb];
% 
%     % -------------------------------------------------
%     % 1) Standardize input using exported scaler
%     % -------------------------------------------------
%     x_mean  = gprS.x_mean(:).';
%     x_scale = gprS.x_scale(:).';
%     x_scaled = (x - x_mean) ./ x_scale;   % 1 x n_features
% 
%     % -------------------------------------------------
%     % 2) Read exported GPR quantities
%     % -------------------------------------------------
%     Xtr = gprS.X_train_scaled;          % n_train x n_features
%     alpha = gprS.alpha(:);              % n_train x 1
%     L = gprS.L;                         % n_train x n_train
% 
%     const_value  = gprS.const_value(1);
%     length_scale = gprS.length_scale(:).';   % 1 x n_features
%     noise_level  = gprS.noise_level(1);
% 
%     % -------------------------------------------------
%     % 3) Cross-kernel k(x, X_train)
%     % Kernel = ConstantKernel * RBF + WhiteKernel
%     % WhiteKernel contributes only on diagonal, not here
%     % -------------------------------------------------
%     DX = (Xtr - x_scaled) ./ length_scale;      % n_train x n_features
%     sqdist = sum(DX.^2, 2);                     % n_train x 1
%     k_trans = const_value * exp(-0.5 * sqdist); % n_train x 1
% 
%     % -------------------------------------------------
%     % 4) Predict mean in normalized y-space
%     % -------------------------------------------------
%     y_mean_norm = k_trans.' * alpha;
% 
%     % -------------------------------------------------
%     % 5) Predict std in normalized y-space
%     % y_var = k(x,x) - v'v , with v = L \ k_trans
%     % -------------------------------------------------
%     v = L \ k_trans;
%     k_xx = const_value + noise_level;
%     y_var_norm = k_xx - (v.' * v);
% 
%     % numerical safety
%     y_var_norm = max(y_var_norm, 0);
% 
%     % -------------------------------------------------
%     % 6) Undo normalize_y if used in sklearn
%     % -------------------------------------------------
%     if isfield(gprS, 'normalize_y') && gprS.normalize_y(1) == 1
%         y_mean = y_mean_norm * gprS.y_train_std(1) + gprS.y_train_mean(1);
%         y_std  = sqrt(y_var_norm) * gprS.y_train_std(1);
%     else
%         y_mean = y_mean_norm;
%         y_std  = sqrt(y_var_norm);
%     end
% 
%     res_hat = y_mean;
%     res_std = y_std;
% end



%%%%%%%%%%%%%%%%%%%%%%%
% ========== NO STD to speed up inference ==========
function res_hat = residualBladeGPR( ...
    gprS, DEL_NN, DD, YD, gamma_upstream, yaw0_deg, ...
    Hs, Tp, wdir, mean_wind, TI_wind, turb)

    % Feature order must match training:
    % ["DEL_NN", "DD", "YD", "gamma_upstream", "yaw0_deg",
    %  "Hs", "Tp", "wvdir", "mean_wind", "TI_wind", "turb"]

    x = [DEL_NN, DD, YD, gamma_upstream, yaw0_deg, ...
         Hs, Tp, wdir, mean_wind, TI_wind, turb];

    % -------------------------------------------------
    % 1) Standardize input using exported scaler
    % -------------------------------------------------
    x_mean  = gprS.x_mean(:).';
    x_scale = gprS.x_scale(:).';
    x_scaled = (x - x_mean) ./ x_scale;   % 1 x n_features

    % -------------------------------------------------
    % 2) Read exported GPR quantities
    % -------------------------------------------------
    Xtr = gprS.X_train_scaled;          % n_train x n_features
    alpha = gprS.alpha(:);              % n_train x 1

    const_value  = gprS.const_value(1);
    length_scale = gprS.length_scale(:).';   % 1 x n_features

    % -------------------------------------------------
    % 3) Cross-kernel k(x, X_train)
    % Kernel = ConstantKernel * RBF + WhiteKernel
    % WhiteKernel contributes only on diagonal, not here
    % -------------------------------------------------
    DX = (Xtr - x_scaled) ./ length_scale;      % n_train x n_features
    sqdist = sum(DX.^2, 2);                     % n_train x 1
    k_trans = const_value * exp(-0.5 * sqdist); % n_train x 1

    % -------------------------------------------------
    % 4) Predict mean in normalized y-space
    % -------------------------------------------------
    y_mean_norm = k_trans.' * alpha;

    % -------------------------------------------------
    % 5) Undo normalize_y if used in sklearn
    % -------------------------------------------------
    if isfield(gprS, 'normalize_y') && gprS.normalize_y(1) == 1
        y_mean = y_mean_norm * gprS.y_train_std(1) + gprS.y_train_mean(1);
    else
        y_mean = y_mean_norm;
    end

    res_hat = y_mean;
end