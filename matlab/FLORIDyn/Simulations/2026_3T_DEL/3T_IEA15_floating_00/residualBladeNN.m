function res_hat = residualBladeNN( ...
    nnS, DEL_NN, DD, YD, gamma_upstream, yaw0_deg, ...
    Hs, Tp, wdir, mean_wind, TI_wind, turb, sin_wind_dir, cos_wind_dir)
% Feature order must match training:
% ["DEL_NN", "DD", "YD", "gamma_upstream", "yaw0_deg",
%  "Hs", "Tp", "wvdir", "mean_wind", "TI_wind", "turb", "sin_wind_dir", "cos_wind_dir"]

    x = [DEL_NN, DD, YD, gamma_upstream, yaw0_deg, ...
         Hs, Tp, wdir, mean_wind, TI_wind, turb, sin_wind_dir, cos_wind_dir];

    % 1) Standardize input
    x_mean  = nnS.x_scaler_mean(:).';
    x_scale = nnS.x_scaler_scale(:).';
    x_std = (x - x_mean) ./ x_scale;

    % 2) Layer 1: Linear + ReLU
    z1 = x_std * nnS.W1.' + nnS.b1(:).';
    a1 = max(z1, 0);

    % 3) Layer 2: Linear + ReLU
    z2 = a1 * nnS.W2.' + nnS.b2(:).';
    a2 = max(z2, 0);

    % 4) Output layer
    y_std = a2 * nnS.Wout.' + nnS.bout(:).';

    % 5) Undo y scaling
    res_hat = y_std * nnS.y_std(1) + nnS.y_mean(1);

end