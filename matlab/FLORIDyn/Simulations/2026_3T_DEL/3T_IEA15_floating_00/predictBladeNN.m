function [mu_real, sigma_real, logvar_std] = predictBladeNN(nnS, yaw, Hs, Tp, wdir, mean_wind, TI_wind, wind_dir)

    % input order must match training:
    % ["yaw", "Hs", "Tp", "wdir", "mean_wind", "TI_wind"]
    x = [yaw, Hs, Tp, wdir, mean_wind, TI_wind, wind_dir];

    x_mean  = nnS.x_scaler_mean(:).';
    x_scale = nnS.x_scaler_scale(:).';
    y_mean  = nnS.y_scaler_mean(1);
    y_scale = nnS.y_scaler_scale(1);

    % standardize input
    x_std = (x - x_mean) ./ x_scale;

    % layer 1: Linear + ReLU
    z1 = x_std * nnS.W1.' + nnS.b1(:).';
    a1 = max(z1, 0);

    % layer 2: Linear + ReLU
    z2 = a1 * nnS.W2.' + nnS.b2(:).';
    a2 = max(z2, 0);

    % heads
    mu_std     = a2 * nnS.Wmu.' + nnS.bmu(:).';
    logvar_std = a2 * nnS.Wlv.' + nnS.blv(:).';

    % clamp logvar like in PyTorch
    logvar_std = min(max(logvar_std, nnS.logvar_min(1)), nnS.logvar_max(1));

    % back to real scale
    mu_real    = mu_std * y_scale + y_mean;
    sigma_real = sqrt(exp(logvar_std)) * y_scale;
end