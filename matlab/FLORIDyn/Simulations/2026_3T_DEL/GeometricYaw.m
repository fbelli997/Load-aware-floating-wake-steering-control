function [yaw_geom_i, fil_wake, upstream_idx] = GeometricYaw(self, x, y)
%GEOMETRICYAW Single-step implementation coherent with Python GeomWFFC logic
% Returns:
%   yaw_geom_i   : geometric yaw command for each turbine [deg]
%   fil_wake     : pairwise wake mask, fil_wake(i,j)=true if turbine i wakes turbine j
%   upstream_idx : closest upstream turbine affecting each turbine (0 if none)

    x = x(:);
    y = y(:);

    ws_eff_i = self.ws_hat(:);
    wd_deg   = self.wd_hat;

    [dx_eq_i, dy_eq_i, fil_wake, dx_ii, ~] = calculate_dxdy_eq_singleStep( ...
        self, x, y, self.diameter, self.ws_rated, wd_deg, ws_eff_i);

    yaw_geom_i = calculate_geometric_yaw(self, dx_eq_i, dy_eq_i);

    % Closest upstream turbine for each turbine
    nT = numel(x);
    upstream_idx = zeros(nT,1);

    for i = 1:nT
        % turbines waking turbine i
        candidates = find(fil_wake(:,i));

        if isempty(candidates)
            upstream_idx(i) = 0;
        else
            dx_candidates = dx_ii(candidates,i);
            [~,kmin] = min(dx_candidates);
            upstream_idx(i) = candidates(kmin);
        end
    end
end


%% ===================================================================== %%
%% Helpers
%% ===================================================================== %%

function [dx_eq_i, dy_eq_i, fil_wake, dx_ii, dy_ii] = ...
    calculate_dxdy_eq_singleStep(self, x, y, diameter, ws_rated, wd_deg, ws_eff_i)

    % Pairwise normalized distances in wake-aligned frame
    [dx_ii, dy_ii] = calculate_dxdy(self, x, y, wd_deg, diameter);

    % Apply effective wind-speed filter
    [dx_fil, dy_fil] = ws_eff_filter(self, dx_ii, dy_ii, ws_eff_i, ws_rated);

    % Wake mask for equivalent farm formulation
    fil_wake = ((dy_fil <= (0.5 + self.k_farmEqFunc .* dx_fil)) & ...
                (dy_fil >= (-0.5 - self.k_farmEqFunc .* dx_fil))) & ...
                isfinite(dx_fil);

    dx_fil(~fil_wake) = inf;
    dy_fil(~fil_wake) = inf;

    % Equivalent downwind distance
    dx_eq_i = min(dx_fil, [], 2);

    % Separate positive and negative lateral offsets
    fil_dy_pos = dy_fil >= 0;
    fil_dy_neg = ~fil_dy_pos;

    dy_pos = dy_fil;
    dy_neg = dy_fil;

    dy_pos(~fil_dy_pos) = inf;
    dy_neg(~fil_dy_neg) = -inf;

    % Delta wrt equivalent dx
    dx_eq_ii = repmat(dx_eq_i, 1, size(dx_fil,2));

    delta_dx_eq = zeros(size(dx_fil));
    delta_dx_eq(fil_wake) = dx_fil(fil_wake) - dx_eq_ii(fil_wake);

    % Equivalent positive / negative dy
    dy_eq_pos_i = min(dy_pos + self.a_farmEqFunc .* delta_dx_eq, [], 2);
    dy_eq_neg_i = -min(abs(dy_neg) + self.a_farmEqFunc .* delta_dx_eq, [], 2);

    % Final equivalent dy
    dy_eq_i = zeros(size(dx_eq_i));
    fil_dominant_side = false(size(dx_eq_i));
    fil_wake_eq = any(fil_wake, 2);

    fil_dominant_side(fil_wake_eq) = ...
        abs(dy_eq_pos_i(fil_wake_eq) - abs(dy_eq_neg_i(fil_wake_eq))) > self.c_farmEqFunc;

    dy_eq_i(fil_dominant_side) = ...
        sign(-(dy_eq_pos_i(fil_dominant_side) - abs(dy_eq_neg_i(fil_dominant_side)))) .* ...
        min(dy_eq_pos_i(fil_dominant_side), abs(dy_eq_neg_i(fil_dominant_side)));
    
    % Safety
    dy_eq_i(~isfinite(dx_eq_i)) = 0;
    
    tol = 1e-3;
    dy_eq_i(abs(dy_eq_i) < tol) = 0;
end


function yaw_geom = calculate_geometric_yaw(self, dx, dy)

    yaw_geom = sign(dy) .* max( ...
        self.yaw_max - self.yaw_m_x .* dx - self.yaw_m_y .* abs(dy), ...
        zeros(size(dx)));

    fil_wake = (dy <= (0.5 + self.k_yaw .* dx)) & ...
               (dy >= (-0.5 - self.k_yaw .* dx));
    yaw_geom(~fil_wake) = 0;

    fil_ub = dy > ( self.yaw_c_center_y + (self.yaw_m_center_x .* dx).^2 );
    fil_lb = dy < (-self.yaw_c_center_y - (self.yaw_m_center_x .* dx).^2);
    yaw_geom((~fil_ub) & (~fil_lb)) = 0;
end


function [dx, dy] = calculate_dxdy(~, x, y, wd_deg, diameter)

    x = x(:);
    y = y(:);

    x_mat_1 = repmat(x, 1, length(x));
    x_mat_2 = repmat(x.', length(x), 1);
    y_mat_1 = repmat(y, 1, length(y));
    y_mat_2 = repmat(y.', length(y), 1);

    d = sqrt((x_mat_1 - x_mat_2).^2 + (y_mat_1 - y_mat_2).^2);
    theta = atan2(y_mat_2 - y_mat_1, x_mat_2 - x_mat_1);
    gamma = wd_deg * (pi/180) - (3/2) * pi + theta;

    dx = d .* cos(gamma) / diameter;
    dy = d .* sin(gamma) / diameter;

    % Keep only downstream interactions
    eps_val = 1e-3;
    dx(dx < eps_val) = inf;
end


function [dx_all_fil, dy_all_fil] = ws_eff_filter(~, dx_all, dy_all, ws_eff_i, ws_rated)

    ws_eff_i = ws_eff_i(:);
    ws_eff_ii = repmat(ws_eff_i, 1, length(ws_eff_i));

    % If turbine i sees rated or above, exclude it from steering logic
    fil_ws_eff = ws_eff_ii >= ws_rated;

    dx_all_fil = dx_all;
    dy_all_fil = dy_all;

    dx_all_fil(fil_ws_eff) = inf;
end