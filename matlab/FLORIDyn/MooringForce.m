function [Fmoor, Fline, dist] = MooringForce(x, y, anchor_radius)
%MOORINGFORCE 2D mooring force (surge/sway) with 3 lines.
%
% Inputs
%   x, y          : platform center position (can be scalar or nT x 1)
%   anchor_radius : neutral mooring length (m) used to place anchors (scalar)
%   LUT_dist      : vector of distances (m) for LUT
%   LUT_force     : vector of tensions/force magnitudes (N) for LUT
%
% Outputs
%   Fmoor : (nT x 2) total mooring force [Fx, Fy]
%   Fline : (nT x 2 x 3) per-line force vectors
%   dist  : (nT x 3) per-line distances

    mooringIntp = MooringSpline();


    % --- Make column vectors ---
    x = x(:); 
    y = y(:);
    nT = numel(x);

    % --- Geometry constants (from OFF code) ---
    mooring_len_neutral = anchor_radius;   
    r_fairlead_2_platform_center = 40.868; 

    % Anchors in global frame (3 anchors at 0, +60, -60 deg), radius = neutral length
    acr1 = [-mooring_len_neutral, 0];
    acr2 = [ mooring_len_neutral*cos(pi/3),  mooring_len_neutral*sin(pi/3)];
    acr3 = [ mooring_len_neutral*cos(pi/3), -mooring_len_neutral*sin(pi/3)];

    % Fairleads in platform frame (relative to platform center)

    fl1_rel = [-r_fairlead_2_platform_center, 0];
    fl2_rel = [ r_fairlead_2_platform_center*cos(pi/3),  r_fairlead_2_platform_center*sin(pi/3)];
    fl3_rel = [ r_fairlead_2_platform_center*cos(pi/3), -r_fairlead_2_platform_center*sin(pi/3)];


    % Fairleads in global frame for each turbine/platform
    fairlead1 = [x + fl1_rel(1), y + fl1_rel(2)];
    fairlead2 = [x + fl2_rel(1), y + fl2_rel(2)];
    fairlead3 = [x + fl3_rel(1), y + fl3_rel(2)];

    % Vectors from fairlead -> anchor (same sign convention as OFF: anchor - fairlead)
    vec1 = [acr1(1) - fairlead1(:,1), acr1(2) - fairlead1(:,2)];
    vec2 = [acr2(1) - fairlead2(:,1), acr2(2) - fairlead2(:,2)];
    vec3 = [acr3(1) - fairlead3(:,1), acr3(2) - fairlead3(:,2)];

    % Distances
    dist1 = sqrt(vec1(:,1).^2 + vec1(:,2).^2);
    dist2 = sqrt(vec2(:,1).^2 + vec2(:,2).^2);
    dist3 = sqrt(vec3(:,1).^2 + vec3(:,2).^2);

    dist = [dist1, dist2, dist3];

    % --- Tension magnitudes from LUT (linear interp) ---
    % Equivalent of mooringIntp(dist) in OFF.
    % 'extrap' mimics extrapolation; you can clamp instead if you prefer.
    T1 = mooringIntp(dist1);
    T2 = mooringIntp(dist2);
    T3 = mooringIntp(dist3);

    % --- Unit vectors (avoid division by zero) ---
    epsL = 1e-6;
    e1 = vec1 ./ max(dist1, epsL);
    e2 = vec2 ./ max(dist2, epsL);
    e3 = vec3 ./ max(dist3, epsL);

    % --- Line forces ---
    force1 = T1 .* e1;     % (nT x 2)
    force2 = T2 .* e2;
    force3 = T3 .* e3;

    % Total
    Fmoor = force1 + force2 + force3;

    % Optional per-line output
    Fline = zeros(nT,2,3);
    Fline(:,:,1) = force1;
    Fline(:,:,2) = force2;
    Fline(:,:,3) = force3;
end
