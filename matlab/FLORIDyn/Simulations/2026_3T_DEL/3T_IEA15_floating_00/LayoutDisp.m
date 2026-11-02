function [DD, YD, row_idx] = LayoutDisp(wd, layoutX, layoutY, D)
    k = 0.1;
    layoutX = layoutX(:);
    layoutY = layoutY(:);
    nT = numel(layoutX);

    % Default values: no upstream turbine found
    DD = 8 * ones(nT,1);
    YD = 1 * ones(nT,1);
    idx_next = zeros(nT,1);
    row_idx = ones(nT,1);  % default: front row = 1

    % Rotate layout into wind reference frame
    layout = [layoutX.'; layoutY.'];
    theta = wd - 270;
    mat_rot = [cosd(theta), -sind(theta);
               sind(theta),  cosd(theta)];
    rot_lay = mat_rot * layout;
    x_rot = rot_lay(1,:).';
    y_rot = rot_lay(2,:).';
    R = D/2;

    for i = 1:nT
        DD_candidates = [];
        YD_candidates = [];
        idx_candidates = [];

        for j = 1:nT
            if j == i
                continue
            end

            DDtmp = x_rot(i) - x_rot(j);
            YDtmp = y_rot(i) - y_rot(j);

            if DDtmp > 0
                Rwake = R + k * DDtmp;
                if abs(YDtmp) <= Rwake
                    DD_candidates(end+1) = DDtmp / D;
                    YD_candidates(end+1) = YDtmp / D;
                    idx_candidates(end+1) = j;
                end
            end
        end

        if ~isempty(DD_candidates)
            [~, imin] = min(DD_candidates);
            DD(i) = round(DD_candidates(imin));
            YD(i) = round(YD_candidates(imin), 1);
            idx_next(i) = idx_candidates(imin);
            row_idx(i) = numel(DD_candidates) + 1;  % number of upstream turbines + 1
        end
    end
end