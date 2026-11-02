function out = plot_vertical_wake_deflection(src, varargin)
%PLOT_VERTICAL_WAKE_DEFLECTION  Stand-alone plot of the (vertical) wake deflection.
%
%   Reads the "Wake Centerline Z [m]" column of an already stored
%   EnKF_measurements.csv and draws, for each turbine, the ensemble mean
%   +/- 1 std versus time. It only READS results: nothing in the simulation /
%   control pipelines is touched, and the figures go in a NEW folder.
%
%   The number of turbines nT and of ensemble members nE are detected from the
%   file (rows are stored as: time step -> nT rows (one per turbine); the block
%   is repeated for each ensemble member), so it works for 3T, 6T, ... layouts.
%
%   Note: the stored value is the deflection (deltaZ [m]) that the upstream wake
%   imposes at the turbine, i.e. 0 for a turbine with no upstream wake
%   (e.g. Turbine 0).
%
%   USAGE
%     plot_vertical_wake_deflection(src)
%     plot_vertical_wake_deflection(src, Name, Value, ...)
%
%   src can be
%     - a case folder   ...\Results\Case2026_10_20_11_34_22
%     - a simulation folder   ...\3T_IEA15_floating_00   (newest Case* with a
%                              EnKF_Measurements\EnKF_measurements.csv is used)
%     - the path of an EnKF_measurements.csv
%     - a table, e.g. from a live workspace:  cat(1, EnKF.M{:})
%
%   NAME-VALUE OPTIONS
%     'Variable'  'Z' (default, vertical) or 'Y' (lateral) wake centerline
%     'Turbines'  turbine indices to plot, 1-based (default: all)
%     'TimeRange' [t0 t1] in s (default: full run)
%     'SharedY'   same y-axis on every panel (default true)
%     'nT'        force number of turbines (default: auto-detected)
%     'OutDir'    output folder (default: <case folder>\Wake_plots; for a table
%                 input nothing is saved unless OutDir is given)
%     'Show'      keep the figure open/visible (default true)
%     'FontSize'  axes/label font size (default: MATLAB default)
%     'Layout'    [rows cols] of the panels (default: automatic)
%     'Position'  figure position [x y w h] in pixels (default: automatic)
%     'SharedLabels' one common x/y label for the whole figure instead of one per panel
%     'FirstIndex' number shown for the first turbine (default 0 -> Turbine 0; use 1 for T1..Tn)
%     'ShowTitle' draw the overall figure title (default true); the per-panel
%                 titles (Turbine n) are always drawn
%
%   EXAMPLES
%     plot_vertical_wake_deflection('D:\...\3T_IEA15_floating_00\Results\Case2026_10_20_11_34_22')
%     plot_vertical_wake_deflection('D:\...\3T_IEA15_floating_00')          % newest case
%     plot_vertical_wake_deflection(csvPath, 'Turbines', 2:3, 'TimeRange', [1200 4200])
%     plot_vertical_wake_deflection(cat(1, EnKF.M{:}))                        % from workspace
%
%   OUTPUT (struct): time, mean, std, nT, nE, csv, files.

p = inputParser;
addParameter(p, 'Variable',  'Z');
addParameter(p, 'Turbines',  []);
addParameter(p, 'TimeRange', []);
addParameter(p, 'SharedY',   true);
addParameter(p, 'nT',        []);
addParameter(p, 'OutDir',    '');
addParameter(p, 'Show',      true);
addParameter(p, 'FontSize',  []);
addParameter(p, 'ShowTitle', true);
addParameter(p, 'FirstIndex', 0);
addParameter(p, 'Layout', []);
addParameter(p, 'Position', []);
addParameter(p, 'SharedLabels', false);
parse(p, varargin{:});
opt = p.Results;
axisName = upper(char(opt.Variable));
if ~ismember(axisName, {'Y', 'Z'})
    error('Variable must be ''Y'' or ''Z''.');
end

%% ---- locate the data ----------------------------------------------------
csvPath = '';
caseDir = '';
if istable(src)
    M = src;
else
    src = char(src);
    if isfolder(src)
        csvPath = find_measurements_csv(src);
    elseif isfile(src)
        csvPath = src;
    else
        error('Not a folder / file: %s', src);
    end
    fprintf('Reading %s\n', csvPath);
    M = readtable(csvPath, 'VariableNamingRule', 'preserve');
    caseDir = fileparts(fileparts(csvPath));          % ...\Case...\
end

vn = M.Properties.VariableNames;
timeCol = find(strncmpi(vn, 'Time', 4), 1);
valCol  = find(~cellfun('isempty', regexpi(vn, ['wake.*centerline\s+' axisName '(\s|$)'])), 1);
if isempty(timeCol) || isempty(valCol)
    error('Columns "Time [s]" / "Wake Centerline %s [m]" not found in the table.', axisName);
end
t_all = M{:, timeCol};
v_all = M{:, valCol};
R = numel(t_all);

%% ---- detect nT and nE ---------------------------------------------------
if isempty(opt.nT)
    nT = find(t_all ~= t_all(1), 1) - 1;             % first run of equal times
    if isempty(nT); nT = R; end
else
    nT = opt.nT;
end
nE = sum(t_all == t_all(1)) / nT;                    % blocks restart at t_all(1)
if nE ~= round(nE) || mod(R, nT * nE) ~= 0
    error(['Cannot split %d rows with nT=%d (nE=%.2f). Pass the right number of ' ...
           'turbines with ''nT'', n.'], R, nT, nE);
end
nPerE  = R / nE;                                     % rows per ensemble member
nSteps = nPerE / nT;
fprintf('Detected nT = %d turbines, nE = %d ensemble members, %d time steps.\n', nT, nE, nSteps);

V = reshape(v_all, nPerE, nE);                       % [nPerE x nE]
timeVec = t_all(1:nT:nPerE);

%% ---- selection ----------------------------------------------------------
turb = opt.Turbines;
if isempty(turb); turb = 1:nT; end
if any(turb < 1 | turb > nT); error('Turbines must be within 1..%d.', nT); end
tmask = true(size(timeVec));
if ~isempty(opt.TimeRange)
    tmask = timeVec >= opt.TimeRange(1) & timeVec <= opt.TimeRange(2);
end

mu = nan(nnz(tmask), numel(turb));
sd = nan(nnz(tmask), numel(turb));
for k = 1:numel(turb)
    rows = turb(k):nT:nPerE;
    Vt = V(rows(tmask), :);
    mu(:, k) = mean(Vt, 2, 'omitnan');
    sd(:, k) = std(Vt, 0, 2, 'omitnan');
end
tt = timeVec(tmask);

%% ---- plot ---------------------------------------------------------------
nP    = numel(turb);
if nP <= 3; nCols = nP; else; nCols = ceil(nP / 2); end
nRows = ceil(nP / nCols);
if ~isempty(opt.Layout); nRows = opt.Layout(1); nCols = opt.Layout(2); end
vis = 'off'; if opt.Show; vis = 'on'; end
figPos = [100 100 max(420*nCols, 700) 300*nRows + 60];
if ~isempty(opt.Position); figPos = opt.Position; end
fig = figure('Visible', vis, 'Color', 'w', 'Position', figPos);
colors = lines(max(nT, 7));
if opt.SharedY
    lo = min(mu(:) - sd(:), [], 'omitnan');
    hi = max(mu(:) + sd(:), [], 'omitnan');
    pad = max(0.1 * (hi - lo), 0.5);
    yl = [lo - pad, hi + pad];
end
txtFS = 12; if ~isempty(opt.FontSize); txtFS = opt.FontSize; end
tl = tiledlayout(fig, nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:nP
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    c = colors(turb(k), :);
    fill(ax, [tt; flipud(tt)], [mu(:, k) + sd(:, k); flipud(mu(:, k) - sd(:, k))], ...
         c, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(ax, tt, mu(:, k), '-', 'Color', c, 'LineWidth', 1.5);
    yline(ax, 0, ':k');
    if ~opt.SharedLabels
        xlabel(ax, 'Time (s)');
        ylabel(ax, sprintf('%s (m)', axisName));
    end
    title(ax, sprintf('Turbine %d', turb(k) - 1 + opt.FirstIndex), 'Color', c);
    if ~isempty(opt.FontSize); ax.FontSize = opt.FontSize; end
    xlim(ax, [tt(1) tt(end)]);
    if opt.SharedY; ylim(ax, yl); end
end
if opt.SharedLabels
    xlabel(tl, 'Time (s)', 'FontSize', txtFS);
    ylabel(tl, sprintf('%s (m)', axisName), 'FontSize', txtFS);
end
if axisName == 'Z'
    ttl = 'Vertical wake deflection (Wake Centerline Z)';
else
    ttl = 'Lateral wake deflection (Wake Centerline Y)';
end
if opt.ShowTitle
    title(tl, sprintf('%s  -  mean \\pm 1\\sigma over %d ensemble members', ttl, nE));
end

%% ---- save (new folder only) --------------------------------------------
files = {};
outDir = opt.OutDir;
if isempty(outDir) && ~isempty(caseDir)
    outDir = fullfile(caseDir, 'Wake_plots');
end
if ~isempty(outDir)
    if ~isfolder(outDir); mkdir(outDir); end
    base = fullfile(outDir, sprintf('WakeCenterline%s_%dT', axisName, nT));
    exportgraphics(fig, [base '.png'], 'Resolution', 200);
    exportgraphics(fig, [base '.pdf'], 'ContentType', 'vector');
    files = {[base '.png'], [base '.pdf']};
    fprintf('Saved:\n  %s\n  %s\n', files{:});
end
if ~opt.Show; close(fig); end

out = struct('time', tt, 'mean', mu, 'std', sd, 'nT', nT, 'nE', nE, ...
             'turbines', turb, 'csv', csvPath, 'files', {files});
end

% ---------------------------------------------------------------------------
function csvPath = find_measurements_csv(folder)
% folder = a Case folder, a Results folder, or a simulation folder.
rel = fullfile('EnKF_Measurements', 'EnKF_measurements.csv');
cand = {fullfile(folder, rel)};
if isfile(cand{1}); csvPath = cand{1}; return; end
% search Case*/ under folder and folder\Results
hits = dir(fullfile(folder, 'Case*', rel));
hits = [hits; dir(fullfile(folder, 'Results', 'Case*', rel))];
if isempty(hits)
    error('No EnKF_Measurements\\EnKF_measurements.csv found under %s', folder);
end
[~, i] = max([hits.datenum]);
csvPath = fullfile(hits(i).folder, hits(i).name);
end
