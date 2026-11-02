%% ===== Estimate damping from digitized CSVs: pitch + surge =====
clear; clc; close all;

% -------- Files --------
files = { ...
    'damping_pitch.csv', ...
    'damping_surge.csv' ...
};

% -------- Common settings --------
col_t = 1;              % time column
col_x = 2;              % signal column
nPeriods = 5;            % spacing for log decrement (5-10 ok)
usePeaks = 'positive';   % 'positive' or 'abs'
doPlot = true;

% Peak detection tuning (set [] for auto)
minPeakProm = [];        % e.g. 0.02
minPeakDist = [];        % e.g. 10   (seconds)

% Offset removal tail fraction
tailFrac = 0.10;         % last 10% of samples

results = struct([]);

for i = 1:numel(files)
    csvFile = files{i};

    % ----- Read -----
    opts = detectImportOptions(csvFile);
    T = readtable(csvFile, opts);

    t = T{:, col_t};
    x = T{:, col_x};

    ok = isfinite(t) & isfinite(x);
    t = t(ok); x = x(ok);

    % Sort / unique time
    [t, I] = sort(t(:));
    x = x(I);
    [t, ia] = unique(t, 'stable');
    x = x(ia);

    % ----- Remove offset (important esp. pitch) -----
    N = numel(x);
    tailN = max(10, round(tailFrac*N));
    x_off = mean(x(end-tailN+1:end), 'omitnan');
    x0 = x - x_off;

    % ----- Peak picking signal -----
    if strcmpi(usePeaks,'abs')
        sig = abs(x0);
    else
        sig = x0;
    end

    % Auto-tune if empty
    if isempty(minPeakDist)
        minPeakDist_use = 0.01*(t(end)-t(1));
    else
        minPeakDist_use = minPeakDist;
    end
    if isempty(minPeakProm)
        minPeakProm_use = 0.05*(max(sig)-min(sig));
    else
        minPeakProm_use = minPeakProm;
    end

    [pks, locs] = findpeaks(sig, t, ...
        'MinPeakProminence', minPeakProm_use, ...
        'MinPeakDistance',   minPeakDist_use);

    if strcmpi(usePeaks,'positive')
        keep = pks > 0;
        pks = pks(keep);
        locs = locs(keep);
    end

    if numel(pks) < (nPeriods+1)
        error('Too few peaks in %s (%d). Tune MinPeakProm/Dist or digitize more points.', csvFile, numel(pks));
    end

    % ----- Log decrement: single-pair (1 and 1+n) -----
    delta_pair = (1/nPeriods) * log(pks(1)/pks(1+nPeriods));
    zeta_pair  = delta_pair / sqrt((2*pi)^2 + delta_pair^2);

    % ----- Log decrement: robust mean over all pairs (i and i+n) -----
    deltas = (1/nPeriods) * log(pks(1:end-nPeriods) ./ pks(1+nPeriods:end));
    delta_mean = mean(deltas);
    delta_std  = std(deltas);
    zeta_mean  = delta_mean / sqrt((2*pi)^2 + delta_mean^2);

    % ----- Damped period/frequency from peak times -----
    % Use first few intervals
    nInt = min(numel(locs)-1, 8);
    Td = mean(diff(locs(1:nInt+1)), 'omitnan');
    wd = 2*pi / Td;
    wn = wd / sqrt(1 - zeta_mean^2); % use mean zeta

    % Store
    results(i).file = csvFile;
    results(i).offset = x_off;
    results(i).delta_pair = delta_pair;
    results(i).zeta_pair = zeta_pair;
    results(i).delta_mean = delta_mean;
    results(i).delta_std = delta_std;
    results(i).zeta_mean = zeta_mean;
    results(i).Td = Td;
    results(i).fd = 1/Td;
    results(i).wn = wn;
    results(i).wd = wd;
    results(i).pks = pks;
    results(i).locs = locs;
    results(i).t = t;
    results(i).x = x;
    results(i).x0 = x0;

    % Print
    fprintf('\n=== %s ===\n', csvFile);
    fprintf('Offset removed: %.6g\n', x_off);
    fprintf('delta (pair)   = %.6g   zeta (pair)   = %.6g\n', delta_pair, zeta_pair);
    fprintf('delta (mean)   = %.6g   (std=%.6g)   zeta (mean) = %.6g\n', delta_mean, delta_std, zeta_mean);
    fprintf('Td = %.6g s   fd = %.6g Hz\n', Td, 1/Td);
    fprintf('wn = %.6g rad/s  (wd = %.6g rad/s)\n', wn, wd);

    % Plot
    if doPlot
        figure('Name', csvFile);
        plot(t, x, '-', 'DisplayName','raw'); hold on;
        plot(t, x0, '-', 'DisplayName','offset-removed');
        if strcmpi(usePeaks,'abs')
            % show peaks on abs signal, but mark at their time
            plot(locs, interp1(t, x0, locs, 'linear','extrap'), 'o', 'DisplayName','peaks (times)');
        else
            plot(locs, pks, 'o', 'DisplayName','peaks used');
        end
        grid on; xlabel('t [s]'); ylabel('x');
        legend('Location','best');
        title(strrep(csvFile,'_','\_'));
    end
end
