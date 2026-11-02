%% Plot OpenFAST/FAST.Farm outputs (multi-turbine) - SINGLE FIGURE with 5 SUBPLOTS
clear; close all; clc;

% ------------------ USER SETTINGS ------------------
rootName   = '13000';         % if turbines>1: files are 13000.T1.out, 13000.T2.out, ...
turbines   = 1;               % number of turbines (set 1 if you have a single file "13000.out")
t_start    = 600;             % start time for analysis/plots [s]
xlim_main  = [600 1200];      % x-limits for plots
xlines_ev  = [1600 2600];     % vertical event lines [s] (will be drawn only if inside xlim_main)
fc_wd      = 0.01;            % [Hz] cutoff for filtered wind direction (Butterworth 2nd order)

colors = cool(turbines);

% ------------------ PREALLOC (as empty) ------------------
Ct_all         = nan(1, turbines);
Cp_all         = nan(1, turbines);
GenPwr_mean    = nan(1, turbines);
Wind_mean      = nan(1, turbines);
PtfmPitch_mean = nan(1, turbines);
Yaw_mean       = nan(1, turbines);

Time_matrix        = [];
GenPwr_matrix       = [];
RotSpeed_matrix     = [];
PtfmPitch_matrix    = [];
Thrust_matrix       = [];
Yaw_matrix          = [];
Wind1VelX_matrix    = [];
WindDir_matrix      = [];
WindDirFilt_matrix  = [];

WindDiskVelX_matrix = [];
WindDiskVelY_matrix = [];
WindDiskVelZ_matrix = [];
DiskRSS_matrix      = [];

BlPitch1_matrix     = [];
BlPitch2_matrix     = [];
BlPitch3_matrix     = [];
BlPitchColl_matrix  = [];

% ------------------ LOOP TURBINES ------------------
t_ref = [];
for i = 1:turbines

    % ---- FILE NAME ----
    if turbines == 1
        % single turbine file: "13000.out"
        fileOut = sprintf('%s.out', rootName);
        legName = 'T1';
    else
        % multi turbine files: "13000.T1.out", ...
        fileOut = sprintf('%s.T%d.out', rootName, i);
        legName = sprintf('T%d', i);
    end

    if ~isfile(fileOut)
        error('File not found: %s', fileOut);
    end

    [ChannelsBase, ChanNameBase, ChanUnit] = ReadFASTtext(fileOut); %#ok<ASGLU>

    % --- TIME ---
    Time = getChan(ChannelsBase, ChanNameBase, {'Time'});
    if all(isnan(Time))
        error('Time channel not found in %s', fileOut);
    end

    idx_start = find(Time >= t_start, 1, 'first');
    if isempty(idx_start); idx_start = 1; end

    t = Time(idx_start:end);

    % --- CHANNELS (robust retrieval) ---
    GenPwr    = getChan(ChannelsBase, ChanNameBase, {'GenPwr'});
    RotSpeed  = getChan(ChannelsBase, ChanNameBase, {'RotSpeed','RotSpeed1'});
    PtfmPitch = getChan(ChannelsBase, ChanNameBase, {'PtfmPitch'});
    Yaw       = getChan(ChannelsBase, ChanNameBase, {'YawPzn','NacYaw','Yaw'});
    wind_u    = getChan(ChannelsBase, ChanNameBase, {'Wind1VelX','Wind1VelXg','Wind1VelX_1'});
    wind_dir  = getChan(ChannelsBase, ChanNameBase, {'WindHubAngXY','WindDir','WindDir_1'});

    RtFldCt   = getChan(ChannelsBase, ChanNameBase, {'RtFldCt'});
    RtFldCp   = getChan(ChannelsBase, ChanNameBase, {'RtFldCp'});

    % Wind at disk
    WindDiskVelX = getChan(ChannelsBase, ChanNameBase, {'WindDiskVelX'});
    WindDiskVelY = getChan(ChannelsBase, ChanNameBase, {'WindDiskVelY'});
    WindDiskVelZ = getChan(ChannelsBase, ChanNameBase, {'WindDiskVelZ'});
    DiskRSS      = sqrt(WindDiskVelX.^2 + WindDiskVelY.^2 + WindDiskVelZ.^2);

    % Blade pitch
    BlPitch1 = getChan(ChannelsBase, ChanNameBase, {'BldPitch1'});
    BlPitch2 = getChan(ChannelsBase, ChanNameBase, {'BldPitch2'});
    BlPitch3 = getChan(ChannelsBase, ChanNameBase, {'BldPitch3'});
    BlPitchColl = mean([BlPitch1, BlPitch2, BlPitch3], 2, 'omitnan');

    % Thrust (candidate names; depends on OutList)
    Thrust = getChan(ChannelsBase, ChanNameBase, { ...
        'RotThrust', ...   % commonly available
        'RtAeroFxh', ...   % hub aero force x
        'RtAeroFx', ...
        'RtFldFx', ...
        'TwrFxt' ...       % tower force x (proxy)
    });

    % --- FILTER wind direction (if available) ---
    wd_seg = wind_dir(idx_start:end);
    wind_dir_filtered = wd_seg;
    if ~all(isnan(wd_seg))
        fs = 1 / mean(diff(t));
        if fc_wd < fs/2
            [b, a] = butter(2, fc_wd/(fs/2));
            wind_dir_filtered = filtfilt(b, a, wd_seg);
        end
    end

    % --- ALIGN LENGTHS (for safety) ---
    sigs = { ...
        GenPwr(idx_start:end), RotSpeed(idx_start:end), PtfmPitch(idx_start:end), Thrust(idx_start:end), ...
        Yaw(idx_start:end), wind_u(idx_start:end), wd_seg, wind_dir_filtered, ...
        WindDiskVelX(idx_start:end), WindDiskVelY(idx_start:end), WindDiskVelZ(idx_start:end), DiskRSS(idx_start:end), ...
        BlPitch1(idx_start:end), BlPitch2(idx_start:end), BlPitch3(idx_start:end), BlPitchColl(idx_start:end), ...
        RtFldCt(idx_start:end), RtFldCp(idx_start:end) ...
    };

    if isempty(t_ref)
        t_ref = t;
    end

    L = min([length(t_ref), length(t)]);
    t_ref = t_ref(1:L);
    for k = 1:numel(sigs)
        if length(sigs{k}) > L
            sigs{k} = sigs{k}(1:L);
        elseif length(sigs{k}) < L
            sigs{k} = [sigs{k}; nan(L - length(sigs{k}), 1)];
        end
    end

    % Unpack
    GenPwr_seg       = sigs{1};
    RotSpeed_seg     = sigs{2};
    PtfmPitch_seg    = sigs{3};
    Thrust_seg       = sigs{4};
    Yaw_seg          = sigs{5};
    wind_u_seg       = sigs{6};
    wd_seg           = sigs{7};
    wd_f_seg         = sigs{8};
    WindDiskVelX_seg = sigs{9};
    WindDiskVelY_seg = sigs{10};
    WindDiskVelZ_seg = sigs{11};
    DiskRSS_seg      = sigs{12};
    BlPitch1_seg     = sigs{13};
    BlPitch2_seg     = sigs{14};
    BlPitch3_seg     = sigs{15};
    BlPitchColl_seg  = sigs{16};
    RtFldCt_seg      = sigs{17};
    RtFldCp_seg      = sigs{18};

    % --- STORE MATRICES ---
    Time_matrix           = t_ref;
    GenPwr_matrix         = [GenPwr_matrix,         GenPwr_seg];
    RotSpeed_matrix       = [RotSpeed_matrix,       RotSpeed_seg];
    PtfmPitch_matrix      = [PtfmPitch_matrix,      PtfmPitch_seg];
    Thrust_matrix         = [Thrust_matrix,         Thrust_seg];
    Yaw_matrix            = [Yaw_matrix,            Yaw_seg];
    Wind1VelX_matrix      = [Wind1VelX_matrix,      wind_u_seg];
    WindDir_matrix        = [WindDir_matrix,        wd_seg];
    WindDirFilt_matrix    = [WindDirFilt_matrix,    wd_f_seg];

    WindDiskVelX_matrix   = [WindDiskVelX_matrix,   WindDiskVelX_seg];
    WindDiskVelY_matrix   = [WindDiskVelY_matrix,   WindDiskVelY_seg];
    WindDiskVelZ_matrix   = [WindDiskVelZ_matrix,   WindDiskVelZ_seg];
    DiskRSS_matrix        = [DiskRSS_matrix,        DiskRSS_seg];

    BlPitch1_matrix       = [BlPitch1_matrix,       BlPitch1_seg];
    BlPitch2_matrix       = [BlPitch2_matrix,       BlPitch2_seg];
    BlPitch3_matrix       = [BlPitch3_matrix,       BlPitch3_seg];
    BlPitchColl_matrix    = [BlPitchColl_matrix,    BlPitchColl_seg];

    % --- MEANS ---
    Ct_all(i)          = mean(RtFldCt_seg, 'omitnan');
    Cp_all(i)          = mean(RtFldCp_seg, 'omitnan');
    GenPwr_mean(i)     = mean(GenPwr_seg, 'omitnan');
    Wind_mean(i)       = mean(wind_u_seg, 'omitnan');
    PtfmPitch_mean(i)  = mean(PtfmPitch_seg, 'omitnan');
    Yaw_mean(i)        = mean(Yaw_seg, 'omitnan');

    fprintf('[%s] Loaded. Samples used: %d\n', legName, L);
end

% ------------------ TOTAL POWER ------------------
GenPwr_total = sum(GenPwr_matrix, 2, 'omitnan');

% ------------------ SINGLE FIGURE: 5 SUBPLOTS ------------------
figure(100); clf

plot_all = @(t, M) arrayfun(@(ii) plot(t, M(:,ii), 'Color', colors(ii,:), ...
    'DisplayName', sprintf('T%d', ii)), 1:size(M,2));

% 1) THRUST
ax1 = subplot(5,1,1); hold on; grid on
plot_all(Time_matrix, Thrust_matrix);
ylabel('Thrust');
title('Thrust / Power / Rotor Speed / Blade Pitch / Platform Pitch');
xlim(xlim_main);
set(gca,'XTickLabel',[]);

% 2) POWER (with total)
ax2 = subplot(5,1,2); hold on; grid on
plot_all(Time_matrix, GenPwr_matrix);
plot(Time_matrix, GenPwr_total, 'k-', 'LineWidth', 2, 'DisplayName','Total');
ylabel('GenPwr');
xlim(xlim_main);
set(gca,'XTickLabel',[]);

% 3) ROTOR SPEED
ax3 = subplot(5,1,3); hold on; grid on
plot_all(Time_matrix, RotSpeed_matrix);
ylabel('RotSpeed');
xlim(xlim_main);
set(gca,'XTickLabel',[]);

% 4) BLADE PITCH (collective)
ax4 = subplot(5,1,4); hold on; grid on
plot_all(Time_matrix, BlPitchColl_matrix);
ylabel('BlPitch coll');
xlim(xlim_main);
set(gca,'XTickLabel',[]);

% 5) PLATFORM PITCH
ax5 = subplot(5,1,5); hold on; grid on
plot_all(Time_matrix, PtfmPitch_matrix);
ylabel('PtfmPitch');
xlabel('Time (s)');
xlim(xlim_main);

% Vertical lines (only if within xlim)
for xv = xlines_ev
    if xv >= xlim_main(1) && xv <= xlim_main(2)
        xline(ax1, xv, 'r--', 'LineWidth', 1.2);
        xline(ax2, xv, 'r--', 'LineWidth', 1.2);
        xline(ax3, xv, 'r--', 'LineWidth', 1.2);
        xline(ax4, xv, 'r--', 'LineWidth', 1.2);
        xline(ax5, xv, 'r--', 'LineWidth', 1.2);
    end
end

% one legend on top panel
legend(ax1, 'Location', 'best');

% link x-axes
linkaxes([ax1 ax2 ax3 ax4 ax5], 'x');

% ------------------ SUMMARY PRINT ------------------
fprintf('\n=== Means (from t >= %.1f s) ===\n', t_start);
for i = 1:turbines
    fprintf('T%d: Ct=%.4f, Cp=%.4f, GenPwr=%.3g, Wind=%.3f, PtfmPitch=%.3f, Yaw=%.3f\n', ...
        i, Ct_all(i), Cp_all(i), GenPwr_mean(i), Wind_mean(i), PtfmPitch_mean(i), Yaw_mean(i));
end
fprintf('Total mean power (sum): %.3g\n\n', mean(GenPwr_total,'omitnan'));

if all(isnan(Thrust_matrix(:)))
    warning(['Thrust channel not found (all NaN). ' ...
             'Add the correct thrust channel to OutList (e.g., RotThrust) or adjust candidates list.']);
end

%% ====================== LOCAL FUNCTIONS ======================
function sig = getChan(ChannelsBase, ChanNameBase, candidates)
% Return channel vector if any candidate exists (case-insensitive).
% If not found, returns NaN vector (same length as time).
    sig = nan(size(ChannelsBase,1),1);
    for ii = 1:numel(candidates)
        col = find(strcmpi(ChanNameBase, candidates{ii}), 1);
        if ~isempty(col)
            sig = ChannelsBase(:, col);
            return
        end
    end
end
