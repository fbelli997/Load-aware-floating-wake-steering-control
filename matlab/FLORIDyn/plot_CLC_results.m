function plot_CLC_results(time_hist, time_hist_control, ws_hist, wd_hist, ws_std_hist, wd_std_hist, ...
                           yaw_hist, DEL_NN_hist_blade, DEL_NN_hist_tower, ...
                           DD_hist, YD_hist, TI_hist, nds_hist, T, Sim, DEL_blade_bl_hist_full, DEL_tower_bl_hist_full)

nCols = 3;
nRows = ceil(T.nT/nCols);

%% Wind estimates
if ~isempty(time_hist)
    figure('Position',[100 100 900 700]);
    subplot(2,1,1); hold on; grid on; box on;
    fill([time_hist; flipud(time_hist)], ...
         [ws_hist+ws_std_hist; flipud(ws_hist-ws_std_hist)], ...
         [0.3010 0.7450 0.9330], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(time_hist, ws_hist, 'b-', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Wind speed (m/s)');
    title('Estimated wind speed history');
    legend({'WS \pm1\sigma','WS mean'}, 'Location', 'best');

    subplot(2,1,2); hold on; grid on; box on;
    fill([time_hist; flipud(time_hist)], ...
         [wd_hist+wd_std_hist; flipud(wd_hist-wd_std_hist)], ...
         [0.8500 0.3250 0.0980], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(time_hist, wd_hist, 'r-', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Wind direction (deg)');
    title('Estimated wind direction history');
    legend({'WD \pm1\sigma','WD mean'}, 'Location', 'best');
    sgtitle('Stored EnKF scalar wind estimates');

    save_folder = fullfile(Sim.PathToSim,'Results','EnKF_WindEstimateHistory');
    if ~exist(save_folder, 'dir'); mkdir(save_folder); end
    saveas(gcf, fullfile(save_folder,'WindEstimateHistory.png'));
end

%% Helper per i subplot ripetuti
plotVar(time_hist_control, yaw_hist,          T.nT, nRows, nCols, 'Yaw [deg]',     'Geometric yaw history',  'Yaw history');
plotVar(time_hist_control, DEL_NN_hist_blade, T.nT, nRows, nCols, 'DEL Blade',     'DEL blade history',      'DEL_ History blade');
plotVar(time_hist_control, DEL_NN_hist_tower, T.nT, nRows, nCols, 'DEL Tower',     'DEL tower history',      'DEL_Tower History');
plotVar(time_hist_control, DD_hist,           T.nT, nRows, nCols, 'DD',            'DD history',             'DD History');
plotVar(time_hist_control, YD_hist,           T.nT, nRows, nCols, 'YD',            'YD history',             'YD History');
plotVar(time_hist_control, TI_hist,           T.nT, nRows, nCols, 'TI',            'TI history',             'TI History');
plotVar(time_hist_control, nds_hist,          T.nT, nRows, nCols, 'nds',           'nds history',            'nds');


%% Combined plot: WD + Yaw + DEL blade + DEL tower per turbina
figure('Name', 'Combined_per_turbine', 'Position', [100 100 1400 900]);

nRows_comb = 4;


% --- limiti globali per il combined plot ---
ylims_wd    = [min(wd_hist - wd_std_hist), max(wd_hist + wd_std_hist)];
if diff(ylims_wd)    == 0; ylims_wd    = ylims_wd    + [-1, 1]; end

ylims_yaw   = [min(yaw_hist(:)), max(yaw_hist(:))];
if diff(ylims_yaw)   == 0; ylims_yaw   = ylims_yaw   + [-1, 1]; end

ylims_blade = [min(DEL_NN_hist_blade(:)), max(DEL_NN_hist_blade(:))];
if diff(ylims_blade) == 0; ylims_blade = ylims_blade + [-1, 1]; end

ylims_tower = [min(DEL_NN_hist_tower(:)), max(DEL_NN_hist_tower(:))];
if diff(ylims_tower) == 0; ylims_tower = ylims_tower + [-1, 1]; end

for iT = 1:T.nT
    % Wind direction
    subplot(nRows_comb, T.nT, iT);
    hold on; grid on; box on;
    fill([time_hist; flipud(time_hist)], ...
         [wd_hist+wd_std_hist; flipud(wd_hist-wd_std_hist)], ...
         [0.8500 0.3250 0.0980], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(time_hist, wd_hist, 'r-', 'LineWidth', 1.5);
    xlim([800, 4600]); ylim(ylims_wd);
    if iT == 1
        ylabel('WD [deg]');
        legend({'WD \pm1\sigma','WD mean'}, 'Location','best');
    end
    title(sprintf('Turbine %d', iT));
    set(gca, 'FontSize', 20);
    set(gca, 'XTickLabel', []);

    % Yaw history
    subplot(nRows_comb, T.nT, T.nT + iT);
    plot(time_hist_control, yaw_hist(:,iT), 'b-', 'LineWidth', 1.5);
    if iT == 1; ylabel('Yaw [deg]'); end
    grid on; box on;
    xlim([800, 4600]); ylim(ylims_yaw);
    set(gca, 'FontSize', 20);
    set(gca, 'XTickLabel', []);

    % DEL blade
    subplot(nRows_comb, T.nT, 2*T.nT + iT);
    plot(time_hist_control, DEL_NN_hist_blade(:,iT), 'b-', 'LineWidth', 1.5);
    if iT == 1; ylabel('DEL Blade'); end
    grid on; box on;
    xlim([800, 4600]); ylim(ylims_blade);
    set(gca, 'FontSize', 20);
    set(gca, 'XTickLabel', []);

    % DEL tower
    subplot(nRows_comb, T.nT, 3*T.nT + iT);
    plot(time_hist_control, DEL_NN_hist_tower(:,iT), 'b-', 'LineWidth', 1.5);
    if iT == 1; ylabel('DEL Tower'); end
    grid on; box on;
    xlim([800, 4600]); ylim(ylims_tower);
    xlabel('Time [s]');
    set(gca, 'FontSize', 20);
end

sgtitle('Per-turbine: WD, Yaw, DEL Blade, DEL Tower');



%% Combined plot: DEL baseline vs yaw per turbina
figure('Name', 'DEL_baseline_vs_yaw', 'Position', [100 100 1400 600]);

ylims_bl_blade = [min([DEL_blade_bl_hist_full(:); DEL_NN_hist_blade(:)]), ...
                  max([DEL_blade_bl_hist_full(:); DEL_NN_hist_blade(:)])];
ylims_bl_tower = [min([DEL_tower_bl_hist_full(:); DEL_NN_hist_tower(:)]), ...
                  max([DEL_tower_bl_hist_full(:); DEL_NN_hist_tower(:)])];
if diff(ylims_bl_blade) == 0; ylims_bl_blade = ylims_bl_blade + [-1, 1]; end
if diff(ylims_bl_tower) == 0; ylims_bl_tower = ylims_bl_tower + [-1, 1]; end

for iT = 1:T.nT
    % DEL blade
    subplot(2, T.nT, iT);
    hold on; grid on; box on;
    plot(time_hist_control, DEL_blade_bl_hist_full(:,iT), 'k-', 'LineWidth', 1.5);
    plot(time_hist_control, DEL_NN_hist_blade(:,iT),      'r-',  'LineWidth', 1.5);
    if iT == 1
        ylabel('DEL Blade');
        legend({'Baseline (yaw=0)', 'With yaw'}, 'Location', 'best');
    end
    title(sprintf('Turbine %d', iT));
    xlim([800, 4600]); ylim(ylims_bl_blade);
    set(gca, 'FontSize', 20);
    set(gca, 'XTickLabel', []);

    % DEL tower
    subplot(2, T.nT, T.nT + iT);
    hold on; grid on; box on;
    plot(time_hist_control, DEL_tower_bl_hist_full(:,iT), 'k-', 'LineWidth', 1.5);
    plot(time_hist_control, DEL_NN_hist_tower(:,iT),      'r-',  'LineWidth', 1.5);
    if iT == 1
        ylabel('DEL Tower');
        legend({'Baseline (yaw=0)', 'With yaw'}, 'Location', 'best');
    end
    xlabel('Time [s]');
    xlim([800, 4600]); ylim(ylims_bl_tower);
    set(gca, 'FontSize', 20);
end

sgtitle('DEL Blade & Tower: Baseline vs Yaw optimized');




end

%% --- helper locale ---
function plotVar(t, data, nT, nRows, nCols, ylbl, stitle, figname)
    figure('Name', figname);
    ylims = [min(data(:)), max(data(:))];
    if diff(ylims) == 0; ylims = ylims + [-1, 1]; end  % evita crash se costante
    for iT = 1:nT
        subplot(nRows, nCols, iT)
        plot(t, data(:,iT), 'LineWidth', 1.5)
        grid on
        xlabel('Time [s]'); ylabel(ylbl)
        xlim([800, 4600]); ylim(ylims)
        title(sprintf('Turbine %d', iT))
        set(gca, 'FontSize', 20)
    end
    sgtitle(stitle)
end