%% run_vertical_deflection_6T.m
% Script da eseguire (tasto Run o F5). Mostra i grafici di deflessione della scia
% delle 6 turbine (caso 6T_IEA15_floating_00, Case2026_04_10_15_26_29).
% Legge solo EnKF_measurements.csv e salva in una cartella NUOVA "Wake_plots":
% non modifica nessuna pipeline.
%
% Figure:
%   1) deflessione VERTICALE (Wake Centerline Z), una turbina per pannello
%   2) deflessione LATERALE  (Wake Centerline Y), una turbina per pannello
%   3) sovrapposizione delle 6 turbine (media dell'ensemble) per Z
% Modifica solo le impostazioni qui sotto.

clear; clc; close all;

%% ===== IMPOSTAZIONI =====================================================
caseDir = ['D:\zero_model\FLORIDyn_Matlab-main_3_floating\FLORIDyn_Matlab-main_3\' ...
           'Simulations\2026_6T_Fede_floating\6T_IEA15_floating_00\Results\' ...
           'Case2026_04_10_15_26_29'];

timeRange = [400 Inf]; % [t0 t1] in s; [] = tutto il run. Qui: parte da 400 s
turbines  = [];     % [] = tutte (0..5); oppure indici da 1, es. [2 3 5 6]
sharedY   = true;   % stessa scala y su tutti i pannelli
plotY     = true;   % disegna anche la deflessione laterale
plotOverlay = true; % disegna la figura con le turbine sovrapposte
%% =========================================================================

addpath(fileparts(mfilename('fullpath')));   % trova la funzione ovunque tu sia

outZ = plot_vertical_wake_deflection(caseDir, 'Variable', 'Z', ...
    'Turbines', turbines, 'TimeRange', timeRange, 'SharedY', sharedY);
files = outZ.files;

if plotY
    outY = plot_vertical_wake_deflection(caseDir, 'Variable', 'Y', ...
        'Turbines', turbines, 'TimeRange', timeRange, 'SharedY', sharedY);
    files = [files, outY.files];
end

if plotOverlay
    fig = figure('Color', 'w', 'Position', [150 150 900 480]);
    ax = axes(fig); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    cols = lines(max(outZ.nT, 7));
    for k = 1:numel(outZ.turbines)
        iT = outZ.turbines(k);
        plot(ax, outZ.time, outZ.mean(:, k), 'LineWidth', 1.3, 'Color', cols(iT, :), ...
             'DisplayName', sprintf('Turbine %d', iT - 1));
    end
    xlabel(ax, 'Time (s)'); ylabel(ax, 'Z (m)');
    title(ax, 'Vertical wake deflection - ensemble mean, all turbines');
    legend(ax, 'Location', 'eastoutside');
    ov = fullfile(caseDir, 'Wake_plots', sprintf('WakeCenterlineZ_%dT_overlay.png', outZ.nT));
    exportgraphics(fig, ov, 'Resolution', 200);
    files = [files, {ov}];
end

disp('Fatto. File salvati:');
disp(files');
