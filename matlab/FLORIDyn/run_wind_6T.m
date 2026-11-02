%% run_wind_6T.m
% Script da eseguire (tasto Run o F5). Risultati per il PAPER del caso
% 6T_IEA15_floating_00 / Case2026_04_10_15_26_29. Legge solo file gia' salvati
% (nessuna pipeline viene toccata) e salva PNG + PDF in outDir.
%
% Figure (font grande; nessun titolo generale, solo il titolo 'Turbine n' di ogni subplot):
%   WindDirection_6T       direzione del vento: misurata vs stimata EnKF
%   WindSpeed_6T           velocita' del vento stimata dall'EnKF
%   Thrust_measured_vs_estimated_6T   spinta misurata / stimata
%   Layout_6T              layout delle 6 turbine in coordinate cartesiane
%   WakeCenterlineZ_6T     deflessione verticale della scia (e Y se plotY = true)
%
% Formato (dedotto dai file): dir.txt / vel.txt [306 x 360] = passi da 15 s x
% (6 turbine x 60 OP); stima per turbina = media sui suoi 60 OP, banda = +-1 sigma
% (radice della varianza media). thrust_predicted/corrected [180x1] = 30 membri x
% 6 turbine (si assume turbina = indice piu' veloce). Modifica solo le impostazioni.

clear; clc; close all;

%% ===== IMPOSTAZIONI =====================================================
caseDir = ['D:\zero_model\FLORIDyn_Matlab-main_3_floating\FLORIDyn_Matlab-main_3\' ...
           'Simulations\2026_6T_Fede_floating\6T_IEA15_floating_00\Results\' ...
           'Case2026_04_10_15_26_29'];
outDir  = ['D:\FLORIDyn_Matlab-main_3_floating\FLORIDyn_Matlab-main_3_floating\' ...
           'FLORIDyn_Matlab-main_3\floating_paper_results'];
nT  = 6;   nOP = 60;   nE = 30;
timeRange = [400 Inf];     % [t0 t1] in s; [] = tutto. Qui: taglia i primi 400 s
fs        = 18;          % font size dei grafici (figure a larghezza di COLONNA: ~8-9 pt nel paper)
figPos    = [50 50 700 820];   % [x y w h] px delle figure a pannelli (stampa a larghezza di colonna)
figRows   = 3;  figCols = 2;   % disposizione dei pannelli (T1 T2 / T3 T4 / T5 T6)
fsLayout  = 18;          % font size del layout
lw        = 2;           % spessore linee
plotY     = false;       % true = salva anche la deflessione laterale (Y)

% Layout (da turbineArrayPropertiesIEA15.m del caso): x, y in m; T0..T5
layoutXY = [  300 1860;  1500 1860;  2700 1860; ...
              300  900;  1500  900;  2700  900];
D_rotor  = 240;          % diametro rotore IEA 15 MW [m]
wd_layout = 270;         % direzione del vento di riferimento [deg] (vento verso +x)
%% =========================================================================

addpath(fileparts(mfilename('fullpath')));
if ~isfolder(outDir); mkdir(outDir); end
set(groot, 'defaultAxesFontSize', fs, 'defaultLegendFontSize', fs-2, ...
    'defaultTextFontSize', fs, 'defaultLineLineWidth', lw);

%% ---- lettura -----------------------------------------------------------
D  = readmatrix(fullfile(caseDir, 'EnKF_States', 'dir.txt'));
DV = readmatrix(fullfile(caseDir, 'EnKF_States', 'dir_var.txt'));
V  = readmatrix(fullfile(caseDir, 'EnKF_States', 'vel.txt'));
VV = readmatrix(fullfile(caseDir, 'EnKF_States', 'vel_var.txt'));
if size(D, 2) ~= nT*nOP || size(V, 2) ~= nT*nOP
    error('dir.txt/vel.txt: colonne diverse da %d x %d.', nT, nOP);
end

files = dir(fullfile(caseDir, 'EnKF_Thrust', '*_dir_measured.txt'));
tm = sort(cellfun(@(s) str2double(s(1:7)), {files.name}).');
nS = min([size(D,1), size(V,1), numel(tm)]);
t = tm(1:nS);  D = D(1:nS,:);  DV = DV(1:nS,:);  V = V(1:nS,:);  VV = VV(1:nS,:);
fprintf('Passi: %d (%g - %g s), turbine: %d\n', nS, t(1), t(end), nT);

Dm = nan(nS, nT);  Tm = nan(nS, nT);
Tp = nan(nS, nT);  Tp_lo = Tp;  Tp_hi = Tp;  Tc = Tp;
for i = 1:nS
    stem = fullfile(caseDir, 'EnKF_Thrust', sprintf('%07d_', t(i)));
    Dm(i,:) = reshape(readmatrix([stem 'dir_measured.txt']), 1, []);
    Tm(i,:) = reshape(readmatrix([stem 'thrust_measured.txt']), 1, []);
    P = reshape(readmatrix([stem 'thrust_predicted.txt']), nT, nE);
    C = reshape(readmatrix([stem 'thrust_corrected.txt']), nT, nE);
    Tp(i,:) = mean(P, 2).';  Tp_lo(i,:) = min(P, [], 2).';  Tp_hi(i,:) = max(P, [], 2).';
    Tc(i,:) = mean(C, 2).';
end

% media sui 60 OP di ogni turbina e sigma
De = nan(nS, nT); Ds = De; Ve = De; Vs = De;
for k = 1:nT
    c = (k-1)*nOP + (1:nOP);
    De(:,k) = mean(D(:,c), 2);   Ds(:,k) = sqrt(mean(DV(:,c), 2));
    Ve(:,k) = mean(V(:,c), 2);   Vs(:,k) = sqrt(mean(VV(:,c), 2));
end

mask = true(size(t));
if ~isempty(timeRange); mask = t >= timeRange(1) & t <= timeRange(2); end
tt = t(mask);
cols = lines(nT);
orange = [0.8 0.4 0.1];

%% ---- Figura 1: direzione del vento -------------------------------------
fig1 = figure('Color', 'w', 'Position', figPos);
tl = tiledlayout(fig1, figRows, figCols, 'TileSpacing', 'compact', 'Padding', 'compact');
lo = min([De(mask,:) - Ds(mask,:); Dm(mask,:)], [], 'all', 'omitnan');
hi = max([De(mask,:) + Ds(mask,:); Dm(mask,:)], [], 'all', 'omitnan');
pad = max(0.1*(hi-lo), 1);
for k = 1:nT
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    fill(ax, [tt; flipud(tt)], [De(mask,k)+Ds(mask,k); flipud(De(mask,k)-Ds(mask,k))], ...
         cols(k,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    h1 = plot(ax, tt, De(mask,k), '-', 'Color', cols(k,:), 'LineWidth', lw);
    h2 = plot(ax, tt, Dm(mask,k), 'k--', 'LineWidth', lw-0.5);
    xlim(ax, [tt(1) tt(end)]); ylim(ax, [lo-pad-14 hi+pad]);
    tag(ax, k, cols(k,:), fs, 'se');
    if k == 1; legend(ax, [h1 h2], {'EnKF', 'Measured'}, 'Location', 'southwest', 'FontSize', fs-4); end
end
xlabel(tl, 'Time (s)', 'FontSize', fs); ylabel(tl, 'Wind direction (deg)', 'FontSize', fs);
savefig_paper(fig1, outDir, 'WindDirection_6T');

%% ---- Figura 2: velocita' del vento -------------------------------------
fig2 = figure('Color', 'w', 'Position', figPos);
tl2 = tiledlayout(fig2, figRows, figCols, 'TileSpacing', 'compact', 'Padding', 'compact');
lo = min(Ve(mask,:) - Vs(mask,:), [], 'all');  hi = max(Ve(mask,:) + Vs(mask,:), [], 'all');
pad = 0.1*(hi-lo);
for k = 1:nT
    ax = nexttile(tl2); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    fill(ax, [tt; flipud(tt)], [Ve(mask,k)+Vs(mask,k); flipud(Ve(mask,k)-Vs(mask,k))], ...
         cols(k,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(ax, tt, Ve(mask,k), 'Color', cols(k,:), 'LineWidth', lw);
    xlim(ax, [tt(1) tt(end)]); ylim(ax, [lo-pad hi+pad]);
    tag(ax, k, cols(k,:), fs, 'ne');
end
xlabel(tl2, 'Time (s)', 'FontSize', fs); ylabel(tl2, 'Wind speed (m/s)', 'FontSize', fs);
savefig_paper(fig2, outDir, 'WindSpeed_6T');

%% ---- Figura 3: spinta misurata / stimata -------------------
fig3 = figure('Color', 'w', 'Position', figPos);
tl3 = tiledlayout(fig3, figRows, figCols, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:nT
    ax = nexttile(tl3); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    h2 = plot(ax, tt, Tc(mask,k), '-', 'Color', cols(k,:), 'LineWidth', lw);
    h3 = plot(ax, tt, Tm(mask,k), 'k--', 'LineWidth', lw-0.5);
    xlim(ax, [tt(1) tt(end)]); ylim(ax, [0.5 2.5]);
    tag(ax, k, cols(k,:), fs, 'sw');
    if k == 1; legend(ax, [h2 h3], {'Estimated', 'Measured'}, 'Location', 'southeast', 'FontSize', fs-4); end
end
xlabel(tl3, 'Time (s)', 'FontSize', fs); ylabel(tl3, 'Thrust (MN)', 'FontSize', fs);
savefig_paper(fig3, outDir, 'Thrust_measured_vs_estimated_6T');

%% ---- Figura 4: layout in coordinate cartesiane -------------------------
fig4 = figure('Color', 'w', 'Position', [100 100 700 430]);   % larghezza di colonna
ax = axes(fig4); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
% rotore = segmento di lunghezza D perpendicolare al vento (vento da wd_layout)
dirFlow = [-sind(wd_layout), -cosd(wd_layout)];       % verso del flusso (270 deg = da ovest -> +x)
nrm = [-dirFlow(2), dirFlow(1)];
for k = 1:nT
    p0 = layoutXY(k,:);
    r = [p0 - nrm*D_rotor/2; p0 + nrm*D_rotor/2];
    plot(ax, r(:,1), r(:,2), '-', 'Color', cols(k,:), 'LineWidth', 6);
    plot(ax, p0(1), p0(2), 'o', 'MarkerSize', 9, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    text(ax, p0(1)+80, p0(2)+110, sprintf('T%d', k), 'Color', cols(k,:), ...
         'FontSize', fsLayout, 'FontWeight', 'bold');
end
axis(ax, 'equal');
xlim(ax, [min(layoutXY(:,1)) - 400, max(layoutXY(:,1)) + 400]);
ylim(ax, [min(layoutXY(:,2)) - 400, max(layoutXY(:,2)) + 400]);
xlabel(ax, 'x (m)'); ylabel(ax, 'y (m)'); ax.FontSize = fsLayout;
savefig_paper(fig4, outDir, 'Layout_6T');

%% ---- Figura 5: deflessione verticale della scia (funzione esistente) ---
outZ = plot_vertical_wake_deflection(caseDir, 'Variable', 'Z', 'TimeRange', timeRange, ...
    'OutDir', outDir, 'FontSize', fs, 'Layout', [figRows figCols], 'Position', figPos, 'SharedLabels', true, 'ShowTitle', false, 'Show', false, 'FirstIndex', 1);
if plotY
    plot_vertical_wake_deflection(caseDir, 'Variable', 'Y', 'TimeRange', timeRange, ...
        'OutDir', outDir, 'FontSize', fs, 'Layout', [figRows figCols], 'Position', figPos, 'SharedLabels', true, 'ShowTitle', false, 'Show', false, 'FirstIndex', 1);
end

disp(['Fatto. Risultati in: ' outDir]);
disp(struct2table(dir(outDir)).name);

%% ========================== funzioni locali ==============================
function tag(ax, k, col, fs, ~)
% Nome della turbina come titolo del subplot.
title(ax, sprintf('Turbine %d', k), 'Color', col, 'FontWeight', 'bold', 'FontSize', fs);
end

function savefig_paper(fig, outDir, name)
exportgraphics(fig, fullfile(outDir, [name '.png']), 'Resolution', 300);
exportgraphics(fig, fullfile(outDir, [name '.pdf']), 'ContentType', 'vector');
end
