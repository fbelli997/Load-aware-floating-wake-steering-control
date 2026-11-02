%% floating_features.m
% Script da eseguire (Run / F5). Solo LETTURA dei risultati: non modifica nessuna
% pipeline. Disegna i moti della piattaforma galleggiante (surge e pitch) e ne
% stima le frequenze proprie.
%
%   Fig 1  Serie temporali di PtfmSurge [m] e PtfmPitch [deg] (FAST.Farm, per turbina)
%   Fig 2  Spettro di potenza (Welch) di surge e pitch + frequenze proprie:
%            - linea continua nera : modello lineare 2-DOF usato da FLORIDyn
%                                    (matrici VolturnUS come in DispActiveFloat.m)
%            - linea tratteggiata  : picco misurato nello spettro FAST.Farm
%   Fig 3  Modi di vibrare del modello 2-DOF (surge-pitch): profilo di spostamento
%          lungo la torre, frequenza propria e smorzamento modale
%
% Dati: <fastfarm_dir>\<porta>.T<n>.out  (colonne PtfmSurge, PtfmPitch, dt = 0.1 s)
% Output: figure a schermo + PNG in  <questa cartella>\floating_features_plots\

clear; clc; close all;

%% ===== IMPOSTAZIONI =====================================================
fastfarm_dir = 'D:\zero_model\fastfarm_simulation';
ports        = 15201;            % una o piu' porte (una serie di figure per porta)
turbines     = 1:3;              % turbine da mostrare
t_range      = [1200 4200];      % finestra analizzata [s] (scarta il transitorio)

psd_seg_s    = 500;              % segmento Welch [s]: risoluzione ~1/seg Hz; piu' corto = spettro piu' liscio
psd_overlap  = 0.5;              % sovrapposizione segmenti
band_surge   = [0.003 0.02];     % banda dove cercare il picco di surge [Hz]
band_pitch   = [0.02  0.06];     % banda dove cercare il picco di pitch [Hz]
f_max_plot   = 0.15;             % limite asse frequenze nei grafici PSD [Hz]
z_hub        = 150;              % quota hub sopra SWL [m] (solo per disegnare i modi)
save_png     = true;
%% =========================================================================

out_dir = fullfile(fileparts(mfilename('fullpath')), 'floating_features_plots');
if save_png && ~isfolder(out_dir); mkdir(out_dir); end

%% ---- Modello lineare 2-DOF (surge, pitch) --------------------------------
% Stesse matrici e correzioni di FLORIDynCL/prepareSimulationFloating_main_ZMQ.m
% + DispActiveFloat.m (added mass a freq. infinita, rigidezza idrostatica e
% gravitazionale sul pitch).
M0 = [ 2.01e7  -4.66e7;  -4.66e7  4.26e10];
B0 = [1.6299e5  0;        0       1.2627e9];
K0 = [ 7.192e4  1.145e6;  1.145e6 2.587e8];
A_inf  = [9.64e6 -1.01e8; -1.01e8 1.16e10];
Chs55  = 2.182166e5 * 1025 * 9.8066;
C_grav = 20706630 * 9.8066 * 1.9;
M2 = M0 + A_inf;
K2 = K0;  K2(2,2) = K2(2,2) + Chs55 + C_grav;
[V, Dm] = eig(K2, M2);
[wn2, ix] = sort(sqrt(diag(Dm)));  V = V(:, ix);
fn_model = wn2 / (2*pi);                                   % [Hz]  (1=surge, 2=pitch)
zeta_model = zeros(2,1);
for k = 1:2
    ph = V(:,k);
    zeta_model(k) = (ph'*B0*ph) / (2*wn2(k)*(ph'*M2*ph));  % smorzamento modale
end
fprintf('\n=== Modello lineare 2-DOF (come DispActiveFloat.m) ===\n');
fprintf('  modo 1 (surge): fn = %.4f Hz (T = %.0f s), zeta = %.3f\n', fn_model(1), 1/fn_model(1), zeta_model(1));
fprintf('  modo 2 (pitch): fn = %.4f Hz (T = %.1f s), zeta = %.3f\n', fn_model(2), 1/fn_model(2), zeta_model(2));

%% ---- Fig 3: modi di vibrare del modello ------------------------------------
fig3 = figure('Color','w','Position',[100 100 900 480]);
zz = linspace(-30, z_hub, 60);
modeNames = {'Modo 1 (surge)', 'Modo 2 (pitch)'};
for k = 1:2
    subplot(1,2,k); hold on; grid on; box on;
    ph = V(:,k);                                    % [surge m ; pitch rad]
    xz = ph(1) + zz*ph(2);                          % spostamento orizzontale lungo l'asse
    s  = max(abs(xz));  xz = xz / s;
    plot(xz, zz, 'b-', 'LineWidth', 2);
    plot(0*zz, zz, 'k:');
    yline(0, 'c-', 'SWL');
    if abs(ph(2)) > 1e-12 && abs(-ph(1)/ph(2)) < 400     % nodo di rotazione
        zn = -ph(1)/ph(2);
        plot(0, zn, 'ro', 'MarkerFaceColor', 'r');
        text(0.05, zn, sprintf(' nodo z = %.0f m', zn), 'Color', 'r');
    end
    xlabel('Spostamento orizzontale normalizzato [-]');
    ylabel('Quota rispetto a SWL [m]');
    title(sprintf('%s\nf_n = %.4f Hz (T = %.0f s), \\zeta = %.3f', modeNames{k}, ...
          fn_model(k), 1/fn_model(k), zeta_model(k)));
    xlim([-1.3 1.3]);
end
sgtitle('Modi di vibrare del modello 2-DOF surge-pitch');
if save_png; exportgraphics(fig3, fullfile(out_dir, 'modes_2DOF.png'), 'Resolution', 200); end

%% ---- Dati FAST.Farm: serie temporali e PSD ----------------------------------
colT = lines(max(turbines));
for port = ports(:)'
    fprintf('\n=== Porta %d ===\n', port);
    D = struct('t', {}, 'surge', {}, 'pitch', {});
    for k = 1:numel(turbines)
        f = fullfile(fastfarm_dir, sprintf('%d.T%d.out', port, turbines(k)));
        [t, sg, pt] = read_ptfm(f);
        m = t >= t_range(1) & t <= t_range(2);
        D(k).t = t(m); D(k).surge = sg(m); D(k).pitch = pt(m);
        fprintf('  T%d: %d campioni, t = %.0f..%.0f s, dt = %.3f s\n', turbines(k), nnz(m), t(find(m,1)), t(find(m,1,'last')), median(diff(t)));
    end

    % ---- Fig 1: serie temporali ----
    fig1 = figure('Color','w','Position',[100 100 1000 600]);
    ax1 = subplot(2,1,1); hold on; grid on; box on;
    ax2 = subplot(2,1,2); hold on; grid on; box on;
    for k = 1:numel(turbines)
        c = colT(turbines(k), :);
        plot(ax1, D(k).t, D(k).surge, 'Color', c, 'DisplayName', sprintf('T%d (media %.1f m)', turbines(k), mean(D(k).surge)));
        plot(ax2, D(k).t, D(k).pitch, 'Color', c, 'DisplayName', sprintf('T%d (media %.2f deg)', turbines(k), mean(D(k).pitch)));
    end
    ylabel(ax1, 'Surge [m]');  legend(ax1, 'Location', 'eastoutside');
    ylabel(ax2, 'Pitch [deg]'); xlabel(ax2, 'Time [s]'); legend(ax2, 'Location', 'eastoutside');
    title(ax1, sprintf('Moti della piattaforma, porta %d', port));
    linkaxes([ax1 ax2], 'x'); xlim(ax1, t_range);
    if save_png; exportgraphics(fig1, fullfile(out_dir, sprintf('platform_motion_%d.png', port)), 'Resolution', 200); end

    % ---- Fig 2: PSD + frequenze proprie ----
    fig2 = figure('Color','w','Position',[100 100 1000 650]);
    axS = subplot(2,1,1); hold on; grid on; box on;
    axP = subplot(2,1,2); hold on; grid on; box on;
    fpk = nan(numel(turbines), 2);
    for k = 1:numel(turbines)
        c  = colT(turbines(k), :);
        fs = 1 / median(diff(D(k).t));
        [fq, Ps] = welch_psd(D(k).surge, fs, psd_seg_s, psd_overlap);
        [~,  Pp] = welch_psd(D(k).pitch, fs, psd_seg_s, psd_overlap);
        if k == 1; PsAll = zeros(numel(fq), numel(turbines)); PpAll = PsAll; end
        PsAll(:,k) = Ps;  PpAll(:,k) = Pp;
        semilogy(axS, fq, Ps, 'Color', c, 'LineWidth', 1.0, 'DisplayName', sprintf('T%d', turbines(k)));
        semilogy(axP, fq, Pp, 'Color', c, 'LineWidth', 1.0, 'DisplayName', sprintf('T%d', turbines(k)));
        fpk(k,1) = peak_in_band(fq, Ps, band_surge);
        fpk(k,2) = peak_in_band(fq, Pp, band_pitch);
    end
    % Picco misurato = picco dello spettro MEDIO delle turbine (meno sensibile al rumore)
    PsM = mean(PsAll, 2);  PpM = mean(PpAll, 2);
    fpk_avg = [peak_in_band(fq, PsM, band_surge), peak_in_band(fq, PpM, band_pitch)];
    semilogy(axS, fq, PsM, 'k', 'LineWidth', 2, 'DisplayName', 'media turbine');
    semilogy(axP, fq, PpM, 'k', 'LineWidth', 2, 'DisplayName', 'media turbine');
    decorate(axS, 'PSD surge [m^2/Hz]', fn_model(1), f_max_plot);
    decorate(axP, 'PSD pitch [deg^2/Hz]', fn_model(2), f_max_plot);
    for k = 1:numel(turbines)          % picco piu' prominente di ogni turbina (linea tratteggiata)
        c = colT(turbines(k), :);
        if isfinite(fpk(k,1)); xline(axS, fpk(k,1), '--', 'Color', c, 'LineWidth', 1.5, 'DisplayName', sprintf('picco T%d: %.4f Hz', turbines(k), fpk(k,1))); end
        if isfinite(fpk(k,2)); xline(axP, fpk(k,2), '--', 'Color', c, 'LineWidth', 1.5, 'DisplayName', sprintf('picco T%d: %.4f Hz', turbines(k), fpk(k,2))); end
    end
    legend(axS, 'Location', 'northeast', 'FontSize', 8); legend(axP, 'Location', 'northeast', 'FontSize', 8);
    xlabel(axP, 'Frequenza [Hz]');
    title(axS, sprintf('Spettro dei moti, porta %d  (finestra %d-%d s, segmenti %d s, risoluzione ~%.4f Hz)', ...
          port, t_range(1), t_range(2), psd_seg_s, 1/psd_seg_s));
    if save_png; exportgraphics(fig2, fullfile(out_dir, sprintf('platform_PSD_%d.png', port)), 'Resolution', 200); end

    fprintf('  Frequenze proprie [Hz]      surge     pitch\n');
    fprintf('    modello 2-DOF          %8.4f  %8.4f\n', fn_model(1), fn_model(2));
    fprintf('    picco misurato (media) %8.4f  %8.4f\n', fpk_avg(1), fpk_avg(2));
    for k = 1:numel(turbines)
        fprintf('    picco misurato T%d      %8.4f  %8.4f\n', turbines(k), fpk(k,1), fpk(k,2));
    end
end
if save_png; fprintf('\nPNG salvati in: %s\n', out_dir); end

%% ===== funzioni locali ======================================================
function [t, surge, pitch] = read_ptfm(file)
% Legge Time, PtfmSurge, PtfmPitch da un .out FAST.Farm (header di testo + riga unita').
fid = fopen(file, 'r');
if fid < 0; error('File non trovato: %s', file); end
cleaner = onCleanup(@() fclose(fid));
hdrIdx = 0; names = {};
for i = 1:15
    ln = fgetl(fid);
    if ~ischar(ln); break; end
    if startsWith(strtrim(ln), 'Time')
        hdrIdx = i; names = strsplit(strtrim(ln)); break;
    end
end
if hdrIdx == 0; error('Header "Time" non trovato in %s', file); end
iT = find(strcmp(names, 'Time'), 1);
iS = find(strcmp(names, 'PtfmSurge'), 1);
iP = find(strcmp(names, 'PtfmPitch'), 1);
if isempty(iS) || isempty(iP); error('PtfmSurge/PtfmPitch non presenti in %s', file); end
fgetl(fid);                                              % riga unita'
fmt = repmat('%*f', 1, numel(names));
% legge solo le 3 colonne: sostituisce %*f (skip) con %f nelle posizioni utili
tok = strsplit(strtrim(repmat('%*f ', 1, numel(names))));
tok([iT iS iP]) = {'%f'};
C = textscan(fid, strjoin(tok, ' '), 'CollectOutput', false);
t = C{1}; surge = C{2}; pitch = C{3};
n = min([numel(t) numel(surge) numel(pitch)]);
t = t(1:n); surge = surge(1:n); pitch = pitch(1:n);
end

function [f, P] = welch_psd(x, fs, segLen_s, overlap)
% PSD monolaterale con il metodo di Welch (finestra di Hann, senza toolbox).
x = x(:) - mean(x);
N = numel(x);
L = min(round(segLen_s * fs), N);
step = max(1, round(L * (1 - overlap)));
w = 0.5 - 0.5*cos(2*pi*(0:L-1)'/(L-1));
nfft = 2^nextpow2(4*L);                                  % zero-padding per picchi piu' lisci
starts = 1:step:(N - L + 1);
P = zeros(nfft/2 + 1, 1);
for s = starts
    seg = x(s:s+L-1);
    seg = seg - mean(seg);
    X = fft(seg .* w, nfft);
    P = P + abs(X(1:nfft/2+1)).^2;
end
P = P / numel(starts) / (fs * sum(w.^2));
P(2:end-1) = 2 * P(2:end-1);
f = (0:nfft/2)' * fs / nfft;
end

function fp = peak_in_band(f, P, band)
% Frequenza del picco piu' PROMINENTE (massimo locale interno alla banda), con
% interpolazione parabolica. La prominenza e' il rapporto tra il picco e la piu'
% alta delle due "valli" che lo separano da un picco piu' alto: evita di scegliere
% il bordo della banda quando lo spettro e' in forte discesa (forzante a bassa
% frequenza del vento).
idx = find(f >= band(1) & f <= band(2));
if numel(idx) < 3; fp = NaN; return; end
Pl = P(idx);
loc = find(Pl(2:end-1) > Pl(1:end-2) & Pl(2:end-1) >= Pl(3:end)) + 1;   % massimi locali interni
if isempty(loc); fp = NaN; return; end
prom = zeros(size(loc));
for j = 1:numel(loc)
    i = loc(j);
    li = find(Pl(1:i-1) > Pl(i), 1, 'last');       % primo punto piu' alto a sinistra
    ri = i + find(Pl(i+1:end) > Pl(i), 1, 'first'); % primo punto piu' alto a destra
    lmin = min(Pl(max(li,1):i)); if isempty(li); lmin = min(Pl(1:i)); end
    if isempty(ri); rmin = min(Pl(i:end)); else; rmin = min(Pl(i:ri)); end
    prom(j) = Pl(i) / max(lmin, rmin);
end
[~, jb] = max(prom); k = idx(loc(jb));
if k > 1 && k < numel(f)
    y1 = P(k-1); y2 = P(k); y3 = P(k+1);
    d = 0.5 * (y1 - y3) / (y1 - 2*y2 + y3);
    if isfinite(d) && abs(d) <= 1; fp = f(k) + d*(f(2)-f(1)); return; end
end
fp = f(k);
end

function decorate(ax, ylab, f_model, f_max)
xline(ax, f_model, 'g-', 'LineWidth', 2.5, 'DisplayName', sprintf('modello 2-DOF: %.4f Hz', f_model));
set(ax, 'YScale', 'log'); xlim(ax, [0 f_max]);
ylabel(ax, ylab);
end
