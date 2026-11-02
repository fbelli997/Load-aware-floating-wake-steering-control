%% load_evaluation.m
% Check manuale del modello surrogato dei DEL (blade + tower) usato dal
% controller in FLORIDyn_step_*.m: NN predittore (DEL "grezzo" da
% wind/wave/yaw) + NN correttore (residuo appreso su FAST.Farm, che tiene
% conto della posizione nella wake).
%
% USO: modifica i valori nel blocco "INPUT" qui sotto ed esegui lo
% script (F5 / Run). I risultati vengono stampati e salvati nella
% struct "DEL" nel workspace.
%
% NOTE SULLE FEATURE (per non confondersi con FLORIDyn_step_*.m):
%   - "wdir" dentro predictBladeNN/predictTowerNN/residualBladeNN/
%     residualTowerNN e' in realta' la WAVE DIRECTION (wvdir), non la
%     wind direction! Il nome e' fuorviante ma e' cosi' nel training.
%   - la wind direction vera e propria entra come feature separata,
%     trasformata: wind_dir_training = 270 - WindDir, e viene anche
%     scomposta in sin/cos per il correttore.
%   - YawUp viene usato sia come "gamma_upstream" sia come "yaw0_deg":
%     nel controller reale (FLORIDyn_step_15001.m righe 400-406) sono
%     sempre lo stesso valore, quindi qui li teniamo uguali di default.
%   - DD/YD di default (8 e 1) corrispondono a "nessuna turbina a monte"
%     (free stream), vedi LayoutDisp.m.
%   - Turb = posizione della turbina nella wake: 1 = prima fila
%     (free stream), 2 = seconda fila, 3 = terza fila, ecc.
%   - BLADE: qui il secondo stadio (correttore NN residuo) e' DISATTIVATO
%     di proposito -> DEL.blade = solo output del primo stadio
%     (predictBladeNN), senza sottrarre il residuo. La TOWER usa invece
%     entrambi gli stadi, come nel controller reale.

clear; 

%% ===== INPUT: MODIFICA QUI I VALORI =====
WS      = 6.35;   % wind speed effettiva alla turbina [m/s]      (mean_wind)
TI      = 0.108;   % turbulence intensity effettiva [-]           (TI_wind)
Hs      = 3;    % wave height significativa [m]
Tp      = 8;    % wave peak period [s]
Yaw     = 0;    % yaw della turbina in esame [deg]
YawUp   = 0;    % yaw della turbina immediatamente a monte [deg] (gamma_upstream / yaw0_deg)
DD      = 5;    % distanza longitudinale normalizzata (in D) dalla turbina a monte piu' vicina
YD      = -0.6;    % scostamento laterale normalizzato (in D) dalla turbina a monte
WindDir = 270;  % direzione del vento [deg], convenzione SOWFA/FLORIS (come wd_hat)
WaveDir = 0.0;    % direzione delle onde [deg]                   (wvdir)
Turb    = 3;      % posizione nella wake: 1=prima fila, 2=seconda fila, ...

%% ===== CARTELLA DEL CASO (dove stanno i .mat dei pesi e le funzioni NN) =====
caseDir = 'D:\zero_model\FLORIDyn_Matlab-main_3_floating\FLORIDyn_Matlab-main_3\Simulations\2026_3T_DEL\3T_IEA15_floating_00';

%% ===== CARICA I MODELLI SURROGATI =====
here = pwd;
cd(caseDir);
cleanupObj = onCleanup(@() cd(here));

NN_predictor_blade = load('del_blade_uncertainty_grouped_wind_dir.mat');
NN_predictor_tower = load('del_tower_uncertainty_grouped_wind_dir.mat');
NN_corrector_blade = load('nn_residual_blade_wind_dir.mat');
NN_corrector_tower = load('nn_residual_tower_wind_dir.mat');

%% ===== FEATURE DERIVATE (come in FLORIDyn_step_*.m) =====
wind_dir_training = 270 - WindDir;
sin_wd = sind(wind_dir_training);
cos_wd = cosd(wind_dir_training);

%% ===== PREDITTORE NN (DEL "grezzo" da yaw/onde/vento) =====
[DEL_NN_blade, sigma_blade, ~] = predictBladeNN( ...
    NN_predictor_blade, Yaw, Hs, Tp, WaveDir, WS, TI, wind_dir_training);

[DEL_NN_tower, sigma_tower, ~] = predictTowerNN( ...
    NN_predictor_tower, Yaw, Hs, Tp, WaveDir, WS, TI, wind_dir_training);

%% ===== CORRETTORE NN (residuo, dipende da wake/layout) =====
% BLADE: secondo stadio disattivato (vedi nota in testa al file) --
% NN_corrector_blade resta caricato ma non viene usato.
res_tower = residualTowerNN( ...
    NN_corrector_tower, DEL_NN_tower, DD, YD, YawUp, YawUp, ...
    Hs, Tp, WaveDir, WS, TI, Turb, sin_wd, cos_wd);

%% ===== DEL FINALE =====
DEL.blade          = DEL_NN_blade;            % solo primo stadio (no correttore)
DEL.tower          = DEL_NN_tower - res_tower; % predittore - correzione, come in FLORIDyn_step
DEL.blade_raw_NN   = DEL_NN_blade;
DEL.tower_raw_NN   = DEL_NN_tower;
DEL.blade_sigma_NN = sigma_blade;
DEL.tower_sigma_NN = sigma_tower;
DEL.tower_residual = res_tower;

%% ===== STAMPA RISULTATI =====
fprintf('\n===== load_evaluation: DEL check =====\n');
fprintf('Input:  WS=%.2f m/s | TI=%.3f | Hs=%.2f m | Tp=%.2f s\n', WS, TI, Hs, Tp);
fprintf('        WindDir=%.1f deg | WaveDir=%.1f deg\n', WindDir, WaveDir);
fprintf('        Yaw=%.2f deg | YawUp=%.2f deg | DD=%.2f D | YD=%.2f D | Turb(row)=%d\n', ...
    Yaw, YawUp, DD, YD, Turb);
fprintf('---------------------------------------\n');
fprintf('Blade DEL:  NN_pred=%.4f (+-%.4f)  [solo 1o stadio]  -->  DEL = %.4f\n', ...
    DEL_NN_blade, sigma_blade, DEL.blade);
fprintf('Tower DEL:  NN_pred=%.4f (+-%.4f)  residual=%.4f  -->  DEL = %.4f\n', ...
    DEL_NN_tower, sigma_tower, res_tower, DEL.tower);
fprintf('=======================================\n\n');
