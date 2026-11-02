%% load_evaluation_GPR.m
% Come load_evaluation.m, ma con il correttore della TOWER = GPR invece del
% NN (stessa chiamata di FLORIDyn_step_15013.m: residualTowerGPR_wd con il
% GPR caricato come in CLC_init_floating_fun_ZMQ_15013.m). Serve a calcolare
% offline il DEL predetto FFNN+GPR dei baseline (15006/15012, 15106/15112,
% 15206/15212), che sono stati simulati solo con il correttore NN.
%
% USO: modifica i valori nel blocco "INPUT" qui sotto ed esegui lo
% script (F5 / Run). I risultati vengono stampati e salvati nella
% struct "DEL" nel workspace.
%
% NOTE SULLE FEATURE (per non confondersi con FLORIDyn_step_*.m):
%   - "wdir" dentro predictBladeNN/predictTowerNN/residualTowerGPR_wd e'
%     in realta' la WAVE DIRECTION (wvdir), non la wind direction! Il nome
%     e' fuorviante ma e' cosi' nel training.
%   - la wind direction vera e propria entra come feature separata,
%     trasformata: wind_dir_training = 270 - WindDir, e viene anche
%     scomposta in sin/cos per il correttore.
%   - YawUp viene usato sia come "gamma_upstream" sia come "yaw0_deg":
%     nel controller reale sono sempre lo stesso valore, quindi qui li
%     teniamo uguali di default. Nel baseline del controller (STEP 4 di
%     FLORIDyn_step_15013.m) sono entrambi 0, come Yaw.
%   - DD/YD di default (8 e 1) corrispondono a "nessuna turbina a monte"
%     (free stream), vedi LayoutDisp.m.
%   - Turb = posizione della turbina nella wake: 1 = prima fila
%     (free stream), 2 = seconda fila, 3 = terza fila, ecc.
%   - BLADE: nessun correttore (solo primo stadio), identico a
%     load_evaluation.m -> il DEL blade e' lo stesso nelle due versioni.
%   - TOWER: DEL = predictTowerNN - residuo GPR. Per confronto viene
%     calcolato anche il residuo NN (DEL.tower_FFNN = valore di
%     load_evaluation.m con gli stessi input).

clear;

%% ===== INPUT: MODIFICA QUI I VALORI =====
WS      = 6.36;   % wind speed effettiva alla turbina [m/s]      (mean_wind)
TI      = 0.108;  % turbulence intensity effettiva [-]           (TI_wind)
Hs      = 4;    % wave height significativa [m]
Tp      = 10;    % wave peak period [s]
Yaw     = 0;    % yaw della turbina in esame [deg]
YawUp   = 0;    % yaw della turbina immediatamente a monte [deg] (gamma_upstream / yaw0_deg)
DD      = 5;    % distanza longitudinale normalizzata (in D) dalla turbina a monte piu' vicina
YD      = 0.6;    % scostamento laterale normalizzato (in D) dalla turbina a monte
WindDir = 270;  % direzione del vento [deg], convenzione SOWFA/FLORIS (come wd_hat)
WaveDir = 0.0;    % direzione delle onde [deg]                   (wvdir)
Turb    = 3;      % posizione nella wake: 1=prima fila, 2=seconda fila, ...

%% ===== CARTELLA DEL CASO (dove stanno i .mat dei pesi e le funzioni NN/GPR) =====
caseDir = 'D:\zero_model\FLORIDyn_Matlab-main_3_floating\FLORIDyn_Matlab-main_3\Simulations\2026_3T_DEL\3T_IEA15_floating_00';
gprFile = 'D:\zero_model\FLORIDyn_Matlab-main_3_floating\FLORIDyn_Matlab-main_3\Simulations\2026_3T_DEL\pts_joblibs\gpr_residual_tower_wind_dir.mat';

%% ===== CARICA I MODELLI SURROGATI =====
here = pwd;
cd(caseDir);
cleanupObj = onCleanup(@() cd(here));

NN_predictor_blade = load('del_blade_uncertainty_grouped_wind_dir.mat');
NN_predictor_tower = load('del_tower_uncertainty_grouped_wind_dir.mat');
NN_corrector_tower = load('nn_residual_tower_wind_dir.mat');   % solo per confronto
GPR_tower = load(gprFile, ...
    'x_mean','x_scale','X_train_scaled','alpha','L','const_value','length_scale', ...
    'noise_level','normalize_y','y_train_mean','y_train_std');

%% ===== FEATURE DERIVATE (come in FLORIDyn_step_*.m) =====
wind_dir_training = 270 - WindDir;
sin_wd = sind(wind_dir_training);
cos_wd = cosd(wind_dir_training);

%% ===== PREDITTORE NN (DEL "grezzo" da yaw/onde/vento) =====
[DEL_NN_blade, sigma_blade, ~] = predictBladeNN( ...
    NN_predictor_blade, Yaw, Hs, Tp, WaveDir, WS, TI, wind_dir_training);

[DEL_NN_tower, sigma_tower, ~] = predictTowerNN( ...
    NN_predictor_tower, Yaw, Hs, Tp, WaveDir, WS, TI, wind_dir_training);

%% ===== CORRETTORE TOWER: GPR (residuo, dipende da wake/layout) =====
res_tower_GPR = residualTowerGPR_wd( ...
    GPR_tower, DEL_NN_tower, DD, YD, YawUp, YawUp, ...
    Hs, Tp, WaveDir, WS, TI, Turb, sin_wd, cos_wd);

% stesso input con il correttore NN, per confronto con load_evaluation.m
res_tower_NN = residualTowerNN( ...
    NN_corrector_tower, DEL_NN_tower, DD, YD, YawUp, YawUp, ...
    Hs, Tp, WaveDir, WS, TI, Turb, sin_wd, cos_wd);

%% ===== DEL FINALE =====
DEL.blade              = DEL_NN_blade;                 % solo primo stadio (no correttore)
DEL.tower              = DEL_NN_tower - res_tower_GPR; % predittore - correzione GPR, come in FLORIDyn_step_15013
DEL.tower_FFNN         = DEL_NN_tower - res_tower_NN;  % valore FFNN+FFNN (= load_evaluation.m)
DEL.blade_raw_NN       = DEL_NN_blade;
DEL.tower_raw_NN       = DEL_NN_tower;
DEL.blade_sigma_NN     = sigma_blade;
DEL.tower_sigma_NN     = sigma_tower;
DEL.tower_residual_GPR = res_tower_GPR;
DEL.tower_residual_NN  = res_tower_NN;

%% ===== STAMPA RISULTATI =====
fprintf('\n===== load_evaluation_GPR: DEL check =====\n');
fprintf('Input:  WS=%.2f m/s | TI=%.3f | Hs=%.2f m | Tp=%.2f s\n', WS, TI, Hs, Tp);
fprintf('        WindDir=%.1f deg | WaveDir=%.1f deg\n', WindDir, WaveDir);
fprintf('        Yaw=%.2f deg | YawUp=%.2f deg | DD=%.2f D | YD=%.2f D | Turb(row)=%d\n', ...
    Yaw, YawUp, DD, YD, Turb);
fprintf('-------------------------------------------\n');
fprintf('Blade DEL:  NN_pred=%.4f (+-%.4f)  [solo 1o stadio]  -->  DEL = %.4f\n', ...
    DEL_NN_blade, sigma_blade, DEL.blade);
fprintf('Tower DEL:  NN_pred=%.4f (+-%.4f)  residual GPR=%.4f  -->  DEL = %.4f\n', ...
    DEL_NN_tower, sigma_tower, res_tower_GPR, DEL.tower);
fprintf('            (confronto FFNN+FFNN: residual NN=%.4f  -->  DEL = %.4f)\n', ...
    res_tower_NN, DEL.tower_FFNN);
fprintf('===========================================\n\n');
