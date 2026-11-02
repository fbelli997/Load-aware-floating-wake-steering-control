%% run_plot_vertical_wake_deflection.m
% Script da eseguire (tasto Run o F5). Disegna la deflessione verticale della
% scia a partire dai risultati gia' salvati. Non modifica nessuna pipeline:
% legge solo EnKF_measurements.csv e salva in una cartella NUOVA "Wake_plots".
%
% Modifica solo le impostazioni qui sotto.

clear; clc; close all;

%% ===== IMPOSTAZIONI =====================================================
% Cartella della simulazione (prende il Case* piu' recente che ha il CSV)
% oppure cartella di un Case preciso, es.:
%   '...\3T_IEA15_floating_00\Results\Case2026_10_20_11_34_22'
src = 'D:\zero_model\FLORIDyn_Matlab-main_3_floating\FLORIDyn_Matlab-main_3\Simulations\2026_3T_DEL\3T_IEA15_floating_00';

variable  = 'Z';      % 'Z' = deflessione verticale, 'Y' = laterale
turbines  = [];       % [] = tutte; oppure es. 2:3 (indici da 1)
timeRange = [];       % [] = tutto il run; oppure es. [1200 4200]
sharedY   = true;     % stessa scala y su tutti i pannelli
%% =========================================================================

addpath(fileparts(mfilename('fullpath')));   % cosi' trova la funzione ovunque tu sia

out = plot_vertical_wake_deflection(src, ...
    'Variable', variable, 'Turbines', turbines, ...
    'TimeRange', timeRange, 'SharedY', sharedY);

disp('Fatto. File salvati:');
disp(out.files');
