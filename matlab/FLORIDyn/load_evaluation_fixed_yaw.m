%% load_evaluation_fixed_yaw.m
% Confronto DEL FAST.Farm (MLife) vs modello surrogato per le simulazioni a
% yaw fisso SENZA controllo attivo: famiglie 156xx / 157xx / 158xx, seriali
% 15..20. Stessa catena surrogata di load_evaluation.m / FLORIDyn_step_*.m:
%   - TOWER: predittore NN - correttore NN residuo
%   - BLADE: solo predittore NN (secondo stadio disattivato, come in
%            load_evaluation.m)
%
% Cosa e' noto dal seriale <fam><nn>:
%   fam 156 -> SeaState.dat    Hs=3 m, Tp=8 s
%   fam 157 -> SeaState_v1.dat Hs=2 m, Tp=7 s
%   fam 158 -> SeaState_v2.dat Hs=4 m, Tp=10 s
%   nn 15/16 -> 8 m/s, 17/18 -> 7 m/s, 19/20 -> 9 m/s (vento nominale)
%   nn dispari -> layout A: y = [-210 -70 70] m, yaw T1/T2/T3 = [+15 +15 0]
%   nn pari    -> layout B: y = [70 -70 -210] m, yaw T1/T2/T3 = [-15 -15 0]
%   x = [300 1500 2700] m per entrambi i layout, D = 240 m
% WS e TI di ogni turbina vengono dall'Excel MLife (MeanWindspeed, TI).
% DD / YD / posizione nella wake con la stessa logica di LayoutDisp.m;
% YawUp = yaw della turbina a monte piu' vicina (0 se free stream).
%
% DEL FAST.Farm = ultima colonna di canale del blocco (Tower_6 / Blade_6)
% del foglio "Short-term DELs".
%
% OUTPUT: Excel con un foglio "Tower" e uno "Blade", una riga per
% seriale x turbina. "sum FF" e "sum surrogate" = somma dei DEL delle 3
% turbine dello stesso seriale (ripetuta sulle 3 righe).
%
% NOTA: i seriali pari di 156xx/157xx fino al 09/10/2026 avevano T1 a +15
% per errore (ElastoDyn sbagliato). Questo script assume gli yaw di
% progetto qui sopra: usarlo con l'Excel MLife delle simulazioni RILANCIATE.

clear;

%% ===== PERCORSI =====
delFile  = 'D:\zero_model\fastfarm_simulation\Mlife\CertTest\FLORIDyn_zeroMQ_fixed_yaw_DEL1_Short-term.xlsx';
delSheet = 'Short-term DELs';
caseDir  = 'D:\zero_model\FLORIDyn_Matlab-main_3_floating\FLORIDyn_Matlab-main_3\Simulations\2026_3T_DEL\3T_IEA15_floating_00';
outFile  = 'D:\zero_model\fastfarm_simulation\Mlife\CertTest\DEL_FF_vs_surrogate_fixed_yaw_DEL1.xlsx';

%% ===== PARAMETRI NOTI =====
seaStates = struct('fam', {156, 157, 158}, 'Hs', {3, 2, 4}, 'Tp', {8, 7, 10});
wsNominal = containers.Map({15, 16, 17, 18, 19, 20}, {8, 8, 7, 7, 9, 9});

layoutX = [300 1500 2700];
layoutY.A = [-210 -70 70];
layoutY.B = [70 -70 -210];
yawSet.A  = [15 15 0];
yawSet.B  = [-15 -15 0];
D = 240;

WindDir = 270;   % convenzione FLORIDyn (WD FAST.Farm = 0 -> 270)
WaveDir = 0.0;

%% ===== CARICA I MODELLI SURROGATI =====
here = pwd;
cd(caseDir);
cleanupObj = onCleanup(@() cd(here));

NN_predictor_blade = load('del_blade_uncertainty_grouped_wind_dir.mat');
NN_predictor_tower = load('del_tower_uncertainty_grouped_wind_dir.mat');
NN_corrector_tower = load('nn_residual_tower_wind_dir.mat');

wind_dir_training = 270 - WindDir;
sin_wd = sind(wind_dir_training);
cos_wd = cosd(wind_dir_training);

%% ===== LEGGI I DEL FAST.Farm (MLife) =====
raw = readcell(delFile, 'Sheet', delSheet);
FF.tower = readMLifeBlock(raw, 'Tower_Base', 'Tower_');
FF.blade = readMLifeBlock(raw, 'Blade_Root', 'Blade_');

%% ===== CALCOLO SURROGATO =====
loadTypes = {'tower', 'blade'};
sheetNames = {'Tower', 'Blade'};

for iL = 1:numel(loadTypes)
    lt  = loadTypes{iL};
    blk = FF.(lt);
    n   = numel(blk.serial);

    seriale = blk.serial;   turbina = blk.turbine;
    layout  = strings(n,1); seastate = strings(n,1);
    wsNom   = nan(n,1);     Hs_v = nan(n,1); Tp_v = nan(n,1);
    Yaw_v   = nan(n,1);     YawUp_v = nan(n,1);
    DD_v    = nan(n,1);     YD_v = nan(n,1);  Turb_v = nan(n,1);
    DEL_sur = nan(n,1);

    for k = 1:n
        fam = floor(seriale(k) / 100);
        nn  = mod(seriale(k), 100);
        t   = turbina(k);

        ss = seaStates([seaStates.fam] == fam);
        if isempty(ss) || ~isKey(wsNominal, nn)
            warning('Seriale %d non previsto, salto.', seriale(k));
            continue
        end
        if mod(nn, 2) == 1, lay = 'A'; else, lay = 'B'; end

        [DD, YD, row, iUp] = layoutDispUp(WindDir, layoutX, layoutY.(lay), D);
        yaw = yawSet.(lay);
        if iUp(t) > 0, YawUp = yaw(iUp(t)); else, YawUp = 0; end

        WS = blk.ws(k);
        TI = blk.ti(k);

        if strcmp(lt, 'tower')
            [DEL_NN, ~, ~] = predictTowerNN(NN_predictor_tower, yaw(t), ...
                ss.Hs, ss.Tp, WaveDir, WS, TI, wind_dir_training);
            res = residualTowerNN(NN_corrector_tower, DEL_NN, DD(t), YD(t), ...
                YawUp, YawUp, ss.Hs, ss.Tp, WaveDir, WS, TI, row(t), sin_wd, cos_wd);
            DEL_sur(k) = DEL_NN - res;
        else
            [DEL_sur(k), ~, ~] = predictBladeNN(NN_predictor_blade, yaw(t), ...
                ss.Hs, ss.Tp, WaveDir, WS, TI, wind_dir_training);
        end

        layout(k)   = lay;
        seastate(k) = sprintf('Hs=%g m, Tp=%g s', ss.Hs, ss.Tp);
        wsNom(k)    = wsNominal(nn);
        Hs_v(k) = ss.Hs;  Tp_v(k) = ss.Tp;
        Yaw_v(k) = yaw(t); YawUp_v(k) = YawUp;
        DD_v(k) = DD(t);  YD_v(k) = YD(t);  Turb_v(k) = row(t);
    end

    % Somme sulla farm (3 turbine dello stesso seriale)
    sumFF  = nan(n,1);
    sumSur = nan(n,1);
    for s = unique(seriale)'
        idx = seriale == s;
        sumFF(idx)  = sum(blk.del(idx));
        sumSur(idx) = sum(DEL_sur(idx));
    end
    mismatch = (DEL_sur - blk.del) ./ blk.del * 100;

    out = table(seriale, turbina, layout, seastate, wsNom, blk.ws, blk.ti, ...
        blk.del, DEL_sur, mismatch, sumFF, sumSur, ...
        Hs_v, Tp_v, Yaw_v, YawUp_v, DD_v, YD_v, Turb_v, ...
        'VariableNames', {'seriale', 'turbina', 'layout', 'seastate', ...
        'wind speed nominal', 'wind speed', 'TI', ...
        'FastFARM DEL', 'surrogate DEL', 'mismatch %', 'sum FF', 'sum surrogate', ...
        'Hs', 'Tp', 'Yaw', 'YawUp', 'DD', 'YD', 'Turb'});
    out = sortrows(out, {'seriale', 'turbina'});

    writetable(out, outFile, 'Sheet', sheetNames{iL});
    fprintf('%s: %d righe scritte nel foglio "%s"\n', lt, height(out), sheetNames{iL});
end

fprintf('\nRisultati salvati in:\n  %s\n', outFile);


%% ===================== FUNZIONI LOCALI =====================

function blk = readMLifeBlock(raw, blockTitle, channelPrefix)
% Legge un blocco del foglio "Short-term DELs" di MLife: per ogni riga
% <seriale>DEL<n>.T<t> prende l'ultima colonna di canale (es. Tower_6),
% MeanWindspeed e TI.
    isTxt = @(c) ischar(c) || isstring(c);
    col1  = raw(:,1);
    titleRow = find(cellfun(@(c) isTxt(c) && startsWith(strtrim(string(c)), blockTitle), col1), 1);
    if isempty(titleRow)
        error('Blocco "%s" non trovato nel foglio.', blockTitle);
    end
    hdrRow = titleRow + 1;
    hdr = raw(hdrRow,:);
    hdrStr = strings(size(hdr));
    for j = 1:numel(hdr)
        if isTxt(hdr{j}), hdrStr(j) = strtrim(string(hdr{j})); end
    end
    cFile = find(hdrStr == "Filename", 1);
    cChan = find(startsWith(hdrStr, channelPrefix), 1, 'last');
    cWS   = find(hdrStr == "MeanWindspeed", 1);
    cTI   = find(hdrStr == "TI", 1);

    serial = []; turbine = []; del = []; ws = []; ti = [];
    started = false;
    for r = hdrRow+1:size(raw,1)
        f = raw{r, cFile};
        tok = [];
        if isTxt(f), tok = regexp(char(f), '^(\d+)DEL(\d+)\.T(\d+)$', 'tokens', 'once'); end
        if isempty(tok)
            if started, break; else, continue; end   % salta L_Ult / L_MF
        end
        started = true;
        serial(end+1,1)  = str2double(tok{1});
        turbine(end+1,1) = str2double(tok{3});
        del(end+1,1)     = raw{r, cChan};
        ws(end+1,1)      = raw{r, cWS};
        ti(end+1,1)      = raw{r, cTI};
    end
    blk = struct('serial', serial, 'turbine', turbine, 'del', del, 'ws', ws, 'ti', ti);
    fprintf('Blocco %s: %d righe, colonna DEL = %s\n', blockTitle, numel(del), hdrStr(cChan));
end

function [DD, YD, row_idx, iUp] = layoutDispUp(wd, layoutX, layoutY, D)
% Copia di Simulations\2026_3T_DEL\3T_IEA15_floating_00\LayoutDisp.m (la
% versione usata dal controller, k = 0.1, raggio di scia in metri), che
% in piu' restituisce l'indice della turbina a monte piu' vicina
% (idx_next; 0 se free stream).
    k = 0.1;
    layoutX = layoutX(:);  layoutY = layoutY(:);
    nT = numel(layoutX);

    DD = 8 * ones(nT,1);
    YD = 1 * ones(nT,1);
    iUp = zeros(nT,1);
    row_idx = ones(nT,1);

    theta = wd - 270;
    mat_rot = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
    rot_lay = mat_rot * [layoutX.'; layoutY.'];
    x_rot = rot_lay(1,:).';
    y_rot = rot_lay(2,:).';
    R = D/2;

    for i = 1:nT
        DD_candidates = []; YD_candidates = []; idx_candidates = [];
        for j = 1:nT
            if j == i, continue; end
            DDtmp = x_rot(i) - x_rot(j);
            YDtmp = y_rot(i) - y_rot(j);
            if DDtmp > 0 && abs(YDtmp) <= R + k * DDtmp
                DD_candidates(end+1)  = DDtmp / D;
                YD_candidates(end+1)  = YDtmp / D;
                idx_candidates(end+1) = j;
            end
        end
        if ~isempty(DD_candidates)
            [~, imin] = min(DD_candidates);
            DD(i) = round(DD_candidates(imin));
            YD(i) = round(YD_candidates(imin), 1);
            iUp(i) = idx_candidates(imin);
            row_idx(i) = numel(DD_candidates) + 1;
        end
    end
end
