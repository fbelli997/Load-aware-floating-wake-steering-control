%% Rewrite SOWFA files to be readable by Becker's importSOWFAFile (space-delimited, 4 cols, 3rd col = string)
clear; clc;

mio_path   = 'D:\FLORIDyn_Matlab-main_3\FLORIDyn_Matlab-main_3\Simulations\2026_Fede\3T_HKN_12_to_15_deg\Data';
files = { 'SOWFA_generatorPower.csv', 'SOWFA_nacelleYaw.csv' };

keepTurb = [0 1 3];   % 3 turbine

for f = 1:numel(files)
    inFile = fullfile(mio_path, files{f});
    if ~isfile(inFile)
        error('Non trovo il file: %s', inFile);
    end

    [folder, base, ext] = fileparts(inFile);
    outFile = fullfile(folder, [base, '_filtered', ext]);

    % Leggi numeri (ignora header e righe vuote)
    M = readmatrix(inFile);
    if isempty(M) || size(M,2) < 4
        error('Non riesco a leggere 4 colonne numeriche da: %s', inFile);
    end

    turb = M(:,1);
    time = M(:,2);
    val  = M(:,4);   % 4a colonna = valore (power o yaw)

    mask = ismember(turb, keepTurb);
    turb = turb(mask);
    time = time(mask);
    val  = val(mask);

    % Ordina come tipico SOWFA (time, poi turbina)
    [~, idx] = sortrows([time turb], [1 2]);
    turb = turb(idx); time = time(idx); val = val(idx);

    % Scrivi header (puoi anche lasciarlo come nel file originale, ma qui è pulito)
    if contains(lower(base), 'nacelleyaw')
        headerLine = '#Turbine Time(s) dt(s) nacelleYaw(deg)';
    elseif contains(lower(base), 'generatorpower')
        headerLine = '#Turbine Time(s) dt(s) rotorGenerator(W)';
    else
        headerLine = '#Turbine Time(s) dt(s) Value';
    end

    % Scrittura: SPACE-delimited, 3a colonna STRINGA ("dt")
    fid = fopen(outFile, 'w');
    if fid < 0, error('Impossibile scrivere: %s', outFile); end

    fprintf(fid, '%s\n', headerLine);

    % Formato: turb (int) time (float) Var3 (string) value (float)
    for i = 1:numel(turb)
        fprintf(fid, '%.0f %.10g %s %.10g\n', turb(i), time(i), 'dt', val(i));
    end

    fclose(fid);

    fprintf('[OK] %s -> %s | righe: %d\n', inFile, outFile, numel(turb));
end
