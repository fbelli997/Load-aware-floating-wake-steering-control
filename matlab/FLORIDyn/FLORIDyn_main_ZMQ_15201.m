% FLORIDyn_main_ZMQ.m — gira su Windows in parallelo a Python
input_path  = 'D:\zero_model\floridyn_input_15201.json';
output_path = 'D:\zero_model\floridyn_output_15201.json';

% ===== WHILE 1: aspetta prime misure e inizializza =====
disp('[FLORIDyn] Waiting for first measurements from FAST.Farm...')
while true
    if ~exist(input_path, 'file'); pause(0.02); continue; end
    info = dir(input_path);
    if info.bytes == 0; pause(0.05); continue; end
    try
        raw  = fileread(input_path);
        data = jsondecode(raw);
    catch
        pause(0.02); continue;
    end
    delete(input_path);
    break;  % ← esce dal primo while
end

WS_init = data.WS;
WD_init = data.WD;
fprintf('[FLORIDyn] First measurements: WS=%.2f m/s  WD=%.2f deg\n', WS_init, WD_init)

% Inizializza FLORIDyn con misure reali
state = FLORIDyn_init_15201(WS_init, WD_init);

% Rispondi con gamma=0 alla prima chiamata
output.gamma_opt = zeros(1, state.T.nT);
fid = fopen(output_path, 'w');
fprintf(fid, '%s', jsonencode(output));
fclose(fid);

disp('[FLORIDyn] Ready, entering control loop...')

% ===== WHILE 2: loop di controllo normale =====
while true
    if ~exist(input_path, 'file'); pause(0.02); continue; end
    info = dir(input_path);
    if info.bytes == 0; pause(0.05); continue; end
    try
        raw  = fileread(input_path);
        data = jsondecode(raw);
    catch
        pause(0.02); continue;
    end
    delete(input_path);

    WS_meas     = data.WS;
    WD_meas     = data.WD;
    Thrust_meas = data.Thrust;
    Power_meas  = data.Power;
    TI_meas     = data.TI;
    fprintf('FASTfarm Wind Speed: %.2f m/s\n', WS_meas)
    fprintf('FASTfarm Wind Dir:   %.2f deg\n', WD_meas)

    t2 = tic;
    [gamma_opt, state] = FLORIDyn_step_15201(state, WS_meas, WD_meas, Thrust_meas, Power_meas, TI_meas);
    fprintf('[FLORIDyn_step] Tempo FLORIDyn: %.2f s\n', toc(t2))

    output.gamma_opt = gamma_opt;
    fid = fopen(output_path, 'w');
    fprintf(fid, '%s', jsonencode(output));
    fclose(fid);
    disp(['[FLORIDyn] gamma_opt: ' num2str(gamma_opt)])
end