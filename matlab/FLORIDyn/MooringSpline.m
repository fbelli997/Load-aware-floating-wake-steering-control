function mooring_force_f = MooringSpline()
    % MooringSpline
    % Legge 'mooring_lookuptable_hori_only.csv' con 2 colonne:
    %   col1 = distance_to_anchor [m]
    %   col2 = force [kN] (assunto) -> convertita in N (x1e3)
    % Ritorna un function handle: F = mooring_force_f(dist)

    filename = 'mooring_lookuptable_hori_only.csv';

    % Leggi CSV (robusto)
    data = readmatrix(filename);  % Nx2
    distance_array = data(:,1);
    force_array    = data(:,2)*1e3;  % kN -> N

    % (opzionale) ordina nel caso il file non sia monotono
    [distance_array, idx] = sort(distance_array);
    force_array = force_array(idx);

    % Interpolante lineare (come CasADi 'linear')
    % Restituisce scalare o vettore (se dist è vettore)
    mooring_force_f = @(dist) interp1(distance_array, force_array, dist, 'linear', 'extrap');
end
