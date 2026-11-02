base_path = 'D:\zero_model\FLORIDyn_Matlab-main_3_floating\FLORIDyn_Matlab-main_3\';

fig_base = openfig(fullfile(base_path, 'zero_yaw.fig'),     'invisible');
fig_ctrl = openfig(fullfile(base_path, 'zero_non_yaw.fig'), 'invisible');

%% Estrai axes ordinati per posizione
ax_ctrl = findall(fig_ctrl, 'Type', 'axes');
ax_base = findall(fig_base, 'Type', 'axes');

pos_ctrl = cell2mat(get(ax_ctrl, 'Position'));
pos_base = cell2mat(get(ax_base, 'Position'));
[~, i_ctrl] = sortrows(pos_ctrl, [2 1]);
[~, i_base] = sortrows(pos_base, [2 1]);
ax_ctrl = ax_ctrl(i_ctrl);
ax_base = ax_base(i_base);

%% Copia le linee dalla baseline nel controlled
nAx = min(numel(ax_ctrl), numel(ax_base));
for k = 1:nAx
    % Copia linee baseline
    lines_base = findall(ax_base(k), 'Type', 'line');
    for l = 1:numel(lines_base)
        h = copyobj(lines_base(l), ax_ctrl(k));
        set(h, 'Color', [0.8500 0.3250 0.0980], ...
               'LineStyle', '--', ...
               'DisplayName', 'Baseline (yaw=0)');
    end

    % Rinomina linee originali
    lines_orig = findall(ax_ctrl(k), 'Type', 'line');
    for l = 1:numel(lines_orig)
        if strcmp(get(lines_orig(l), 'DisplayName'), '')
            set(lines_orig(l), 'DisplayName', 'Controlled');
        end
    end

    % Legend solo sul primo axes, piccola
    if k == 1
        lg = legend(ax_ctrl(k), 'show', 'Location', 'best');
        lg.FontSize = 8;
    else
        legend(ax_ctrl(k), 'off');
    end
end

%% Uniforma ylim per riga (WD, Yaw, DEL blade, DEL tower)
nT    = 3;
nRows = 4;

for row = 1:nRows
    idx_row = (row-1)*nT + (1:nT);

    % Raccogli tutti i dati Y di questa riga
    all_y_row = [];
    for k = idx_row
        if k > nAx; break; end
        all_lines = findall(ax_ctrl(k), 'Type', 'line');
        for l = 1:numel(all_lines)
            all_y_row = [all_y_row; get(all_lines(l), 'YData')'];
        end
    end

    % Calcola e applica limiti globali per la riga
    if ~isempty(all_y_row)
        ymin = min(all_y_row);
        ymax = max(all_y_row);
        if ymax - ymin < 1e-6
            ylims_row = [ymin-1, ymax+1];
        else
            margin = 0.05 * (ymax - ymin);
            ylims_row = [ymin - margin, ymax + margin];
        end
        for k = idx_row
            if k > nAx; break; end
            ylim(ax_ctrl(k), ylims_row);
        end
    end
end

fig_ctrl.Name    = 'Combined comparison: Controlled vs Baseline';
fig_ctrl.Visible = 'on';

%% Salva
saveas(fig_ctrl, fullfile(base_path, 'comparison_combined.fig'));
saveas(fig_ctrl, fullfile(base_path, 'comparison_combined.png'));
disp('Saved.')