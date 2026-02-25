%% COMPARE_CLASSIFIER_ACCURACY
% Comparacion exhaustiva de la precision de clasificacion de faltas
% de los tres metodos (GA, Clarke, PE3D) con diferentes tamanos de ventana.
%
% Metrica principal: porcentaje de acierto DURANTE la ventana de falla.
%
% Salidas:
%   - Tablas LaTeX (.tex) en results/latex/
%   - Figuras (.png/.fig) en results/figures/
%   - Tablas formateadas en consola
%   - Variables guardadas en results/classifier_accuracy_results.mat
%
% Autor: Jorge - Doctorado
% Fecha: 2026-02

clear; clc; close all;

% Agregar rutas a clasificadores y utilidades
addpath(fullfile('..', 'classifiers'));
addpath(fullfile('..', 'utils'));

%% =====================================================================
%  PARAMETROS DEL SISTEMA
%  =====================================================================
Ts = 2.5e-4;                       % Periodo de muestreo (s)
fs = 1/Ts;                          % 4000 Hz
f  = 50;                            % Frecuencia del sistema (Hz)
samples_per_cycle = round(fs / f);  % 80 muestras/ciclo

%% =====================================================================
%  PARAMETROS DE VENTANA
%  =====================================================================
window_sizes  = [80, 60, 40, 20];   % muestras: 1, 3/4, 1/2, 1/4 ciclo
window_labels = {'1', '3/4', '1/2', '1/4'};
window_names  = {'1 ciclo', '3/4 ciclo', '1/2 ciclo', '1/4 ciclo'};
overlap_pct   = 75;
stride_fraction = 1 - overlap_pct / 100;  % 0.25

n_windows = length(window_sizes);

%% =====================================================================
%  MAPA DE ARCHIVOS DE FALLAS (14 archivos, igual que en el GUI)
%  =====================================================================
fault_file_map = {
    'DataGridRL_AB.mat',   'AB',   'Bifasica A-B (RL)';
    'DataGridRL_AC.mat',   'CA',   'Bifasica C-A (RL)';
    'DataGridRL_BC.mat',   'BC',   'Bifasica B-C (RL)';
    'DataGrid_AG.mat',     'AG',   'Monofasica A-G';
    'DataGrid_BG.mat',     'BG',   'Monofasica B-G';
    'DataGrid_CG.mat',     'CG',   'Monofasica C-G';
    'DataGrid_ABG.mat',    'ABG',  'Doble linea A-B-G';
    'DataGrid_ACG.mat',    'CAG',  'Doble linea C-A-G';
    'DataGrid_BCG.mat',    'BCG',  'Doble linea B-C-G';
    'DataGrid_AB.mat',     'AB',   'Bifasica A-B';
    'DataGrid_BC.mat',     'BC',   'Bifasica B-C';
    'DataGrid_AC.mat',     'CA',   'Bifasica C-A';
};
n_files = size(fault_file_map, 1);

%% =====================================================================
%  CONVERTIDORES Y METODOS
%  =====================================================================
converters       = {'DG1', 'Pex'};
converter_labels = {'DG1 (Grid Forming)', 'Pex (Grid Following)'};
n_converters     = length(converters);

methods = {
    @(V) ga_fault_classifier(V, fs, f), 'GA';
    @(V) clarke_fault_classifier(V),    'Clarke';
    @(V) pe3d_fault_classifier(V),      'PE3D';
};
method_names = methods(:, 2)';
n_methods    = size(methods, 1);

%% =====================================================================
%  INICIALIZAR ALMACENAMIENTO
%  =====================================================================
% accuracy(file, window, method) = porcentaje de acierto en region de falla
accuracy_DG1 = NaN(n_files, n_windows, n_methods);
accuracy_Pex = NaN(n_files, n_windows, n_methods);

n_win_DG1     = zeros(n_files, n_windows, n_methods);
n_correct_DG1 = zeros(n_files, n_windows, n_methods);
n_win_Pex     = zeros(n_files, n_windows, n_methods);
n_correct_Pex = zeros(n_files, n_windows, n_methods);

%% =====================================================================
%  BUCLE PRINCIPAL DE PROCESAMIENTO
%  =====================================================================
fprintf('================================================================\n');
fprintf('  COMPARACION DE PRECISION DE CLASIFICADORES DE FALTAS\n');
fprintf('================================================================\n\n');
fprintf('Configuracion:\n');
fprintf('  Archivos:       %d\n', n_files);
fprintf('  Convertidores:  %s\n', strjoin(converters, ', '));
fprintf('  Metodos:        %s\n', strjoin(method_names, ', '));
fprintf('  Ventanas:       %s\n', strjoin(window_names, ', '));
fprintf('  Overlap:        %d%%\n', overlap_pct);
total_configs = n_files * n_converters * n_windows * n_methods;
fprintf('  Total configs:  %d\n\n', total_configs);

tic;

for file_idx = 1:n_files
    filename       = fault_file_map{file_idx, 1};
    expected_fault = fault_file_map{file_idx, 2};
    fault_desc     = fault_file_map{file_idx, 3};

    fprintf('[%02d/%02d] %-25s (esperado: %-4s) ', ...
        file_idx, n_files, filename, expected_fault);

    % Cargar archivo
    data_path = fullfile('..', 'data', filename);
    try
        file_data = load(data_path);
    catch ME
        fprintf('ERROR: %s\n', ME.message);
        continue;
    end

    for conv_idx = 1:n_converters
        converter_name = converters{conv_idx};

        % Extraer datos del convertidor
        try
            [V_abc_full, time_full] = extract_converter_data(file_data, converter_name);
        catch
            continue;
        end

        % Determinar region de falla (misma logica que el GUI)
        t_duration = time_full(end) - time_full(1);
        if t_duration >= 4.5
            fault_start = 3.0;
            fault_end   = 3.16;
        elseif t_duration >= 2.0
            fault_start = 1.5;
            fault_end   = 1.66;
        else
            t_center    = (time_full(1) + time_full(end)) / 2;
            fault_start = t_center - 0.08;
            fault_end   = t_center + 0.08;
        end

        % Indices de la region de falla
        idx_fault_start = find(time_full >= fault_start, 1, 'first');
        idx_fault_end   = find(time_full <= fault_end, 1, 'last');

        if isempty(idx_fault_start) || isempty(idx_fault_end)
            continue;
        end

        for win_idx = 1:n_windows
            window_size = window_sizes(win_idx);
            stride      = max(1, round(window_size * stride_fraction));

            for met_idx = 1:n_methods
                classifier_func = methods{met_idx, 1};

                correct_count = 0;
                total_count   = 0;
                start_sample  = idx_fault_start;

                % Ventana deslizante SOLO dentro de la region de falla
                while start_sample + window_size - 1 <= idx_fault_end
                    end_sample = start_sample + window_size - 1;
                    V_window   = V_abc_full(start_sample:end_sample, :);

                    try
                        [fault_type_raw, ~] = classifier_func(V_window);
                        fault_type_norm = normalize_fault_name(fault_type_raw);
                    catch
                        fault_type_norm = 'Unclassified';
                    end

                    if strcmpi(fault_type_norm, expected_fault)
                        correct_count = correct_count + 1;
                    end
                    total_count = total_count + 1;

                    start_sample = start_sample + stride;
                end

                % Calcular accuracy
                if total_count > 0
                    acc = 100 * correct_count / total_count;
                else
                    acc = NaN;
                end

                % Almacenar
                if conv_idx == 1
                    accuracy_DG1(file_idx, win_idx, met_idx)  = acc;
                    n_win_DG1(file_idx, win_idx, met_idx)     = total_count;
                    n_correct_DG1(file_idx, win_idx, met_idx)  = correct_count;
                else
                    accuracy_Pex(file_idx, win_idx, met_idx)  = acc;
                    n_win_Pex(file_idx, win_idx, met_idx)     = total_count;
                    n_correct_Pex(file_idx, win_idx, met_idx)  = correct_count;
                end
            end
        end
    end

    fprintf('OK\n');
end

elapsed_time = toc;
fprintf('\nProcesamiento completado en %.1f segundos.\n\n', elapsed_time);

%% =====================================================================
%  AGREGAR POR TIPO DE FALLA UNICO
%  =====================================================================
all_fault_types = fault_file_map(:, 2);
[~, first_idx]  = unique(all_fault_types, 'first');
unique_faults   = all_fault_types(sort(first_idx));
n_unique        = length(unique_faults);

% agg(unique_fault, window, method) = media de accuracy
agg_DG1 = NaN(n_unique, n_windows, n_methods);
agg_Pex = NaN(n_unique, n_windows, n_methods);

for uf = 1:n_unique
    mask = strcmp(all_fault_types, unique_faults{uf});
    for w = 1:n_windows
        for m = 1:n_methods
            agg_DG1(uf, w, m) = mean(accuracy_DG1(mask, w, m), 'omitnan');
            agg_Pex(uf, w, m) = mean(accuracy_Pex(mask, w, m), 'omitnan');
        end
    end
end

%% Detectar y reportar NaN residuales en la agregacion
nan_dg1 = isnan(agg_DG1);
nan_pex = isnan(agg_Pex);
if any(nan_dg1(:)) || any(nan_pex(:))
    fprintf('\n[AVISO] Se encontraron NaN en los resultados agregados:\n');
    for uf = 1:n_unique
        for w = 1:n_windows
            for m = 1:n_methods
                if nan_dg1(uf, w, m)
                    fprintf('  DG1 | Falta=%-5s | Ventana=%-8s | Metodo=%-8s -> NaN\n', ...
                        unique_faults{uf}, window_labels{w}, method_names{m});
                end
                if nan_pex(uf, w, m)
                    fprintf('  Pex | Falta=%-5s | Ventana=%-8s | Metodo=%-8s -> NaN\n', ...
                        unique_faults{uf}, window_labels{w}, method_names{m});
                end
            end
        end
    end
else
    fprintf('  Sin NaN en resultados agregados.\n');
end

%% =====================================================================
%  TABLAS EN CONSOLA
%  =====================================================================
for conv_idx = 1:n_converters
    if conv_idx == 1
        agg_data = agg_DG1;
    else
        agg_data = agg_Pex;
    end

    for w = 1:n_windows
        fprintf('\n=== %s | Ventana: %s ===\n', converter_labels{conv_idx}, window_names{w});
        fprintf('%-8s | %8s | %8s | %8s\n', 'Falla', 'GA', 'Clarke', 'PE3D');
        fprintf('%s\n', repmat('-', 1, 44));

        for uf = 1:n_unique
            fprintf('%-8s |', unique_faults{uf});
            for m = 1:n_methods
                val = agg_data(uf, w, m);
                if isnan(val)
                    fprintf(' %7s |', '--');
                else
                    fprintf(' %6.1f%% |', val);
                end
            end
            fprintf('\n');
        end

        fprintf('%s\n', repmat('-', 1, 44));
        fprintf('%-8s |', 'MEDIA');
        for m = 1:n_methods
            fprintf(' %6.1f%% |', mean(agg_data(:, w, m), 'omitnan'));
        end
        fprintf('\n');
    end
end

%% =====================================================================
%  GENERAR TABLAS LATEX
%  =====================================================================
latex_dir = fullfile('..', 'results', 'latex');
if ~exist(latex_dir, 'dir')
    mkdir(latex_dir);
end

% --- Tablas por convertidor (una tabla por ventana) ---
for conv_idx = 1:n_converters
    if conv_idx == 1
        agg_data = agg_DG1;
        conv_name = 'DG1';
    else
        agg_data = agg_Pex;
        conv_name = 'Pex';
    end

    tex_file = fullfile(latex_dir, sprintf('table_accuracy_%s.tex', conv_name));
    fid = fopen(tex_file, 'w');

    for w = 1:n_windows
        fprintf(fid, '%% Tabla: %s - Ventana %s ciclo\n', conv_name, window_labels{w});
        fprintf(fid, '\\begin{table}[htbp]\n');
        fprintf(fid, '\\centering\n');
        fprintf(fid, '\\caption{Precision de clasificacion (\\%%) para convertidor %s, ventana = %s ciclo.}\n', ...
            conv_name, window_labels{w});
        fprintf(fid, '\\label{tab:acc_%s_%s}\n', lower(conv_name), ...
            strrep(strrep(window_labels{w}, '/', ''), ' ', ''));
        fprintf(fid, '\\begin{tabular}{lSSS}\n');
        fprintf(fid, '\\toprule\n');
        fprintf(fid, '{Tipo de Falla} & {GA} & {Clarke} & {PE3D} \\\\\n');
        fprintf(fid, '\\midrule\n');

        for uf = 1:n_unique
            fprintf(fid, '%s', unique_faults{uf});
            for m = 1:n_methods
                val = agg_data(uf, w, m);
                if isnan(val)
                    fprintf(fid, ' & {--}');
                else
                    fprintf(fid, ' & %.1f', val);
                end
            end
            fprintf(fid, ' \\\\\n');
        end

        fprintf(fid, '\\midrule\n');
        fprintf(fid, '\\textbf{Media}');
        for m = 1:n_methods
            fprintf(fid, ' & \\textbf{%.1f}', mean(agg_data(:, w, m), 'omitnan'));
        end
        fprintf(fid, ' \\\\\n');
        fprintf(fid, '\\bottomrule\n');
        fprintf(fid, '\\end{tabular}\n');
        fprintf(fid, '\\end{table}\n\n');
    end

    fclose(fid);
    fprintf('Tabla LaTeX guardada: %s\n', tex_file);
end

% --- Tabla resumen (metodo x ventana, promedio global) ---
tex_summary = fullfile(latex_dir, 'table_accuracy_summary.tex');
fid = fopen(tex_summary, 'w');

fprintf(fid, '\\begin{table}[htbp]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Precision media de clasificacion (\\%%) por metodo y tamano de ventana.}\n');
fprintf(fid, '\\label{tab:acc_summary}\n');
fprintf(fid, '\\begin{tabular}{ll');
for w = 1:n_windows
    fprintf(fid, 'S');
end
fprintf(fid, '}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, '{Convertidor} & {Metodo}');
for w = 1:n_windows
    fprintf(fid, ' & {%s ciclo}', window_labels{w});
end
fprintf(fid, ' \\\\\n');
fprintf(fid, '\\midrule\n');

for conv_idx = 1:n_converters
    if conv_idx == 1
        agg_data = agg_DG1;
    else
        agg_data = agg_Pex;
    end

    for m = 1:n_methods
        if m == 1
            fprintf(fid, '\\multirow{%d}{*}{%s}', n_methods, converters{conv_idx});
        end
        fprintf(fid, ' & %s', method_names{m});
        for w = 1:n_windows
            fprintf(fid, ' & %.1f', mean(agg_data(:, w, m), 'omitnan'));
        end
        fprintf(fid, ' \\\\\n');
    end

    if conv_idx < n_converters
        fprintf(fid, '\\midrule\n');
    end
end

fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
fprintf('Tabla resumen LaTeX guardada: %s\n', tex_summary);

%% =====================================================================
%  FIGURAS
%  =====================================================================
fig_dir = fullfile('..', 'results', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% --- Figura 1: Barras agrupadas (accuracy media vs ventana) ---
fig1 = figure('Position', [100, 100, 900, 500], 'Color', 'w');

% Promediar ambos convertidores
combined_avg = NaN(n_unique, n_windows, n_methods);
for uf = 1:n_unique
    for w = 1:n_windows
        for m = 1:n_methods
            vals = [agg_DG1(uf, w, m), agg_Pex(uf, w, m)];
            combined_avg(uf, w, m) = mean(vals, 'omitnan');
        end
    end
end

bar_data = zeros(n_windows, n_methods);
for w = 1:n_windows
    for m = 1:n_methods
        bar_data(w, m) = mean(combined_avg(:, w, m), 'omitnan');
    end
end

b = bar(bar_data, 'grouped');
colors = [0.2 0.4 0.8; 0.8 0.2 0.2; 0.2 0.7 0.3];
for m = 1:n_methods
    b(m).FaceColor = colors(m, :);
end

set(gca, 'XTickLabel', window_names, 'FontSize', 11);
xlabel('Tamano de Ventana', 'FontSize', 12);
ylabel('Precision Media (%)', 'FontSize', 12);
title('Precision de Clasificacion vs. Tamano de Ventana', 'FontSize', 13);
legend(method_names, 'Location', 'southwest', 'FontSize', 11);
grid on;
ylim([0 105]);

% Etiquetas de valor sobre las barras
for m = 1:n_methods
    for w = 1:n_windows
        text(b(m).XEndPoints(w), b(m).YEndPoints(w) + 1.5, ...
            sprintf('%.1f', bar_data(w, m)), ...
            'HorizontalAlignment', 'center', 'FontSize', 8);
    end
end

saveas(fig1, fullfile(fig_dir, 'accuracy_vs_window_bar.png'));
savefig(fig1, fullfile(fig_dir, 'accuracy_vs_window_bar.fig'));
fprintf('Figura 1 guardada: accuracy_vs_window_bar.png\n');

% --- Figura 2: Heatmaps por convertidor y metodo ---
fig2 = figure('Position', [50, 50, 1400, 600], 'Color', 'w');

for conv_idx = 1:n_converters
    if conv_idx == 1
        agg_data = agg_DG1;
    else
        agg_data = agg_Pex;
    end

    for m = 1:n_methods
        subplot(n_converters, n_methods, (conv_idx - 1) * n_methods + m);

        hm_data = agg_data(:, :, m);  % [n_unique x n_windows]

        imagesc(hm_data);
        cb = colorbar;
        cb.Label.String = 'Precision (%)';
        caxis([0, 100]);
        colormap(gca, 'parula');

        set(gca, 'XTick', 1:n_windows, 'XTickLabel', window_labels);
        set(gca, 'YTick', 1:n_unique, 'YTickLabel', unique_faults);
        xlabel('Ventana (ciclos)');

        if m == 1
            ylabel(converters{conv_idx}, 'FontWeight', 'bold', 'FontSize', 12);
        end
        title(method_names{m}, 'FontWeight', 'bold');

        % Etiquetas numericas en cada celda
        for uf = 1:n_unique
            for ww = 1:n_windows
                val = hm_data(uf, ww);
                if ~isnan(val)
                    if val > 50
                        txt_color = 'k';
                    else
                        txt_color = 'w';
                    end
                    text(ww, uf, sprintf('%.0f', val), ...
                        'HorizontalAlignment', 'center', 'FontSize', 7, ...
                        'Color', txt_color, 'FontWeight', 'bold');
                end
            end
        end
    end
end

sgtitle('Precision de Clasificacion por Tipo de Falla y Metodo', 'FontSize', 14);
saveas(fig2, fullfile(fig_dir, 'accuracy_heatmap.png'));
savefig(fig2, fullfile(fig_dir, 'accuracy_heatmap.fig'));
fprintf('Figura 2 guardada: accuracy_heatmap.png\n');

%% =====================================================================
%  GUARDAR RESULTADOS
%  =====================================================================
save(fullfile('..', 'results', 'classifier_accuracy_results.mat'), ...
    'accuracy_DG1', 'accuracy_Pex', ...
    'agg_DG1', 'agg_Pex', ...
    'fault_file_map', 'unique_faults', ...
    'window_sizes', 'window_labels', 'window_names', ...
    'method_names', 'converters', ...
    'n_win_DG1', 'n_correct_DG1', ...
    'n_win_Pex', 'n_correct_Pex');

fprintf('\nResultados guardados en: results/classifier_accuracy_results.mat\n');
fprintf('Tablas LaTeX en: results/latex/\n');
fprintf('Figuras en: results/figures/\n');

fprintf('\n================================================================\n');
fprintf('  COMPARACION COMPLETADA\n');
fprintf('================================================================\n');

%% =====================================================================
%  FUNCIONES LOCALES
%  =====================================================================

function canonical = normalize_fault_name(raw_name)
% Convierte cualquier nombre de falta al formato canonico.
% Cubre las tres convenciones: GA, Clarke y PE3D.
    mapping = {
        % PE3D / canonicos
        'Normal',       'Normal';
        'ABCG',         'ABCG';
        'AG',           'AG';
        'BG',           'BG';
        'CG',           'CG';
        'AB',           'AB';
        'BC',           'BC';
        'CA',           'CA';
        'ABG',          'ABG';
        'BCG',          'BCG';
        'CAG',          'CAG';
        % GA y Clarke (mismo convenio de nombres con guiones)
        'A-B-C',        'Normal';
        'A-B-C-G',      'ABCG';
        'A-G',          'AG';
        'B-G',          'BG';
        'C-G',          'CG';
        'A-B',          'AB';
        'B-C',          'BC';
        'C-A',          'CA';
        'A-B-G',        'ABG';
        'B-C-G',        'BCG';
        'C-A-G',        'CAG';
        'L-L-G',        'LLG_unknown';
        'L-L',          'LL_unknown';
        % Fallback
        'Unclassified', 'Unclassified';
    };

    canonical = 'Unclassified';
    for k = 1:size(mapping, 1)
        if strcmpi(strtrim(raw_name), mapping{k, 1})
            canonical = mapping{k, 2};
            return;
        end
    end
end
