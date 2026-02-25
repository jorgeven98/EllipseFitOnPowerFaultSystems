%% TEST DEL EXPERIMENTO - VERSIÓN REDUCIDA
% Prueba con solo 2 archivos para verificar funcionamiento
%
% Fecha: 2026-02-10

clear; clc; close all;

% Agregar rutas a clasificadores, utilidades y experimentos
addpath(fullfile('..', 'classifiers'));
addpath(fullfile('..', 'utils'));
addpath(fullfile('..', 'experiments'));

%% CONFIGURACIÓN REDUCIDA
fprintf('========================================\n');
fprintf('TEST DEL EXPERIMENTO (versión reducida)\n');
fprintf('========================================\n\n');

% Parámetros temporales
Ts = 2.5e-4;
fs = 1/Ts;
f = 50;
fault_start = 3.0;
fault_end = 3.16;

% Solo 2 tamaños de ventana para prueba rápida
window_sizes = [80, 40];  % 1 ciclo y 1/2 ciclo
window_names = {'1 ciclo', '1/2 ciclo'};
stride_fraction = 0.25;

% Convertidores
converters = {'DG1', 'Pex'};
converter_types = {'Grid Forming', 'Grid Following'};

% Solo 2 métodos para test rápido (GA y Clarke)
methods = {
    @(V) ga_fault_classifier(V, fs, f), 'GA', 'Geometric Algebra';
    @(V) clarke_fault_classifier(V), 'Clarke', 'Clarke Transform';
};

% Solo 2 archivos para test
fault_files = {
    'DataGrid_AG.mat', 'A-G';
    'DataGridRL_AB.mat', 'A-B';
};

fprintf('Test con:\n');
fprintf('  - %d archivos\n', size(fault_files, 1));
fprintf('  - %d métodos\n', size(methods, 1));
fprintf('  - %d convertidores\n', length(converters));
fprintf('  - %d ventanas\n', length(window_sizes));
fprintf('  = %d configuraciones\n\n', size(fault_files,1) * size(methods,1) * length(converters) * length(window_sizes));


% Crear carpetas de resultados si no existen
if ~exist('results', 'dir')
    mkdir('results');
end
if ~exist('results/figures', 'dir')
    mkdir('results/figures');
end

%% EJECUTAR TEST
all_results = {};
result_counter = 1;

tic;

for file_idx = 1:size(fault_files, 1)
    fault_file = fault_files{file_idx, 1};
    expected_fault = fault_files{file_idx, 2};

    fprintf('[%d/%d] %s (%s)\n', file_idx, size(fault_files, 1), fault_file, expected_fault);

    % Cargar datos
    try
        data_path = fullfile('data', fault_file);
        data = load(data_path);
        fprintf('  ✓ Archivo cargado\n');
    catch ME
        fprintf('  ✗ Error: %s\n', ME.message);
        continue;
    end

    % Por cada método
    for method_idx = 1:size(methods, 1)
        classifier_func = methods{method_idx, 1};
        method_name = methods{method_idx, 2};
        method_desc = methods{method_idx, 3};

        % Por cada convertidor
        for conv_idx = 1:length(converters)
            converter_name = converters{conv_idx};

            fprintf('  %s-%s: ', method_name, converter_name);

            % Extraer voltajes
            try
                [V_abc_full, time_full] = extract_converter_data(data, converter_name);
                fprintf('datos OK, ');
            catch ME
                fprintf('ERROR - %s\n', ME.message);
                continue;
            end

            % Por cada ventana
            for win_idx = 1:length(window_sizes)
                window_size = window_sizes(win_idx);
                window_name = window_names{win_idx};
                stride = round(window_size * stride_fraction);

                % Ventana móvil
                classifications = {};
                time_centers = [];
                start_idx = 1;

                while start_idx + window_size - 1 <= length(V_abc_full)
                    end_idx = start_idx + window_size - 1;
                    V_window = V_abc_full(start_idx:end_idx, :);

                    % Clasificar con método actual
                    try
                        [fault_type, ~] = classifier_func(V_window);
                    catch ME
                        fault_type = 'Error';
                    end

                    center_idx = round((start_idx + end_idx) / 2);
                    time_center = time_full(center_idx);
                    classifications{end+1} = fault_type;
                    time_centers(end+1) = time_center;

                    start_idx = start_idx + stride;
                end

                fprintf('%s (%d), ', window_name, length(classifications));

                % Calcular métricas
                try
                    metrics = compute_detection_metrics(...
                        classifications, time_centers, expected_fault, ...
                        fault_start, fault_end);
                catch ME
                    fprintf('ERROR métricas - %s\n', ME.message);
                    continue;
                end

                % Guardar
                all_results{result_counter} = struct(...
                    'fault_type', expected_fault, ...
                    'method', method_name, ...
                    'method_description', method_desc, ...
                    'converter', converter_name, ...
                    'converter_type', converter_types{conv_idx}, ...
                    'window_size', window_size, ...
                    'window_name', window_name, ...
                    'classifications', {classifications}, ...
                    'time_vector', time_centers, ...
                    'metrics', metrics);

                result_counter = result_counter + 1;
            end

            fprintf('OK\n');
        end
    end

    fprintf('\n');
end

elapsed = toc;

%% RESUMEN
fprintf('========================================\n');
fprintf('TEST COMPLETADO\n');
fprintf('========================================\n\n');
fprintf('Tiempo: %.1f segundos\n', elapsed);
fprintf('Configuraciones: %d\n\n', length(all_results));

% Mostrar algunas métricas
if ~isempty(all_results)
    fprintf('Resultados de ejemplo:\n');
    for i = 1:min(6, length(all_results))
        r = all_results{i};
        fprintf('  %s | %s | %s | %s: Precisión=%.1f%%, Det=%.0f ms\n', ...
            r.fault_type, r.method, r.converter, r.window_name, ...
            r.metrics.fault_accuracy * 100, ...
            r.metrics.detection_time * 1000);
    end
end

%% PRUEBA DE EXPORTACIÓN
fprintf('\nProbando exportación...\n');
try
    export_results_to_csv(all_results, 'results/test_metrics.csv');
catch ME
    fprintf('✗ Error en export: %s\n', ME.message);
end

%% PRUEBA DE GRÁFICOS
fprintf('\nProbando gráficos...\n');
try
    figs = generate_comparison_plots(all_results);
    fprintf('✓ %d gráficos generados\n', length(figs));

    % Cerrar figuras
    for i = 1:length(figs)
        close(figs(i));
    end
catch ME
    fprintf('✗ Error en gráficos: %s\n', ME.message);
end

%% PRUEBA DE TIMELINE
fprintf('\nProbando timeline...\n');
try
    if ~isempty(all_results)
        fig = plot_detection_timeline(all_results{1}, fault_start, fault_end);
        fprintf('✓ Timeline generado\n');
        close(fig);
    end
catch ME
    fprintf('✗ Error en timeline: %s\n', ME.message);
end

fprintf('\n========================================\n');
fprintf('✓ TEST FINALIZADO\n');
fprintf('========================================\n');
