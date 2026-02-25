%% EXPERIMENTO DE DETECCIÓN DE FALTAS CON VENTANA MÓVIL
% Script para analizar el rendimiento del clasificador GA_ellipse_fault_classifier
% usando ventanas móviles de diferentes tamaños sobre datos reales de simulación
%
% Autor: Jorge - Experimento Doctorado
% Fecha: 2026-02-10

clear; clc; close all;

% Agregar rutas a clasificadores y utilidades
addpath(fullfile('..', 'classifiers'));
addpath(fullfile('..', 'utils'));

%% CONFIGURACIÓN DEL EXPERIMENTO
fprintf('========================================\n');
fprintf('EXPERIMENTO DE DETECCIÓN CON VENTANA MÓVIL\n');
fprintf('========================================\n\n');

% Parámetros temporales (según datos reales)
Ts = 2.5e-4;              % Periodo de muestreo (s)
fs = 1/Ts;                % 4000 Hz
f = 50;                   % Frecuencia del sistema (Hz)
fault_start = 3.0;        % Inicio de falta (s)
fault_end = 3.16;         % Fin de falta (s)

% Tamaños de ventana
% f=50Hz → T=20ms → 1 ciclo = 80 muestras @ 4kHz
window_sizes = [80, 60, 40, 20];  % muestras
window_names = {'1 ciclo', '3/4 ciclo', '1/2 ciclo', '1/4 ciclo'};

% Paso de ventana (solapamiento 75%)
stride_fraction = 0.25;

% Convertidores a analizar
converters = {'DG1', 'Pex'};
converter_types = {'Grid Forming', 'Grid Following'};

% Métodos de clasificación a comparar
methods = {
    @(V) ga_fault_classifier(V, fs, f), 'GA', 'Geometric Algebra';
    @(V) clarke_fault_classifier(V), 'Clarke', 'Clarke Transform';
    @(V) pe3d_fault_classifier(V), 'PE3D', 'Polarization Ellipse 3D';
};

fprintf('Parámetros del experimento:\n');
fprintf('  Frecuencia muestreo: %.0f Hz\n', fs);
fprintf('  Frecuencia sistema: %.0f Hz\n', f);
fprintf('  Inicio falta: %.2f s\n', fault_start);
fprintf('  Fin falta: %.2f s\n', fault_end);
fprintf('  Tamaños ventana: %s\n', strjoin(window_names, ', '));
fprintf('  Convertidores: %s\n', strjoin(converters, ', '));
fprintf('  Métodos: %s\n\n', strjoin(methods(:,2)', ', '));

%% DEFINIR ARCHIVOS A PROCESAR

% Archivos principales (con impedancia RL)
fault_files_RL = {
    'DataGridRL_AB.mat', 'A-B';
    'DataGridRL_AC.mat', 'C-A';
    'DataGridRL_BC.mat', 'B-C';
    'DataGridRL_ABC.mat', 'A-B-C';
    'DataGrid_AG.mat', 'A-G';
    'DataGrid_BG.mat', 'B-G';
    'DataGrid_CG.mat', 'C-G';
    'DataGrid_ABG.mat', 'A-B-G';
    'DataGrid_ACG.mat', 'C-A-G';
    'DataGrid_BCG.mat', 'B-C-G';
    'DataGrid_ABCG.mat', 'A-B-C-G';
};

% Archivos especiales (sin impedancia RL - forman líneas)
fault_files_noRL = {
    'DataGrid_AB.mat', 'A-B';
    'DataGrid_BC.mat', 'B-C';
    'DataGrid_AC.mat', 'C-A';
};

% Usar solo archivos RL por defecto (cambiar a fault_files_noRL para casos especiales)
fault_files = fault_files_RL;

fprintf('Archivos a procesar: %d\n', size(fault_files, 1));
fprintf('Convertidores: %d\n', length(converters));
fprintf('Ventanas: %d\n', length(window_sizes));
fprintf('Métodos: %d\n', size(methods, 1));
fprintf('Total configuraciones: %d\n\n', size(fault_files, 1) * length(converters) * length(window_sizes) * size(methods, 1));

%% INICIALIZAR ALMACENAMIENTO DE RESULTADOS
all_results = {};
result_counter = 1;

% Añadir carpeta utils al path
addpath('utils');

% Crear carpetas de resultados si no existen
if ~exist('results', 'dir')
    mkdir('results');
end
if ~exist('results/figures', 'dir')
    mkdir('results/figures');
end

%% LOOP PRINCIPAL DEL EXPERIMENTO
fprintf('Iniciando experimento...\n');
fprintf('========================================\n\n');

tic;  % Iniciar cronómetro

for file_idx = 1:size(fault_files, 1)
    fault_file = fault_files{file_idx, 1};
    expected_fault = fault_files{file_idx, 2};

    fprintf('[%d/%d] Procesando: %s (Falta: %s)\n', file_idx, size(fault_files, 1), ...
        fault_file, expected_fault);

    % Cargar datos
    try
        data_path = fullfile('data', fault_file);
        data = load(data_path);
    catch ME
        fprintf('  ✗ Error cargando archivo: %s\n', ME.message);
        continue;
    end

    % Por cada método de clasificación
    for method_idx = 1:size(methods, 1)
        classifier_func = methods{method_idx, 1};
        method_name = methods{method_idx, 2};
        method_desc = methods{method_idx, 3};

        fprintf('  Método: %s (%s)\n', method_name, method_desc);

        % Por cada convertidor
        for conv_idx = 1:length(converters)
            converter_name = converters{conv_idx};

            fprintf('    Convertidor: %s (%s)\n', converter_name, converter_types{conv_idx});

            % Extraer voltajes
            try
                [V_abc_full, time_full] = extract_converter_data(data, converter_name);
            catch ME
                fprintf('      ✗ Error extrayendo datos: %s\n', ME.message);
                continue;
            end

            % Por cada tamaño de ventana
            for win_idx = 1:length(window_sizes)
                window_size = window_sizes(win_idx);
                window_name = window_names{win_idx};
                stride = round(window_size * stride_fraction);

                fprintf('      Ventana: %s (%d muestras, stride=%d)\n', ...
                    window_name, window_size, stride);

                % Inicializar vectores de resultados
                classifications = {};
                time_centers = [];

                % VENTANA MÓVIL
                start_idx = 1;
                n_windows = 0;

                while start_idx + window_size - 1 <= length(V_abc_full)
                    end_idx = start_idx + window_size - 1;

                    % Extraer ventana
                    V_window = V_abc_full(start_idx:end_idx, :);

                    % Clasificar usando el método actual
                    try
                        [fault_type, params] = classifier_func(V_window);
                    catch ME
                        fault_type = 'Error';
                    end

                    % Tiempo central de ventana
                    center_idx = round((start_idx + end_idx) / 2);
                    time_center = time_full(center_idx);

                    % Guardar resultado
                    classifications{end+1} = fault_type;
                    time_centers(end+1) = time_center;

                    n_windows = n_windows + 1;

                    % Avanzar ventana
                    start_idx = start_idx + stride;
                end

                fprintf('        Procesadas %d ventanas\n', n_windows);

                % Calcular métricas
                try
                    metrics = compute_detection_metrics(...
                        classifications, time_centers, expected_fault, ...
                        fault_start, fault_end);
                catch ME
                    fprintf('        ✗ Error calculando métricas: %s\n', ME.message);
                    continue;
                end

                % Guardar resultados
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

                fprintf('        ✓ Precisión: %.1f%%, Detección: %.0f ms\n', ...
                    metrics.fault_accuracy * 100, metrics.detection_time * 1000);
            end
        end
    end

    fprintf('\n');
end

elapsed_time = toc;

%% RESUMEN DEL EXPERIMENTO
fprintf('========================================\n');
fprintf('EXPERIMENTO COMPLETADO\n');
fprintf('========================================\n\n');
fprintf('Tiempo total: %.1f minutos\n', elapsed_time / 60);
fprintf('Configuraciones procesadas: %d\n\n', length(all_results));

%% GUARDAR RESULTADOS
fprintf('Guardando resultados...\n');

% Guardar estructura completa
save('results/all_results.mat', 'all_results', 'fault_files', 'converters', ...
    'window_sizes', 'window_names', 'fault_start', 'fault_end', 'fs', 'f');
fprintf('✓ Resultados guardados en: results/all_results.mat\n');

% Exportar a CSV
export_results_to_csv(all_results, 'results/detection_metrics.csv');

% Generar tabla comparativa de métodos
fprintf('\nGenerando tabla comparativa...\n');
generate_comparison_table(all_results, 'results/comparison_table');

%% GENERAR VISUALIZACIONES
fprintf('\nGenerando visualizaciones...\n');

try
    figs = generate_comparison_plots(all_results);

    % Guardar figuras
    for i = 1:length(figs)
        filename = sprintf('results/figures/comparison_plot_%d.png', i);
        saveas(figs(i), filename);
        fprintf('  ✓ Figura %d guardada\n', i);
    end
catch ME
    fprintf('  ✗ Error generando gráficos: %s\n', ME.message);
end

%% GENERAR ALGUNOS TIMELINES EJEMPLO
fprintf('\nGenerando timelines de ejemplo...\n');

% Seleccionar algunos casos interesantes
example_indices = [1, length(all_results)/2, length(all_results)];

for i = 1:min(3, length(all_results))
    idx = round(example_indices(i));
    if idx > 0 && idx <= length(all_results)
        try
            fig = plot_detection_timeline(all_results{idx}, fault_start, fault_end);
            filename = sprintf('results/figures/timeline_%s_%s_%s.png', ...
                all_results{idx}.fault_type, ...
                all_results{idx}.converter, ...
                all_results{idx}.window_name);
            % Reemplazar caracteres inválidos en nombre de archivo
            filename = strrep(filename, '-', '_');
            filename = strrep(filename, ' ', '_');
            saveas(fig, filename);
            close(fig);
            fprintf('  ✓ Timeline guardado: %s\n', filename);
        catch ME
            fprintf('  ✗ Error en timeline: %s\n', ME.message);
        end
    end
end

%% FINALIZAR
fprintf('\n========================================\n');
fprintf('EXPERIMENTO FINALIZADO EXITOSAMENTE\n');
fprintf('========================================\n\n');
fprintf('Resultados disponibles en:\n');
fprintf('  - results/all_results.mat\n');
fprintf('  - results/detection_metrics.csv\n');
fprintf('  - results/detection_metrics_summary.txt\n');
fprintf('  - results/figures/\n\n');
