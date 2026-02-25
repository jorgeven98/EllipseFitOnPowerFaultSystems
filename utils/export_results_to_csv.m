function export_results_to_csv(all_results, filename)
% EXPORT_RESULTS_TO_CSV Exporta resultados a archivo CSV
%
% Input:
%   all_results: Cell array con todas las estructuras de resultados
%   filename: Nombre del archivo CSV de salida
%
% Autor: Experiment Script
% Fecha: 2026-02-10

if isempty(all_results)
    warning('No hay resultados para exportar');
    return;
end

%% Preparar datos para tabla
n_results = length(all_results);

% Pre-alocar arrays
fault_type = cell(n_results, 1);
method = cell(n_results, 1);
converter = cell(n_results, 1);
window_name = cell(n_results, 1);
window_size = zeros(n_results, 1);

overall_accuracy = zeros(n_results, 1);
detection_time = zeros(n_results, 1);
false_positive_rate = zeros(n_results, 1);
fault_accuracy = zeros(n_results, 1);
recovery_time = zeros(n_results, 1);

n_windows_total = zeros(n_results, 1);
n_windows_pre = zeros(n_results, 1);
n_windows_during = zeros(n_results, 1);
n_windows_post = zeros(n_results, 1);

%% Extraer datos de cada resultado
for i = 1:n_results
    result = all_results{i};

    fault_type{i} = result.fault_type;
    method{i} = result.method;
    converter{i} = result.converter;
    window_name{i} = result.window_name;
    window_size(i) = result.window_size;

    overall_accuracy(i) = result.metrics.overall_accuracy;
    detection_time(i) = result.metrics.detection_time;
    false_positive_rate(i) = result.metrics.false_positive_rate;
    fault_accuracy(i) = result.metrics.fault_accuracy;
    recovery_time(i) = result.metrics.recovery_time;

    n_windows_total(i) = result.metrics.n_windows_total;
    n_windows_pre(i) = result.metrics.n_windows_pre;
    n_windows_during(i) = result.metrics.n_windows_during;
    n_windows_post(i) = result.metrics.n_windows_post;
end

%% Crear tabla
T = table(fault_type, method, converter, window_name, window_size, ...
    overall_accuracy, detection_time, false_positive_rate, ...
    fault_accuracy, recovery_time, ...
    n_windows_total, n_windows_pre, n_windows_during, n_windows_post);

% Renombrar columnas para mejor legibilidad
T.Properties.VariableNames = {...
    'Fault_Type', 'Method', 'Converter', 'Window_Name', 'Window_Size_Samples', ...
    'Overall_Accuracy', 'Detection_Time_s', 'False_Positive_Rate', ...
    'Fault_Accuracy', 'Recovery_Time_s', ...
    'Total_Windows', 'Pre_Fault_Windows', 'During_Fault_Windows', 'Post_Fault_Windows'};

%% Exportar a CSV
try
    writetable(T, filename);
    fprintf('✓ Resultados exportados a: %s\n', filename);
    fprintf('  - %d filas exportadas\n', height(T));
catch ME
    warning('Error al exportar CSV: %s', ME.message);
end

%% Generar resumen estadístico adicional
[filepath, name, ~] = fileparts(filename);
summary_filename = fullfile(filepath, [name, '_summary.txt']);

try
    fid = fopen(summary_filename, 'w');

    fprintf(fid, '=== RESUMEN DEL EXPERIMENTO ===\n\n');
    fprintf(fid, 'Total de configuraciones: %d\n', n_results);
    fprintf(fid, 'Métodos analizados: %s\n', strjoin(unique(method), ', '));
    fprintf(fid, 'Convertidores analizados: %s\n', strjoin(unique(converter), ', '));
    fprintf(fid, 'Tipos de falta: %s\n', strjoin(unique(fault_type), ', '));
    fprintf(fid, 'Tamaños de ventana: %s\n\n', strjoin(unique(window_name), ', '));

    fprintf(fid, '--- MÉTRICAS GLOBALES ---\n');
    fprintf(fid, 'Precisión promedio: %.2f%%\n', mean(fault_accuracy, 'omitnan') * 100);
    fprintf(fid, 'Tiempo detección promedio: %.1f ms\n', mean(detection_time, 'omitnan') * 1000);
    fprintf(fid, 'Tasa falsos positivos promedio: %.2f%%\n\n', mean(false_positive_rate, 'omitnan') * 100);

    fprintf(fid, '--- POR MÉTODO ---\n');
    unique_methods = unique(method);
    for m = 1:length(unique_methods)
        met = unique_methods{m};
        idx = strcmp(method, met);

        fprintf(fid, '\n%s:\n', met);
        fprintf(fid, '  Precisión promedio: %.2f%%\n', mean(fault_accuracy(idx), 'omitnan') * 100);
        fprintf(fid, '  Tiempo detección promedio: %.1f ms\n', mean(detection_time(idx), 'omitnan') * 1000);
    end

    fprintf(fid, '\n--- POR CONVERTIDOR ---\n');
    unique_converters = unique(converter);
    for c = 1:length(unique_converters)
        conv = unique_converters{c};
        idx = strcmp(converter, conv);

        fprintf(fid, '\n%s:\n', conv);
        fprintf(fid, '  Precisión promedio: %.2f%%\n', mean(fault_accuracy(idx), 'omitnan') * 100);
        fprintf(fid, '  Tiempo detección promedio: %.1f ms\n', mean(detection_time(idx), 'omitnan') * 1000);
    end

    fprintf(fid, '\n--- POR TAMAÑO DE VENTANA ---\n');
    unique_windows = unique(window_name);
    for w = 1:length(unique_windows)
        win = unique_windows{w};
        idx = strcmp(window_name, win);

        fprintf(fid, '\n%s:\n', win);
        fprintf(fid, '  Precisión promedio: %.2f%%\n', mean(fault_accuracy(idx), 'omitnan') * 100);
        fprintf(fid, '  Tiempo detección promedio: %.1f ms\n', mean(detection_time(idx), 'omitnan') * 1000);
    end

    fclose(fid);
    fprintf('✓ Resumen exportado a: %s\n', summary_filename);
catch ME
    warning('Error al generar resumen: %s', ME.message);
    if fid > 0
        fclose(fid);
    end
end

end
