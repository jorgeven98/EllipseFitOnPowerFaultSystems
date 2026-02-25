function generate_comparison_table(all_results, output_filename)
% GENERATE_COMPARISON_TABLE Genera tabla comparativa de métodos
%
% Crea una tabla con:
% - Filas: Tipos de falta
% - Columnas: Método × Tamaño de Ventana × Convertidor
% - Valores: Porcentaje de acierto durante falta
%
% Input:
%   all_results: Cell array con resultados del experimento
%   output_filename: Nombre del archivo de salida (opcional)
%
% Output:
%   Archivos generados:
%   - Tabla CSV con precisión por falta/método/ventana/convertidor
%   - Tabla TXT formateada para lectura humana
%
% Autor: Experiment Script
% Fecha: 2026-02-10

if nargin < 2
    output_filename = 'results/comparison_table';
end

if isempty(all_results)
    warning('No hay resultados para generar tabla');
    return;
end

%% Extraer información única
fault_types = unique(cellfun(@(x) x.fault_type, all_results, 'UniformOutput', false));
methods = unique(cellfun(@(x) x.method, all_results, 'UniformOutput', false));
converters = unique(cellfun(@(x) x.converter, all_results, 'UniformOutput', false));
windows = unique(cellfun(@(x) x.window_name, all_results, 'UniformOutput', false));

%% Crear tabla de precisión durante falta
n_faults = length(fault_types);
n_combos = length(methods) * length(windows) * length(converters);

% Matriz de datos
data_matrix = zeros(n_faults, n_combos);

% Generar nombres de columnas
column_names = {};
for m = 1:length(methods)
    for w = 1:length(windows)
        for c = 1:length(converters)
            col_name = sprintf('%s_%s_%s', methods{m}, strrep(windows{w}, ' ', ''), converters{c});
            column_names{end+1} = col_name;
        end
    end
end

% Llenar matriz con precisiones
for f = 1:n_faults
    fault = fault_types{f};
    col_idx = 1;

    for m = 1:length(methods)
        for w = 1:length(windows)
            for c = 1:length(converters)
                % Buscar resultado
                idx = find_result(all_results, fault, methods{m}, windows{w}, converters{c});

                if ~isempty(idx)
                    data_matrix(f, col_idx) = all_results{idx}.metrics.fault_accuracy * 100;
                else
                    data_matrix(f, col_idx) = NaN;
                end

                col_idx = col_idx + 1;
            end
        end
    end
end

%% Crear tabla MATLAB
T = array2table(data_matrix, 'VariableNames', column_names, 'RowNames', fault_types);

%% Exportar a CSV
csv_file = [output_filename, '.csv'];
writetable(T, csv_file, 'WriteRowNames', true);
fprintf('✓ Tabla comparativa guardada en: %s\n', csv_file);

%% Generar tabla formateada para texto
txt_file = [output_filename, '.txt'];
fid = fopen(txt_file, 'w');

fprintf(fid, '================================================================================\n');
fprintf(fid, 'TABLA COMPARATIVA: PRECISIÓN POR MÉTODO, VENTANA Y CONVERTIDOR\n');
fprintf(fid, '================================================================================\n\n');

fprintf(fid, 'Valores: Precisión durante falta (%%)\n');
fprintf(fid, 'Métodos: %s\n', strjoin(methods, ', '));
fprintf(fid, 'Convertidores: %s\n', strjoin(converters, ', '));
fprintf(fid, 'Ventanas: %s\n\n', strjoin(windows, ', '));

% Tabla por método
for m = 1:length(methods)
    fprintf(fid, '--- MÉTODO: %s ---\n\n', methods{m});

    % Encabezado
    fprintf(fid, '%-10s', 'Falta');
    for w = 1:length(windows)
        for c = 1:length(converters)
            fprintf(fid, ' | %-12s', sprintf('%s-%s', strrep(windows{w}, ' ', ''), converters{c}));
        end
    end
    fprintf(fid, '\n');
    fprintf(fid, repmat('-', 1, 10 + length(windows)*length(converters)*15));
    fprintf(fid, '\n');

    % Datos
    for f = 1:n_faults
        fprintf(fid, '%-10s', fault_types{f});

        for w = 1:length(windows)
            for c = 1:length(converters)
                idx = find_result(all_results, fault_types{f}, methods{m}, windows{w}, converters{c});

                if ~isempty(idx)
                    val = all_results{idx}.metrics.fault_accuracy * 100;
                    fprintf(fid, ' | %10.1f%%', val);
                else
                    fprintf(fid, ' | %12s', 'N/A');
                end
            end
        end
        fprintf(fid, '\n');
    end

    fprintf(fid, '\n');
end

%% Tabla de ranking por tipo de falta
fprintf(fid, '\n================================================================================\n');
fprintf(fid, 'RANKING: MEJOR MÉTODO POR TIPO DE FALTA (promedio sobre ventanas/convertidores)\n');
fprintf(fid, '================================================================================\n\n');

fprintf(fid, '%-10s | %-10s | %-10s | %-10s | %-10s\n', ...
    'Falta', '1º Lugar', 'Precisión', '2º Lugar', 'Precisión');
fprintf(fid, repmat('-', 1, 65));
fprintf(fid, '\n');

for f = 1:n_faults
    fault = fault_types{f};
    method_scores = zeros(1, length(methods));

    % Calcular promedio por método
    for m = 1:length(methods)
        scores = [];
        for w = 1:length(windows)
            for c = 1:length(converters)
                idx = find_result(all_results, fault, methods{m}, windows{w}, converters{c});
                if ~isempty(idx)
                    scores = [scores, all_results{idx}.metrics.fault_accuracy * 100];
                end
            end
        end
        method_scores(m) = mean(scores, 'omitnan');
    end

    % Ordenar
    [sorted_scores, sort_idx] = sort(method_scores, 'descend');

    fprintf(fid, '%-10s | %-10s | %9.1f%% | %-10s | %9.1f%%\n', ...
        fault, methods{sort_idx(1)}, sorted_scores(1), ...
        methods{sort_idx(2)}, sorted_scores(2));
end

%% Resumen global por método
fprintf(fid, '\n================================================================================\n');
fprintf(fid, 'RESUMEN GLOBAL POR MÉTODO\n');
fprintf(fid, '================================================================================\n\n');

fprintf(fid, '%-10s | %-15s | %-15s | %-15s\n', ...
    'Método', 'Precisión Media', 'Mejor Caso', 'Peor Caso');
fprintf(fid, repmat('-', 1, 70));
fprintf(fid, '\n');

for m = 1:length(methods)
    all_scores = [];

    for f = 1:n_faults
        for w = 1:length(windows)
            for c = 1:length(converters)
                idx = find_result(all_results, fault_types{f}, methods{m}, windows{w}, converters{c});
                if ~isempty(idx)
                    all_scores = [all_scores, all_results{idx}.metrics.fault_accuracy * 100];
                end
            end
        end
    end

    fprintf(fid, '%-10s | %14.1f%% | %14.1f%% | %14.1f%%\n', ...
        methods{m}, mean(all_scores, 'omitnan'), max(all_scores), min(all_scores));
end

fclose(fid);
fprintf('✓ Tabla formateada guardada en: %s\n', txt_file);

end

function idx = find_result(all_results, fault, method, window, converter)
% Buscar índice de resultado específico

for i = 1:length(all_results)
    if strcmp(all_results{i}.fault_type, fault) && ...
       strcmp(all_results{i}.method, method) && ...
       strcmp(all_results{i}.window_name, window) && ...
       strcmp(all_results{i}.converter, converter)
        idx = i;
        return;
    end
end
idx = [];
end
