function fig = plot_detection_timeline(result, fault_start, fault_end)
% PLOT_DETECTION_TIMELINE Genera gráfico de timeline de clasificación
%
% Input:
%   result: Estructura con resultados de un experimento
%   fault_start: Tiempo de inicio de falta (s)
%   fault_end: Tiempo de fin de falta (s)
%
% Output:
%   fig: Handle de la figura generada
%
% Autor: Experiment Script
% Fecha: 2026-02-10

% Crear figura
fig = figure('Position', [100, 100, 1200, 400]);

% Convertir clasificaciones a números para graficar
fault_types_list = {'A-B-C', 'A-G', 'B-G', 'C-G', 'A-B', 'B-C', 'C-A', ...
                    'A-B-G', 'B-C-G', 'C-A-G', 'A-B-C-G', 'Error', 'Unclassified'};
n_types = length(fault_types_list);

% Mapear clasificaciones a números
y_values = zeros(size(result.time_vector));
for i = 1:length(result.classifications)
    idx = find(strcmp(fault_types_list, result.classifications{i}), 1);
    if isempty(idx)
        % Si no se encuentra, usar "Unclassified"
        idx = find(strcmp(fault_types_list, 'Unclassified'));
    end
    y_values(i) = idx;
end

% Graficar
hold on;

% Áreas de fondo para regiones
% Pre-falta (azul claro)
fill([0, fault_start, fault_start, 0], [0, 0, n_types+1, n_types+1], ...
    [0.9, 0.9, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% Durante falta (rojo claro)
fill([fault_start, fault_end, fault_end, fault_start], [0, 0, n_types+1, n_types+1], ...
    [1, 0.9, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% Post-falta (verde claro)
fill([fault_end, max(result.time_vector), max(result.time_vector), fault_end], ...
    [0, 0, n_types+1, n_types+1], [0.9, 1, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

% Líneas verticales en inicio y fin de falta
plot([fault_start, fault_start], [0, n_types+1], 'r--', 'LineWidth', 2);
plot([fault_end, fault_end], [0, n_types+1], 'g--', 'LineWidth', 2);

% Graficar clasificaciones
plot(result.time_vector, y_values, 'bo-', 'MarkerSize', 4, 'LineWidth', 1.5);

% Marcar falta esperada
expected_idx = find(strcmp(fault_types_list, result.fault_type), 1);
if ~isempty(expected_idx)
    plot([0, max(result.time_vector)], [expected_idx, expected_idx], ...
        'k--', 'LineWidth', 1.5);
end

hold off;

% Configurar ejes
xlim([0, max(result.time_vector)]);
ylim([0, n_types+1]);
yticks(1:n_types);
yticklabels(fault_types_list);

xlabel('Tiempo (s)', 'FontSize', 12);
ylabel('Tipo de Falta Detectada', 'FontSize', 12);

% Título
title_str = sprintf('%s - %s - Ventana: %s\nPrecisión: %.1f%%, Detección: %.0f ms', ...
    result.fault_type, result.converter, result.window_name, ...
    result.metrics.fault_accuracy * 100, result.metrics.detection_time * 1000);
title(title_str, 'FontSize', 13, 'FontWeight', 'bold');

% Leyenda
legend({'Pre-falta', 'Durante falta', 'Post-falta', 'Inicio falta', 'Fin falta', ...
    'Clasificación', 'Falta esperada'}, 'Location', 'eastoutside', 'FontSize', 10);

grid on;
set(gca, 'FontSize', 10);

end
