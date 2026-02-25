function figs = generate_comparison_plots(all_results)
% GENERATE_COMPARISON_PLOTS Genera gráficos comparativos del experimento
%
% Input:
%   all_results: Cell array con todas las estructuras de resultados
%
% Output:
%   figs: Array de handles de figuras generadas
%
% Autor: Experiment Script
% Fecha: 2026-02-10

figs = [];

%% Extraer información para gráficos
n_results = length(all_results);
if n_results == 0
    warning('No hay resultados para graficar');
    return;
end

% Extraer datos organizados
fault_types = unique(cellfun(@(x) x.fault_type, all_results, 'UniformOutput', false));
converters = unique(cellfun(@(x) x.converter, all_results, 'UniformOutput', false));
window_names = unique(cellfun(@(x) x.window_name, all_results, 'UniformOutput', false));

% Extraer métodos (puede que algunos resultados no tengan método)
methods = {};
for i = 1:length(all_results)
    if isfield(all_results{i}, 'method')
        methods{end+1} = all_results{i}.method;
    end
end
if ~isempty(methods)
    methods = unique(methods);
else
    methods = {'GA'};  % Default si no hay método
end

%% 1. GRÁFICO DE BARRAS: Precisión por Ventana y Convertidor
fig1 = figure('Position', [100, 100, 1400, 800]);
figs = [figs, fig1];

n_faults = length(fault_types);
n_rows = ceil(sqrt(n_faults));
n_cols = ceil(n_faults / n_rows);

for f = 1:n_faults
    subplot(n_rows, n_cols, f);

    fault = fault_types{f};
    data_matrix = [];
    labels = {};

    for w = 1:length(window_names)
        row = [];
        for c = 1:length(converters)
            % Buscar resultado correspondiente
            idx = find_result_idx(all_results, fault, converters{c}, window_names{w});
            if ~isempty(idx)
                row = [row, all_results{idx}.metrics.fault_accuracy * 100];
            else
                row = [row, 0];
            end
        end
        data_matrix = [data_matrix; row];
        labels{end+1} = window_names{w};
    end

    bar(data_matrix);
    title(sprintf('Precisión: %s', fault), 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('Tamaño de Ventana');
    ylabel('Precisión (%)');
    set(gca, 'XTickLabel', labels);
    legend(converters, 'Location', 'best');
    grid on;
    ylim([0, 105]);
end

sgtitle('Precisión de Detección por Falta, Ventana y Convertidor', 'FontSize', 14, 'FontWeight', 'bold');

%% 2. GRÁFICO DE BARRAS: Tiempo de Detección
fig2 = figure('Position', [150, 150, 1400, 800]);
figs = [figs, fig2];

for f = 1:n_faults
    subplot(n_rows, n_cols, f);

    fault = fault_types{f};
    data_matrix = [];

    for w = 1:length(window_names)
        row = [];
        for c = 1:length(converters)
            idx = find_result_idx(all_results, fault, converters{c}, window_names{w});
            if ~isempty(idx) && ~isnan(all_results{idx}.metrics.detection_time)
                row = [row, all_results{idx}.metrics.detection_time * 1000]; % ms
            else
                row = [row, NaN];
            end
        end
        data_matrix = [data_matrix; row];
    end

    bar(data_matrix);
    title(sprintf('Tiempo Detección: %s', fault), 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('Tamaño de Ventana');
    ylabel('Tiempo (ms)');
    set(gca, 'XTickLabel', labels);
    legend(converters, 'Location', 'best');
    grid on;
end

sgtitle('Tiempo de Detección por Falta, Ventana y Convertidor', 'FontSize', 14, 'FontWeight', 'bold');

%% 3. HEATMAP: Precisión Global
fig3 = figure('Position', [200, 200, 1000, 600]);
figs = [figs, fig3];

for c = 1:length(converters)
    subplot(1, length(converters), c);

    % Matriz: filas=faltas, columnas=ventanas
    heatmap_data = zeros(length(fault_types), length(window_names));

    for f = 1:length(fault_types)
        for w = 1:length(window_names)
            idx = find_result_idx(all_results, fault_types{f}, converters{c}, window_names{w});
            if ~isempty(idx)
                heatmap_data(f, w) = all_results{idx}.metrics.fault_accuracy * 100;
            end
        end
    end

    imagesc(heatmap_data);
    colorbar;
    colormap(jet);
    caxis([0, 100]);

    title(sprintf('Precisión: %s', converters{c}), 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Tamaño de Ventana');
    ylabel('Tipo de Falta');
    set(gca, 'XTick', 1:length(window_names), 'XTickLabel', window_names);
    set(gca, 'YTick', 1:length(fault_types), 'YTickLabel', fault_types);
    xtickangle(45);

    % Agregar valores en celdas
    for f = 1:length(fault_types)
        for w = 1:length(window_names)
            text(w, f, sprintf('%.0f%%', heatmap_data(f, w)), ...
                'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold');
        end
    end
end

sgtitle('Mapa de Calor: Precisión de Detección', 'FontSize', 14, 'FontWeight', 'bold');

%% 4. COMPARACIÓN DG1 vs Pex
if length(converters) >= 2
    fig4 = figure('Position', [250, 250, 1200, 600]);
    figs = [figs, fig4];

    % Por cada ventana, comparar DG1 vs Pex
    for w = 1:length(window_names)
        subplot(2, 2, w);

        dg1_accuracies = [];
        pex_accuracies = [];

        for f = 1:length(fault_types)
            idx_dg1 = find_result_idx(all_results, fault_types{f}, converters{1}, window_names{w});
            idx_pex = find_result_idx(all_results, fault_types{f}, converters{2}, window_names{w});

            if ~isempty(idx_dg1)
                dg1_accuracies = [dg1_accuracies, all_results{idx_dg1}.metrics.fault_accuracy * 100];
            else
                dg1_accuracies = [dg1_accuracies, 0];
            end

            if ~isempty(idx_pex)
                pex_accuracies = [pex_accuracies, all_results{idx_pex}.metrics.fault_accuracy * 100];
            else
                pex_accuracies = [pex_accuracies, 0];
            end
        end

        plot(1:length(fault_types), dg1_accuracies, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
        hold on;
        plot(1:length(fault_types), pex_accuracies, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
        hold off;

        title(sprintf('Ventana: %s', window_names{w}), 'FontSize', 11, 'FontWeight', 'bold');
        xlabel('Tipo de Falta');
        ylabel('Precisión (%)');
        set(gca, 'XTick', 1:length(fault_types), 'XTickLabel', fault_types);
        xtickangle(45);
        legend(converters, 'Location', 'best');
        grid on;
        ylim([0, 105]);
    end

    sgtitle('Comparación DG1 vs Pex por Tamaño de Ventana', 'FontSize', 14, 'FontWeight', 'bold');
end

%% 5. COMPARACIÓN DE MÉTODOS (solo si hay múltiples métodos)
if length(methods) > 1
    fig5 = figure('Position', [300, 300, 1400, 800]);
    figs = [figs, fig5];

    % Por cada convertidor, comparar métodos
    for c = 1:length(converters)
        subplot(1, length(converters), c);

        % Matriz: filas=ventanas, columnas=métodos
        method_data = zeros(length(window_names), length(methods));

        for w = 1:length(window_names)
            for m = 1:length(methods)
                % Promediar precisión sobre todas las faltas
                accuracies = [];
                for f = 1:length(fault_types)
                    idx = find_result_idx_with_method(all_results, fault_types{f}, ...
                        converters{c}, window_names{w}, methods{m});
                    if ~isempty(idx)
                        accuracies = [accuracies, all_results{idx}.metrics.fault_accuracy * 100];
                    end
                end
                method_data(w, m) = mean(accuracies, 'omitnan');
            end
        end

        bar(method_data);
        title(sprintf('Convertidor: %s', converters{c}), 'FontSize', 12, 'FontWeight', 'bold');
        xlabel('Tamaño de Ventana');
        ylabel('Precisión Promedio (%)');
        set(gca, 'XTickLabel', window_names);
        xtickangle(45);
        legend(methods, 'Location', 'best');
        grid on;
        ylim([0, 105]);
    end

    sgtitle('Comparación de Métodos: Precisión Promedio por Ventana', 'FontSize', 14, 'FontWeight', 'bold');
end

%% 6. RANKING DE MÉTODOS POR FALTA (solo si hay múltiples métodos)
if length(methods) > 1
    fig6 = figure('Position', [350, 350, 1200, 800]);
    figs = [figs, fig6];

    n_rows = ceil(sqrt(n_faults));
    n_cols = ceil(n_faults / n_rows);

    for f = 1:n_faults
        subplot(n_rows, n_cols, f);

        fault = fault_types{f};
        method_scores = zeros(1, length(methods));

        % Calcular precisión promedio por método (sobre ventanas y convertidores)
        for m = 1:length(methods)
            scores = [];
            for w = 1:length(window_names)
                for c = 1:length(converters)
                    idx = find_result_idx_with_method(all_results, fault, ...
                        converters{c}, window_names{w}, methods{m});
                    if ~isempty(idx)
                        scores = [scores, all_results{idx}.metrics.fault_accuracy * 100];
                    end
                end
            end
            method_scores(m) = mean(scores, 'omitnan');
        end

        bar(method_scores);
        title(sprintf('Falta: %s', fault), 'FontSize', 11, 'FontWeight', 'bold');
        xlabel('Método');
        ylabel('Precisión Promedio (%)');
        set(gca, 'XTickLabel', methods);
        xtickangle(45);
        grid on;
        ylim([0, 105]);

        % Agregar valores sobre barras
        for m = 1:length(methods)
            text(m, method_scores(m) + 2, sprintf('%.1f%%', method_scores(m)), ...
                'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold');
        end
    end

    sgtitle('Ranking de Métodos por Tipo de Falta', 'FontSize', 14, 'FontWeight', 'bold');
end

fprintf('✓ %d gráficos comparativos generados\n', length(figs));

end

function idx = find_result_idx(all_results, fault, converter, window)
% Buscar índice de resultado específico (sin filtro de método)

for i = 1:length(all_results)
    if strcmp(all_results{i}.fault_type, fault) && ...
       strcmp(all_results{i}.converter, converter) && ...
       strcmp(all_results{i}.window_name, window)
        idx = i;
        return;
    end
end
idx = [];
end

function idx = find_result_idx_with_method(all_results, fault, converter, window, method)
% Buscar índice de resultado específico incluyendo método

for i = 1:length(all_results)
    % Verificar si tiene campo method
    result_method = 'GA';  % Default
    if isfield(all_results{i}, 'method')
        result_method = all_results{i}.method;
    end

    if strcmp(all_results{i}.fault_type, fault) && ...
       strcmp(all_results{i}.converter, converter) && ...
       strcmp(all_results{i}.window_name, window) && ...
       strcmp(result_method, method)
        idx = i;
        return;
    end
end
idx = [];
end
