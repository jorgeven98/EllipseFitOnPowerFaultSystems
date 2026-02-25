%% comparacion_metodos_ventana.m
% Comparacion de la efectividad de los tres clasificadores de fallas
% en funcion del tamano de la ventana de datos.
%
% Clasificadores:
%   1) pe3d_fault_classifier       (Polarization Ellipse 3D)
%   2) clarke_fault_classifier (Clarke + Ellipse)
%   3) ga_fault_classifier (Geometric Algebra + Ellipse)
%
% Genera multiples escenarios de falla realistas (modelo de red electrica
% con impedancias de Thevenin, linea, falta y tierra) variando:
%   - Tipo de falta (11 tipos)
%   - Distancia de falta (m)
%   - Impedancia de falta (Zf)
%   - Tamano de ventana (fraccion de ciclo)
%
% Metodologia de generacion de fallas: identica a prueba.m

clear; clc; close all;

% Agregar rutas a clasificadores
addpath(fullfile('..', 'classifiers'));

fprintf('================================================================\n');
fprintf('  COMPARACION DE CLASIFICADORES vs. TAMANO DE VENTANA DE DATOS\n');
fprintf('================================================================\n\n');

%% =====================================================================
%  PARAMETROS DEL SISTEMA
%  =====================================================================
fs = 10000;          % Frecuencia de muestreo (Hz)
f0 = 50;             % Frecuencia fundamental (Hz)
T  = 1/f0;           % Periodo (s)
samples_per_cycle = round(fs / f0);  % Muestras por ciclo

% Operador de rotacion 120 grados
a = exp(1j*2*pi/3);

%% =====================================================================
%  PARAMETROS DE RED (del paper, Seccion II-C) - Base de prueba.m
%  =====================================================================
Vth = 1.0;                         % Voltaje Thevenin (pu)
Zth = 0.01 + 1j*0.183;            % Impedancia Thevenin
Zg  = 1j*0.047;                   % Impedancia de tierra
ZL  = 0.01435 + 1j*0.03102;       % Impedancia de linea por km

%% =====================================================================
%  DEFINIR VARIACIONES DE PARAMETROS PARA GENERAR CASOS DE PRUEBA
%  =====================================================================

% Distancias de falta (fraccion de la linea)
m_values = [0.1, 0.3, 0.5, 0.7, 0.9];

% Impedancias de falta
Zf_values = [0.001, 0.01, 0.02, 0.05, 0.1];

% Tamanos de ventana (en fracciones de ciclo)
window_fractions = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0];
num_windows = length(window_fractions);

% Tipos de falta con sus etiquetas
fault_names = {'Normal', 'ABCG', 'AG', 'BG', 'CG', 'AB', 'BC', 'CA', 'ABG', 'BCG', 'CAG'};
num_fault_types = length(fault_names);

%% =====================================================================
%  EJECUTAR COMPARACION
%  =====================================================================

% Numero total de casos por tipo de falta y ventana
num_m  = length(m_values);
num_Zf = length(Zf_values);
cases_per_fault = num_m * num_Zf;  % 25 combinaciones por tipo

fprintf('Configuracion del experimento:\n');
fprintf('  Tipos de falta:       %d\n', num_fault_types);
fprintf('  Distancias de falta:  %d valores (%.1f a %.1f)\n', num_m, m_values(1), m_values(end));
fprintf('  Impedancias de falta: %d valores\n', num_Zf);
fprintf('  Tamanos de ventana:   %d valores (%.2f a %.1f ciclos)\n', num_windows, window_fractions(1), window_fractions(end));
fprintf('  Casos por tipo/ventana: %d\n', cases_per_fault);
fprintf('  Total de clasificaciones: %d\n\n', num_fault_types * cases_per_fault * num_windows * 3);

% Matrices de resultados: [num_fault_types x num_windows] por metodo
results_PE3D   = zeros(num_fault_types, num_windows);
results_Clarke = zeros(num_fault_types, num_windows);
results_GA     = zeros(num_fault_types, num_windows);

% Generar senal base larga (para cubrir ventana maxima)
max_cycles = max(window_fractions) + 1;
t_full = (0 : 1/fs : max_cycles*T - 1/fs)';
gen_signal = @(V_phasor) real(V_phasor * sqrt(2) * exp(1j*2*pi*f0*t_full));

% Verificar disponibilidad del toolbox GA
ga_available = true;
try
    gafulInit;
catch
    ga_available = false;
    fprintf('AVISO: GA-FuL Toolbox no disponible. El metodo GA se omitira.\n\n');
end

%% Bucle principal
fprintf('Ejecutando pruebas...\n');
fprintf('---------------------------------------------------------------------\n');

for w = 1:num_windows
    wfrac = window_fractions(w);
    window_samples = round(wfrac * samples_per_cycle);

    % Asegurar minimo de muestras
    if window_samples < 6
        window_samples = 6;
    end

    fprintf('\n>> Ventana: %.2f ciclos (%d muestras)\n', wfrac, window_samples);

    for ft = 1:num_fault_types
        correct_PE3D   = 0;
        correct_Clarke = 0;
        correct_GA     = 0;
        total_cases    = 0;

        for im = 1:num_m
            m_val = m_values(im);

            for iz = 1:num_Zf
                Zf = Zf_values(iz);
                total_cases = total_cases + 1;

                % Generar fasores (metodologia identica a prueba.m)
                [Va_ph, Vb_ph, Vc_ph] = generate_fault_phasors( ...
                    ft, Vth, a, Zth, ZL, Zf, Zg, m_val);

                % Generar senales temporales completas
                Va_full = gen_signal(Va_ph);
                Vb_full = gen_signal(Vb_ph);
                Vc_full = gen_signal(Vc_ph);

                % Extraer ventana de datos
                start_idx = 1;
                end_idx = min(start_idx + window_samples - 1, length(t_full));
                V_abc = [Va_full(start_idx:end_idx), ...
                         Vb_full(start_idx:end_idx), ...
                         Vc_full(start_idx:end_idx)];

                expected = fault_names{ft};

                % ---- Clasificador PE3D ----
                try
                    [res_PE3D, ~] = pe3d_fault_classifier(V_abc);
                    res_PE3D_norm = normalize_fault_name(res_PE3D);
                    if strcmpi(res_PE3D_norm, expected)
                        correct_PE3D = correct_PE3D + 1;
                    end
                catch
                    % Error en clasificador: no se cuenta como acierto
                end

                % ---- Clasificador Clarke ----
                try
                    [res_Clarke, ~] = clarke_fault_classifier(V_abc);
                    res_Clarke_norm = normalize_fault_name(res_Clarke);
                    if strcmpi(res_Clarke_norm, expected)
                        correct_Clarke = correct_Clarke + 1;
                    end
                catch
                    % Error en clasificador
                end

                % ---- Clasificador GA ----
                if ga_available
                    try
                        [res_GA, ~] = ga_fault_classifier(V_abc, fs, f0);
                        res_GA_norm = normalize_fault_name(res_GA);
                        if strcmpi(res_GA_norm, expected)
                            correct_GA = correct_GA + 1;
                        end
                    catch
                        % Error en clasificador
                    end
                end

            end  % Zf
        end  % m

        % Calcular porcentajes
        results_PE3D(ft, w)   = 100 * correct_PE3D / total_cases;
        results_Clarke(ft, w) = 100 * correct_Clarke / total_cases;
        if ga_available
            results_GA(ft, w) = 100 * correct_GA / total_cases;
        else
            results_GA(ft, w) = NaN;
        end

    end  % fault type

    % Resumen de la ventana actual
    fprintf('  Acierto promedio - PE3D: %.1f%% | Clarke: %.1f%%', ...
        mean(results_PE3D(:, w)), mean(results_Clarke(:, w)));
    if ga_available
        fprintf(' | GA: %.1f%%', mean(results_GA(:, w)));
    end
    fprintf('\n');
end

fprintf('\n---------------------------------------------------------------------\n');
fprintf('Pruebas completadas.\n\n');

%% =====================================================================
%  RESULTADOS: TABLA DETALLADA POR TIPO DE FALTA Y VENTANA
%  =====================================================================

fprintf('================================================================\n');
fprintf('  TABLA 1: PORCENTAJE DE ACIERTO POR TIPO DE FALTA Y VENTANA\n');
fprintf('================================================================\n\n');

methods_names = {'PE3D', 'Clarke', 'GA'};
all_results = cat(3, results_PE3D, results_Clarke, results_GA);

for met = 1:3
    if met == 3 && ~ga_available
        fprintf('\n[Metodo GA no disponible - omitido]\n');
        continue;
    end

    fprintf('--- Metodo: %s ---\n', methods_names{met});
    fprintf('%-8s', 'Falta');
    for w = 1:num_windows
        fprintf(' %6.2fc', window_fractions(w));
    end
    fprintf('  | Media\n');
    fprintf('%s\n', repmat('-', 1, 8 + num_windows*7 + 9));

    for ft = 1:num_fault_types
        fprintf('%-8s', fault_names{ft});
        for w = 1:num_windows
            val = all_results(ft, w, met);
            if isnan(val)
                fprintf('   N/A ');
            else
                fprintf(' %5.1f%%', val);
            end
        end
        row_mean = mean(all_results(ft, :, met), 'omitnan');
        fprintf('  | %5.1f%%', row_mean);
        fprintf('\n');
    end

    fprintf('%s\n', repmat('-', 1, 8 + num_windows*7 + 9));
    fprintf('%-8s', 'MEDIA');
    for w = 1:num_windows
        col_vals = all_results(:, w, met);
        fprintf(' %5.1f%%', mean(col_vals, 'omitnan'));
    end
    all_vals = all_results(:, :, met);
    fprintf('  | %5.1f%%', mean(all_vals(:), 'omitnan'));
    fprintf('\n\n');
end

%% =====================================================================
%  TABLA 2: RESUMEN COMPARATIVO (PROMEDIO GLOBAL POR VENTANA)
%  =====================================================================

fprintf('================================================================\n');
fprintf('  TABLA 2: RESUMEN COMPARATIVO - ACIERTO GLOBAL (%%) POR VENTANA\n');
fprintf('================================================================\n\n');

num_methods_active = 2 + ga_available;
fprintf('%-14s', 'Ventana');
for met = 1:3
    if met == 3 && ~ga_available, continue; end
    fprintf(' %10s', methods_names{met});
end
fprintf('\n');
fprintf('%s\n', repmat('-', 1, 14 + 10*num_methods_active));

for w = 1:num_windows
    fprintf('%-14s', sprintf('%.2f ciclos', window_fractions(w)));
    for met = 1:3
        if met == 3 && ~ga_available, continue; end
        col_vals = all_results(:, w, met);
        fprintf(' %9.1f%%', mean(col_vals, 'omitnan'));
    end
    fprintf('\n');
end

fprintf('%s\n', repmat('-', 1, 14 + 10*num_methods_active));
fprintf('%-14s', 'GLOBAL');
for met = 1:3
    if met == 3 && ~ga_available, continue; end
    all_vals = all_results(:, :, met);
    fprintf(' %9.1f%%', mean(all_vals(:), 'omitnan'));
end
fprintf('\n\n');

%% =====================================================================
%  TABLA 3: MEJOR METODO POR TIPO DE FALTA
%  =====================================================================

fprintf('================================================================\n');
fprintf('  TABLA 3: MEJOR METODO POR TIPO DE FALTA (ventana = 1 ciclo)\n');
fprintf('================================================================\n\n');

% Buscar columna de 1 ciclo
[~, idx_1cycle] = min(abs(window_fractions - 1.0));

fprintf('%-8s  %-8s  %-8s  %-8s  | Mejor\n', 'Falta', 'PE3D', 'Clarke', 'GA');
fprintf('%s\n', repmat('-', 1, 52));

for ft = 1:num_fault_types
    vals = [results_PE3D(ft, idx_1cycle), ...
            results_Clarke(ft, idx_1cycle), ...
            results_GA(ft, idx_1cycle)];

    if ~ga_available
        vals(3) = -1;
    end

    [~, best_idx] = max(vals);
    if ~ga_available && best_idx == 3
        [~, best_idx] = max(vals(1:2));
    end
    best_name = methods_names{best_idx};

    fprintf('%-8s  %6.1f%%   %6.1f%%   ', fault_names{ft}, vals(1), vals(2));
    if ga_available
        fprintf('%6.1f%%', vals(3));
    else
        fprintf('  N/A  ');
    end
    fprintf('   | %s\n', best_name);
end

fprintf('%s\n', repmat('-', 1, 52));

%% =====================================================================
%  GRAFICAS
%  =====================================================================

% Grafica 1: Acierto global vs. ventana
figure('Name', 'Comparacion Global', 'Position', [100 100 800 500]);
plot(window_fractions, mean(results_PE3D, 1), 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(window_fractions, mean(results_Clarke, 1), 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
if ga_available
    plot(window_fractions, mean(results_GA, 1, 'omitnan'), 'g-^', 'LineWidth', 2, 'MarkerSize', 8);
    legend('PE3D', 'Clarke', 'GA', 'Location', 'southeast');
else
    legend('PE3D', 'Clarke', 'Location', 'southeast');
end
hold off;
xlabel('Tamano de ventana (ciclos)');
ylabel('Acierto promedio (%)');
title('Efectividad de clasificacion vs. tamano de ventana');
grid on;
ylim([0 105]);
set(gca, 'FontSize', 12);

% Grafica 2: Mapa de calor por tipo de falta
figure('Name', 'Mapa de calor por metodo', 'Position', [150 150 900 500]);
subplot(1,3,1);
imagesc(results_PE3D);
colorbar;
caxis([0 100]);
set(gca, 'XTick', 1:num_windows, 'XTickLabel', ...
    arrayfun(@(x) sprintf('%.2f', x), window_fractions, 'UniformOutput', false));
set(gca, 'YTick', 1:num_fault_types, 'YTickLabel', fault_names);
xlabel('Ventana (ciclos)');
title('PE3D (%)');
colormap(gca, 'parula');
xtickangle(45);

subplot(1,3,2);
imagesc(results_Clarke);
colorbar;
caxis([0 100]);
set(gca, 'XTick', 1:num_windows, 'XTickLabel', ...
    arrayfun(@(x) sprintf('%.2f', x), window_fractions, 'UniformOutput', false));
set(gca, 'YTick', 1:num_fault_types, 'YTickLabel', fault_names);
xlabel('Ventana (ciclos)');
title('Clarke (%)');
colormap(gca, 'parula');
xtickangle(45);

subplot(1,3,3);
if ga_available
    imagesc(results_GA);
else
    imagesc(NaN(num_fault_types, num_windows));
end
colorbar;
caxis([0 100]);
set(gca, 'XTick', 1:num_windows, 'XTickLabel', ...
    arrayfun(@(x) sprintf('%.2f', x), window_fractions, 'UniformOutput', false));
set(gca, 'YTick', 1:num_fault_types, 'YTickLabel', fault_names);
xlabel('Ventana (ciclos)');
title('GA (%)');
colormap(gca, 'parula');
xtickangle(45);

sgtitle('Mapa de calor: Acierto por tipo de falta y ventana', 'FontSize', 14);

% Grafica 3: Acierto por categoria de falta
figure('Name', 'Acierto por categoria', 'Position', [200 200 900 400]);

% Agrupar por categoria
cat_names_list = {'Normal/ABCG', 'L-G (AG,BG,CG)', 'L-L (AB,BC,CA)', 'L-L-G (ABG,BCG,CAG)'};
cat_indices = {[1,2], [3,4,5], [6,7,8], [9,10,11]};

for c = 1:4
    subplot(1,4,c);
    idx = cat_indices{c};
    plot(window_fractions, mean(results_PE3D(idx,:), 1), 'b-o', 'LineWidth', 1.5);
    hold on;
    plot(window_fractions, mean(results_Clarke(idx,:), 1), 'r-s', 'LineWidth', 1.5);
    if ga_available
        plot(window_fractions, mean(results_GA(idx,:), 1, 'omitnan'), 'g-^', 'LineWidth', 1.5);
    end
    hold off;
    xlabel('Ventana (ciclos)');
    ylabel('Acierto (%)');
    title(cat_names_list{c});
    grid on;
    ylim([0 105]);
    if c == 4
        if ga_available
            legend('PE3D', 'Clarke', 'GA', 'Location', 'southeast');
        else
            legend('PE3D', 'Clarke', 'Location', 'southeast');
        end
    end
end

sgtitle('Efectividad por categoria de falta', 'FontSize', 14);

fprintf('================================================================\n');
fprintf('  ANALISIS COMPLETADO\n');
fprintf('  Resultados guardados en workspace:\n');
fprintf('    results_PE3D, results_Clarke, results_GA\n');
fprintf('    fault_names, window_fractions\n');
fprintf('================================================================\n');

%% =====================================================================
%  FUNCIONES LOCALES (deben ir al final del script)
%  =====================================================================

function canonical = normalize_fault_name(raw_name)
% Convierte cualquier nombre de falta al formato canonico
    mapping = {
        % PE3D names (ya canonicos)
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
        % Clarke names
        'Balanced',     'Normal';
        '3-Phase',      'ABCG';
        'L-G-A',        'AG';
        'L-G-B',        'BG';
        'L-G-C',        'CG';
        'L-L-AB',       'AB';
        'L-L-BC',       'BC';
        'L-L-CA',       'CA';
        'L-L-G-AB',     'ABG';
        'L-L-G-BC',     'BCG';
        'L-L-G-CA',     'CAG';
        'L-L-G',        'LLG_unknown';
        'L-L',          'LL_unknown';
        % GA names
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

function [Va_ph, Vb_ph, Vc_ph] = generate_fault_phasors(fault_id, Vth, a, Zth, ZL, Zf, Zg, m)
% Genera fasores de voltaje para cada tipo de falta
% usando el modelo de red electrica de prueba.m

    switch fault_id
        case 1  % Normal
            Va_ph = Vth;
            Vb_ph = Vth * a^2;
            Vc_ph = Vth * a;

        case 2  % ABCG (Simetrica)
            h = (m*ZL + Zf) / (Zth + m*ZL + Zf);
            Va_ph = h * Vth;
            Vb_ph = h * Vth * a^2;
            Vc_ph = h * Vth * a;

        case 3  % AG
            h = (m*ZL + Zf) / (Zth + m*ZL + Zf + 3*Zg);
            Va_ph = (1 - h) * Vth;
            Vb_ph = Vth * a^2;
            Vc_ph = Vth * a;

        case 4  % BG
            h = (m*ZL + Zf) / (Zth + m*ZL + Zf + 3*Zg);
            Va_ph = Vth;
            Vb_ph = (1 - h) * Vth * a^2;
            Vc_ph = Vth * a;

        case 5  % CG
            h = (m*ZL + Zf) / (Zth + m*ZL + Zf + 3*Zg);
            Va_ph = Vth;
            Vb_ph = Vth * a^2;
            Vc_ph = (1 - h) * Vth * a;

        case 6  % AB (linea-linea)
            h = (m*ZL + Zf) / (Zth + m*ZL + Zf);
            Va_ph = Vth - h*Vth*(1 - a^2)/2;
            Vb_ph = Vth*a^2 - h*Vth*(a^2 - 1)/2;
            Vc_ph = Vth * a;

        case 7  % BC (linea-linea)
            h = (m*ZL + Zf) / (Zth + m*ZL + Zf);
            Va_ph = Vth;
            Vb_ph = Vth*a^2 - h*Vth*(a^2 - a)/2;
            Vc_ph = Vth*a - h*Vth*(a - a^2)/2;

        case 8  % CA (linea-linea)
            h = (m*ZL + Zf) / (Zth + m*ZL + Zf);
            Va_ph = Vth - h*Vth*(1 - a)/2;
            Vb_ph = Vth * a^2;
            Vc_ph = Vth*a - h*Vth*(a - 1)/2;

        case 9  % ABG (doble linea a tierra)
            Z1 = Zth + m*ZL;
            Z2 = Zth + m*ZL;
            Z0 = Zth + m*ZL + 3*Zg;
            Zeq = Z1 + (Z2 * Z0)/(Z2 + Z0);
            If1 = Vth / (Zeq + Zf);
            Va_ph = Vth - If1*Z1;
            Vb_ph = Vth*a^2 - If1*Z1*a^2 - If1*(Z2*Z0/(Z2+Z0))*(a^2-1);
            Vc_ph = Vth*a;

        case 10  % BCG (doble linea a tierra)
            Z1 = Zth + m*ZL;
            Z2 = Zth + m*ZL;
            Z0 = Zth + m*ZL + 3*Zg;
            Zeq = Z1 + (Z2 * Z0)/(Z2 + Z0);
            If1 = Vth / (Zeq + Zf);
            Va_ph = Vth;
            Vb_ph = Vth*a^2 - If1*Z1*a^2 - If1*(Z2*Z0/(Z2+Z0))*(a^2-a);
            Vc_ph = Vth*a - If1*Z1*a - If1*(Z2*Z0/(Z2+Z0))*(a-a^2);

        case 11  % CAG (doble linea a tierra)
            Z1 = Zth + m*ZL;
            Z2 = Zth + m*ZL;
            Z0 = Zth + m*ZL + 3*Zg;
            Zeq = Z1 + (Z2 * Z0)/(Z2 + Z0);
            If1 = Vth / (Zeq + Zf);
            Va_ph = Vth - If1*Z1 - If1*(Z2*Z0/(Z2+Z0))*(1-a);
            Vb_ph = Vth*a^2;
            Vc_ph = Vth*a - If1*Z1*a - If1*(Z2*Z0/(Z2+Z0))*(a-1);

        otherwise
            error('Tipo de falta no reconocido: %d', fault_id);
    end
end
