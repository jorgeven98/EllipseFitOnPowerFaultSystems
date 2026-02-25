%% TEST_GA_ELLIPSE_FAULT_CLASSIFIER
% Script de prueba para la función ga_fault_classifier.m
% Genera diferentes tipos de faltas y verifica la clasificación
%
% Autor: Test Script
% Fecha: 2026-01-22

clear; clc; close all;

% Agregar rutas a clasificadores y utilidades
addpath(fullfile('..', 'classifiers'));
addpath(fullfile('..', 'utils'));

%% PARÁMETROS DE SIMULACIÓN
fs = 10000;                 % Frecuencia de muestreo (Hz)
f = 50;                     % Frecuencia del sistema (Hz)
T = 1/f;                    % Periodo
N_cycles = 2;               % Número de ciclos a simular
t = 0:1/fs:N_cycles*T-1/fs; % Vector de tiempo
N = length(t);

% Voltaje nominal (fase-neutro)
V_nom = 1.0;  % Normalizado a 1 p.u.

fprintf('=== TEST GA_ELLIPSE_FAULT_CLASSIFIER ===\n\n');
fprintf('Parámetros de simulación:\n');
fprintf('  Frecuencia: %d Hz\n', f);
fprintf('  Frecuencia de muestreo: %d Hz\n', fs);
fprintf('  Ciclos: %d\n', N_cycles);
fprintf('  Puntos: %d\n\n', N);

% Función auxiliar para verificar restricción de Kirchhoff
verify_kirchhoff = @(V_abc) max(abs(sum(V_abc, 2))) < 1e-10;

%% TEST 1: SISTEMA BALANCEADO (A-B-C)
fprintf('TEST 1: Sistema Balanceado (A-B-C)\n');
fprintf('------------------------------------\n');

Va = V_nom * sin(2*pi*f*t);
Vb = V_nom * sin(2*pi*f*t - 2*pi/3);
Vc = V_nom * sin(2*pi*f*t + 2*pi/3);

V_abc_balanced = [Va', Vb', Vc'];

% Verificar restricción de Kirchhoff
if verify_kirchhoff(V_abc_balanced)
    fprintf('  ✓ Restricción de Kirchhoff satisfecha (Va+Vb+Vc=0)\n');
else
    fprintf('  ✗ ADVERTENCIA: Suma de tensiones no es cero!\n');
end

try
    [fault_type, params] = ga_fault_classifier(V_abc_balanced, fs, f);
    fprintf('  Tipo detectado: %s\n', fault_type);
    fprintf('  Bivector: [%.4f, %.4f, %.4f]\n', params.bivector);
    fprintf('  Semiejes: a=%.4f, b=%.4f\n', params.semi_major, params.semi_minor);
    fprintf('  Ángulo: %.2f°\n', params.angle_deg);
    fprintf('  Redondez: %.4f\n', params.roundness);
    if params.is_line
        fprintf('  ⚠ Curva clasificada como LÍNEA RECTA\n');
    end
    fprintf('  ✓ Test 1 completado\n\n');
catch ME
    fprintf('  ✗ Error: %s\n\n', ME.message);
end

%% TEST 2: FALTA MONOFÁSICA A TIERRA (A-G)
fprintf('TEST 2: Falta Monofásica A-G\n');
fprintf('-----------------------------\n');

Va_ag = 0.1 * V_nom * sin(2*pi*f*t);  % Fase A colapsa
Vb_ag = V_nom * sin(2*pi*f*t - 2*pi/3);
Vc_ag = V_nom * sin(2*pi*f*t + 2*pi/3);

V_abc_ag = [Va_ag', Vb_ag', Vc_ag'];

% Verificar restricción de Kirchhoff
if verify_kirchhoff(V_abc_ag)
    fprintf('  ✓ Restricción de Kirchhoff satisfecha (Va+Vb+Vc=0)\n');
else
    fprintf('  ✗ ADVERTENCIA: Suma de tensiones = %.2e (no cero)\n', max(abs(sum(V_abc_ag, 2))));
end

try
    [fault_type, params] = ga_fault_classifier(V_abc_ag, fs, f);
    fprintf('  Tipo detectado: %s\n', fault_type);
    fprintf('  Bivector: [%.4f, %.4f, %.4f]\n', params.bivector);
    fprintf('  Semiejes: a=%.4f, b=%.4f\n', params.semi_major, params.semi_minor);
    fprintf('  Ángulo: %.2f°\n', params.angle_deg);
    fprintf('  Redondez: %.4f\n', params.roundness);
    if params.is_line
        fprintf('  ⚠ Curva clasificada como LÍNEA RECTA\n');
    end
    fprintf('  ✓ Test 2 completado\n\n');
catch ME
    fprintf('  ✗ Error: %s\n\n', ME.message);
end

%% TEST 3: FALTA MONOFÁSICA B TIERRA (B-G)
fprintf('TEST 3: Falta Monofásica B-G\n');
fprintf('-----------------------------\n');

Va_bg = V_nom * sin(2*pi*f*t);
Vb_bg = 0.1 * V_nom * sin(2*pi*f*t - 2*pi/3);  % Fase B colapsa
Vc_bg = V_nom * sin(2*pi*f*t + 2*pi/3);

V_abc_bg = [Va_bg', Vb_bg', Vc_bg'];

% Verificar restricción de Kirchhoff
if verify_kirchhoff(V_abc_bg)
    fprintf('  ✓ Restricción de Kirchhoff satisfecha (Va+Vb+Vc=0)\n');
else
    fprintf('  ✗ ADVERTENCIA: Suma de tensiones = %.2e (no cero)\n', max(abs(sum(V_abc_bg, 2))));
end

try
    [fault_type, params] = ga_fault_classifier(V_abc_bg, fs, f);
    fprintf('  Tipo detectado: %s\n', fault_type);
    fprintf('  Bivector: [%.4f, %.4f, %.4f]\n', params.bivector);
    fprintf('  Semiejes: a=%.4f, b=%.4f\n', params.semi_major, params.semi_minor);
    fprintf('  Ángulo: %.2f°\n', params.angle_deg);
    fprintf('  Redondez: %.4f\n', params.roundness);
    if params.is_line
        fprintf('  ⚠ Curva clasificada como LÍNEA RECTA\n');
    end
    fprintf('  ✓ Test 3 completado\n\n');
catch ME
    fprintf('  ✗ Error: %s\n\n', ME.message);
end

%% TEST 4: FALTA MONOFÁSICA C TIERRA (C-G)
fprintf('TEST 4: Falta Monofásica C-G\n');
fprintf('-----------------------------\n');

Va_cg = V_nom * sin(2*pi*f*t);
Vb_cg = V_nom * sin(2*pi*f*t - 2*pi/3);
Vc_cg = 0.1 * V_nom * sin(2*pi*f*t + 2*pi/3);  % Fase C colapsa

V_abc_cg = [Va_cg', Vb_cg', Vc_cg'];

% Verificar restricción de Kirchhoff
if verify_kirchhoff(V_abc_cg)
    fprintf('  ✓ Restricción de Kirchhoff satisfecha (Va+Vb+Vc=0)\n');
else
    fprintf('  ✗ ADVERTENCIA: Suma de tensiones = %.2e (no cero)\n', max(abs(sum(V_abc_cg, 2))));
end

try
    [fault_type, params] = ga_fault_classifier(V_abc_cg, fs, f);
    fprintf('  Tipo detectado: %s\n', fault_type);
    fprintf('  Bivector: [%.4f, %.4f, %.4f]\n', params.bivector);
    fprintf('  Semiejes: a=%.4f, b=%.4f\n', params.semi_major, params.semi_minor);
    fprintf('  Ángulo: %.2f°\n', params.angle_deg);
    fprintf('  Redondez: %.4f\n', params.roundness);
    if params.is_line
        fprintf('  ⚠ Curva clasificada como LÍNEA RECTA\n');
    end
    fprintf('  ✓ Test 4 completado\n\n');
catch ME
    fprintf('  ✗ Error: %s\n\n', ME.message);
end

%% TEST 5: FALTA LÍNEA-LÍNEA A-B
fprintf('TEST 5: Falta Línea-Línea A-B\n');
fprintf('------------------------------\n');

% Durante falta A-B: Va = Vb, y Va + Vb + Vc = 0
% Por lo tanto: 2*Va + Vc = 0 => Vc = -2*Va
% Usamos Vab como voltaje común entre A y B
Vab_common = 0.5 * V_nom * sin(2*pi*f*t);  % Voltaje común de fases en cortocircuito
Va_ab = Vab_common;
Vb_ab = Vab_common;
Vc_ab = -2 * Vab_common;  % Asegurar que Va + Vb + Vc = 0

V_abc_ab = [Va_ab', Vb_ab', Vc_ab'];

% Verificar restricción de Kirchhoff
if verify_kirchhoff(V_abc_ab)
    fprintf('  ✓ Restricción de Kirchhoff satisfecha (Va+Vb+Vc=0)\n');
else
    fprintf('  ✗ ADVERTENCIA: Suma de tensiones = %.2e (no cero)\n', max(abs(sum(V_abc_ab, 2))));
end
[fault_type, params] = ga_fault_classifier(V_abc_ab, fs, f);
try
    [fault_type, params] = ga_fault_classifier(V_abc_ab, fs, f);
    fprintf('  Tipo detectado: %s\n', fault_type);
    fprintf('  Bivector: [%.4f, %.4f, %.4f]\n', params.bivector);
    fprintf('  Semiejes: a=%.4f, b=%.4f\n', params.semi_major, params.semi_minor);
    fprintf('  Ángulo: %.2f°\n', params.angle_deg);
    fprintf('  Redondez: %.4f\n', params.roundness);
    if params.is_line
        fprintf('  ⚠ Curva clasificada como LÍNEA RECTA\n');
    end
    fprintf('  ✓ Test 5 completado\n\n');
catch ME
    fprintf('  ✗ Error: %s\n\n', ME.message);
end

%% TEST 6: FALTA LÍNEA-LÍNEA B-C
fprintf('TEST 6: Falta Línea-Línea B-C\n');
fprintf('------------------------------\n');

% Durante falta B-C: Vb = Vc, y Va + Vb + Vc = 0
% Por lo tanto: Va + 2*Vb = 0 => Va = -2*Vb
Vbc_common = 0.5 * V_nom * sin(2*pi*f*t - 2*pi/3);  % Voltaje común de fases en cortocircuito
Va_bc = -2 * Vbc_common;  % Asegurar que Va + Vb + Vc = 0
Vb_bc = Vbc_common;
Vc_bc = Vbc_common;

V_abc_bc = [Va_bc', Vb_bc', Vc_bc'];

% Verificar restricción de Kirchhoff
if verify_kirchhoff(V_abc_bc)
    fprintf('  ✓ Restricción de Kirchhoff satisfecha (Va+Vb+Vc=0)\n');
else
    fprintf('  ✗ ADVERTENCIA: Suma de tensiones = %.2e (no cero)\n', max(abs(sum(V_abc_bc, 2))));
end

try
    [fault_type, params] = ga_fault_classifier(V_abc_bc, fs, f);
    fprintf('  Tipo detectado: %s\n', fault_type);
    fprintf('  Bivector: [%.4f, %.4f, %.4f]\n', params.bivector);
    fprintf('  Semiejes: a=%.4f, b=%.4f\n', params.semi_major, params.semi_minor);
    fprintf('  Ángulo: %.2f°\n', params.angle_deg);
    fprintf('  Redondez: %.4f\n', params.roundness);
    if params.is_line
        fprintf('  ⚠ Curva clasificada como LÍNEA RECTA\n');
    end
    fprintf('  ✓ Test 6 completado\n\n');
catch ME
    fprintf('  ✗ Error: %s\n\n', ME.message);
end

%% TEST 7: FALTA LÍNEA-LÍNEA C-A
fprintf('TEST 7: Falta Línea-Línea C-A\n');
fprintf('------------------------------\n');

% Durante falta C-A: Vc = Va, y Va + Vb + Vc = 0
% Por lo tanto: 2*Va + Vb = 0 => Vb = -2*Va
Vca_common = 0.5 * V_nom * sin(2*pi*f*t + 2*pi/3);  % Voltaje común de fases en cortocircuito
Va_ca = Vca_common;
Vb_ca = -2 * Vca_common;  % Asegurar que Va + Vb + Vc = 0
Vc_ca = Vca_common;

V_abc_ca = [Va_ca', Vb_ca', Vc_ca'];

% Verificar restricción de Kirchhoff
if verify_kirchhoff(V_abc_ca)
    fprintf('  ✓ Restricción de Kirchhoff satisfecha (Va+Vb+Vc=0)\n');
else
    fprintf('  ✗ ADVERTENCIA: Suma de tensiones = %.2e (no cero)\n', max(abs(sum(V_abc_ca, 2))));
end

try
    [fault_type, params] = ga_fault_classifier(V_abc_ca, fs, f);
    fprintf('  Tipo detectado: %s\n', fault_type);
    fprintf('  Bivector: [%.4f, %.4f, %.4f]\n', params.bivector);
    fprintf('  Semiejes: a=%.4f, b=%.4f\n', params.semi_major, params.semi_minor);
    fprintf('  Ángulo: %.2f°\n', params.angle_deg);
    fprintf('  Redondez: %.4f\n', params.roundness);
    if params.is_line
        fprintf('  ⚠ Curva clasificada como LÍNEA RECTA\n');
    end
    fprintf('  ✓ Test 7 completado\n\n');
catch ME
    fprintf('  ✗ Error: %s\n\n', ME.message);
end

%% TEST 8: FALTA DOBLE LÍNEA-TIERRA A-B-G
fprintf('TEST 8: Falta Doble Línea-Tierra A-B-G\n');
fprintf('---------------------------------------\n');

Va_abg =  0.3*V_nom * sin(2*pi*f*t);
Vb_abg =  0.3* V_nom * sin(2*pi*f*t - 2*pi/3);
Vc_abg =  V_nom * sin(2*pi*f*t + 2*pi/3);

V_abc_abg = [Va_abg', Vb_abg', Vc_abg'];

% Verificar restricción de Kirchhoff
if verify_kirchhoff(V_abc_abg)
    fprintf('  ✓ Restricción de Kirchhoff satisfecha (Va+Vb+Vc=0)\n');
else
    fprintf('  ✗ ADVERTENCIA: Suma de tensiones = %.2e (no cero)\n', max(abs(sum(V_abc_abg, 2))));
end

try
    [fault_type, params] = ga_fault_classifier(V_abc_abg, fs, f);
    fprintf('  Tipo detectado: %s\n', fault_type);
    fprintf('  Bivector: [%.4f, %.4f, %.4f]\n', params.bivector);
    fprintf('  Semiejes: a=%.4f, b=%.4f\n', params.semi_major, params.semi_minor);
    fprintf('  Ángulo: %.2f°\n', params.angle_deg);
    fprintf('  Redondez: %.4f\n', params.roundness);
    if params.is_line
        fprintf('  ⚠ Curva clasificada como LÍNEA RECTA\n');
    end
    fprintf('  ✓ Test 8 completado\n\n');
catch ME
    fprintf('  ✗ Error: %s\n\n', ME.message);
end

%% VISUALIZACIÓN DE RESULTADOS
fprintf('=== GENERANDO GRÁFICAS ===\n\n');

figure('Name', 'Resultados de Clasificación', 'Position', [100, 100, 1200, 800]);

test_cases = {
    'Balanceado (A-B-C)', V_abc_balanced;
    'A-G', V_abc_ag;
    'B-G', V_abc_bg;
    'C-G', V_abc_cg;
    'A-B', V_abc_ab;
    'B-C', V_abc_bc;
    'C-A', V_abc_ca;
    'A-B-G', V_abc_abg;
};

for i = 1:size(test_cases, 1)
    try
        [fault_type, params] = ga_fault_classifier(test_cases{i, 2}, fs, f);

        subplot(3, 3, i);
        plot(params.x_T(:,1), params.x_T(:,2), 'b.', 'MarkerSize', 8);
        hold on;

        % Dibujar elipse ajustada si existe
        if params.semi_major > 0
            theta_ellipse = linspace(0, 2*pi, 100);
            x_ellipse = params.semi_major * cos(theta_ellipse);
            y_ellipse = params.semi_minor * sin(theta_ellipse);

            % Rotar elipse
            R_mat = [cos(deg2rad(params.angle_deg)), -sin(deg2rad(params.angle_deg));
                     sin(deg2rad(params.angle_deg)), cos(deg2rad(params.angle_deg))];
            ellipse_rot = R_mat * [x_ellipse; y_ellipse];

            % Trasladar al centro
            ellipse_rot(1,:) = ellipse_rot(1,:) + params.center(1);
            ellipse_rot(2,:) = ellipse_rot(2,:) + params.center(2);

            plot(ellipse_rot(1,:), ellipse_rot(2,:), 'r-', 'LineWidth', 2);
        end

        grid on;
        axis equal;
        title(sprintf('%s → %s', test_cases{i,1}, fault_type));
        xlabel('\sigma_1');
        ylabel('\sigma_2');
        legend('Datos rotados', 'Elipse ajustada', 'Location', 'best');

    catch ME
        subplot(3, 3, i);
        text(0.5, 0.5, sprintf('Error: %s', ME.message), ...
            'HorizontalAlignment', 'center');
        title(test_cases{i,1});
    end
end

fprintf('✓ Gráficas generadas\n\n');
fprintf('=== TESTS COMPLETADOS ===\n');
