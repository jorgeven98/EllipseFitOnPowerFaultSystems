%% SMOKE_TEST — verifica que el repo corre tal cual se clona
% Ejecutar desde la raíz del repo. No requiere toolboxes ni librerías GA:
% si GA-FuL no está instalado, ga_fault_classifier delega automáticamente
% en el port MATLAB-puro (paridad verificada 3696/3696 ventanas).
clear; clc;
addpath('classifiers'); addpath('utils'); addpath('core');

fprintf('GA-FuL instalado: %s\n', string(~isempty(which('gafulInit'))));

data = load(fullfile('data', 'DataGridRL_AB.mat'));
[V, t] = extract_converter_data(data, 'Pex');
i0 = find(t >= 3.05, 1);
Vw = V(i0:i0+19, :);                       % ventana de 1/4 de ciclo (20 muestras @4 kHz)

[ft, p] = ga_fault_classifier(Vw, 4000, 50);
fprintf('Clasificación ventana de falta A-B: %s (roundness=%.3f, theta=%.1f deg)\n', ...
    ft, p.roundness, p.angle_deg);
assert(strcmp(ft, 'A-B'), 'smoke test FALLIDO: se esperaba A-B, se obtuvo %s', ft);

% ventana sana
iH = find(t >= 2.5, 1);
ftH = ga_fault_classifier(V(iH:iH+79, :), 4000, 50);   % 1 ciclo para estabilidad
fprintf('Clasificación ventana sana (1 ciclo): %s\n', ftH);

fprintf('\nsmoke_test: OK — el repositorio funciona en esta instalación.\n');
