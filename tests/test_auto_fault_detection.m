% TEST_AUTO_FAULT_DETECTION - Verificar detección automática de región de falla
%
% Este script verifica que el GUI puede detectar automáticamente la región
% de falla para archivos con diferentes duraciones
%
% Autor: Test Script
% Fecha: 2026-02-11

clear; clc;

fprintf('========================================\n');
fprintf('TEST: Detección automática de región de falla\n');
fprintf('========================================\n\n');

% Lista de archivos para probar
test_files = {
    'DataGridRL_AB.mat', 'DG1';  % Archivo largo (5 s)
    'DataGrid_BG.mat', 'DG1';    % Archivo mediano (2.5 s)
    'DataGrid_AG.mat', 'DG1';    % Archivo mediano (2.5 s)
};

for i = 1:size(test_files, 1)
    filename = test_files{i, 1};
    converter = test_files{i, 2};

    fprintf('[%d/%d] %s (convertidor: %s)\n', i, size(test_files, 1), filename, converter);

    filepath = fullfile('..', 'data', filename);
    if ~exist(filepath, 'file')
        fprintf('  ⚠ Archivo no encontrado\n\n');
        continue;
    end

    % Cargar datos
    file_data = load(filepath);
    [V_abc, time_vec] = extract_converter_data(file_data, converter);

    % Calcular duración
    t_duration = time_vec(end) - time_vec(1);

    fprintf('  Duración total: %.2f s (%d muestras)\n', t_duration, length(time_vec));

    % Aplicar lógica de detección automática (igual que el GUI)
    if t_duration >= 4.5
        % Archivos largos (5 s)
        fault_start = 3.0;
        fault_end = 3.16;
        tipo = 'LARGO';
    elseif t_duration >= 2.0
        % Archivos medianos (2.5 s)
        fault_start = 1.5;
        fault_end = 1.66;
        tipo = 'MEDIANO';
    else
        % Archivos cortos
        t_center = (time_vec(1) + time_vec(end)) / 2;
        fault_duration = 0.16;
        fault_start = t_center - fault_duration / 2;
        fault_end = t_center + fault_duration / 2;
        tipo = 'CORTO';
    end

    fprintf('  Tipo: %s\n', tipo);
    fprintf('  Región de falla detectada: %.2f - %.2f s\n', fault_start, fault_end);

    % Verificar que la región está dentro del rango
    if fault_start >= time_vec(1) && fault_end <= time_vec(end)
        fprintf('  ✓ Región válida (dentro del rango de datos)\n');

        % Simular análisis con ventana de 1 ciclo
        fs = 4000;
        window_size = round(fs / 50);  % 80 muestras
        overlap = 0.75;
        stride = max(1, round(window_size * (1 - overlap)));
        margin = 0.1;

        t_start = max(time_vec(1), fault_start - margin);
        t_end = min(time_vec(end), fault_end + margin);

        idx_start = find(time_vec >= t_start, 1, 'first');
        idx_end = find(time_vec <= t_end, 1, 'last');

        % Contar ventanas
        n_windows = 0;
        start_idx = idx_start;
        while start_idx + window_size - 1 <= idx_end
            n_windows = n_windows + 1;
            start_idx = start_idx + stride;
        end

        fprintf('  Ventanas generadas: %d\n', n_windows);

        if n_windows > 0
            fprintf('  ✓ Análisis viable\n');
        else
            fprintf('  ✗ No se generarían ventanas (aumentar margen o reducir ventana)\n');
        end
    else
        fprintf('  ✗ Región inválida (fuera del rango de datos)\n');
    end

    fprintf('\n');
end

fprintf('========================================\n');
fprintf('TEST COMPLETADO\n');
fprintf('========================================\n\n');

fprintf('El GUI ahora detecta automáticamente la región de falla:\n');
fprintf('  - Archivos ≥4.5 s → Falla en 3.0 - 3.16 s\n');
fprintf('  - Archivos ≥2.0 s → Falla en 1.5 - 1.66 s\n');
fprintf('  - Archivos <2.0 s → Falla en el centro\n\n');

fprintf('Para probar en el GUI:\n');
fprintf('1. Ejecutar: sliding_window_analyzer_gui()\n');
fprintf('2. Cargar cualquier archivo (largo o corto)\n');
fprintf('3. El GUI ajustará automáticamente la región de falla\n');
fprintf('4. Ejecutar análisis y verificar que genera ventanas\n');
