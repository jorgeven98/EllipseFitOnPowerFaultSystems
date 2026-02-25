function sliding_window_analysis_gui()
% SLIDING_WINDOW_ANALYSIS_GUI - Interfaz para análisis de detección de fallas con ventana móvil
%
% Esta GUI permite:
%   - Seleccionar método de clasificación (GA, Clarke, PE3D)
%   - Seleccionar archivo de datos de falla
%   - Ajustar tamaño de ventana
%   - Visualizar el ratio de acierto al desplazar la ventana
%   - Ver la curva 3D y proyectada según el método
%   - Mostrar parámetros importantes de cada método
%
% Basado en: experiment_sliding_window_fault_detection.m
%
% Autor: Script generado para análisis de ventana deslizante
% Fecha: 2024

    %% Crear figura principal
    fig = figure('Name', 'Análisis de Ventana Deslizante - Detección de Fallas', ...
                 'NumberTitle', 'off', ...
                 'Position', [30, 30, 1500, 900], ...
                 'MenuBar', 'none', ...
                 'ToolBar', 'figure', ...
                 'Resize', 'on', ...
                 'Color', [0.94 0.94 0.94], ...
                 'CloseRequestFcn', @close_gui);

    %% Parámetros del sistema
    data.Ts = 2.5e-4;              % Periodo de muestreo (s)
    data.fs = 1/data.Ts;           % 4000 Hz
    data.f = 50;                   % Frecuencia del sistema (Hz)
    data.fault_start = 3.0;        % Inicio de falta (s)
    data.fault_end = 3.16;         % Fin de falta (s)
    data.samples_per_cycle = data.fs / data.f;  % 80 muestras/ciclo

    % Métodos disponibles
    data.methods = {'GA', 'Clarke', 'PE3D'};
    data.method_descriptions = {
        'Geometric Algebra - Bivector + Rotor', ...
        'Clarke Transform + PCA Ellipse', ...
        'Polarization Ellipse 3D (Stokes)'
    };

    % Convertidores disponibles
    data.converters = {'DG1', 'Pex'};
    data.converter_types = {'Grid Forming', 'Grid Following'};

    % Tamaños de ventana (en ciclos)
    data.window_cycles_options = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0];
    data.window_names = {'1/4 ciclo', '1/2 ciclo', '3/4 ciclo', '1 ciclo', '1.5 ciclos', '2 ciclos'};

    % Mapeo de archivos a tipos de falla
    data.fault_file_map = {
        'DataGridRL_AB.mat', 'AB', 'Bifásica A-B (con RL)';
        'DataGridRL_AC.mat', 'CA', 'Bifásica C-A (con RL)';
        'DataGridRL_BC.mat', 'BC', 'Bifásica B-C (con RL)';
        'DataGridRL_ABC.mat', 'ABCG', 'Trifásica (con RL)';
        'DataGrid_AG.mat', 'AG', 'Monofásica A-G';
        'DataGrid_BG.mat', 'BG', 'Monofásica B-G';
        'DataGrid_CG.mat', 'CG', 'Monofásica C-G';
        'DataGrid_ABG.mat', 'ABG', 'Doble línea A-B-G';
        'DataGrid_ACG.mat', 'CAG', 'Doble línea C-A-G';
        'DataGrid_BCG.mat', 'BCG', 'Doble línea B-C-G';
        'DataGrid_ABCG.mat', 'ABCG', 'Trifásica a tierra';
        'DataGrid_AB.mat', 'AB', 'Bifásica A-B (sin RL)';
        'DataGrid_BC.mat', 'BC', 'Bifásica B-C (sin RL)';
        'DataGrid_AC.mat', 'CA', 'Bifásica C-A (sin RL)';
    };

    % Estado inicial
    data.V_abc_full = [];
    data.time_full = [];
    data.current_file = '';
    data.expected_fault = '';
    data.classifications = {};
    data.time_centers = [];
    data.current_window_idx = 1;
    data.is_playing = false;
    data.play_timer = [];

    %% Panel de controles (izquierda)
    panel_ctrl = uipanel('Parent', fig, ...
                         'Title', 'Configuración', ...
                         'FontSize', 11, ...
                         'FontWeight', 'bold', ...
                         'Position', [0.005 0.01 0.18 0.98], ...
                         'BackgroundColor', [0.94 0.94 0.94]);

    % Selector de método
    uicontrol('Parent', panel_ctrl, 'Style', 'text', ...
              'String', 'Método de Clasificación:', ...
              'Units', 'normalized', 'Position', [0.05 0.94 0.9 0.025], ...
              'HorizontalAlignment', 'left', 'FontSize', 10, ...
              'BackgroundColor', [0.94 0.94 0.94]);

    data.popup_method = uicontrol('Parent', panel_ctrl, 'Style', 'popupmenu', ...
                                  'String', data.methods, ...
                                  'Units', 'normalized', 'Position', [0.05 0.905 0.9 0.035], ...
                                  'FontSize', 10, 'Value', 1, ...
                                  'Callback', @(~,~) on_method_change(fig));

    data.text_method_desc = uicontrol('Parent', panel_ctrl, 'Style', 'text', ...
                                      'String', data.method_descriptions{1}, ...
                                      'Units', 'normalized', 'Position', [0.05 0.87 0.9 0.03], ...
                                      'HorizontalAlignment', 'left', 'FontSize', 8, ...
                                      'ForegroundColor', [0.3 0.3 0.6], ...
                                      'BackgroundColor', [0.94 0.94 0.94]);

    % Separador
    uicontrol('Parent', panel_ctrl, 'Style', 'text', ...
              'String', '─────────────────────────', ...
              'Units', 'normalized', 'Position', [0.05 0.85 0.9 0.015], ...
              'BackgroundColor', [0.94 0.94 0.94]);

    % Selector de archivo
    uicontrol('Parent', panel_ctrl, 'Style', 'text', ...
              'String', 'Archivo de Datos:', ...
              'Units', 'normalized', 'Position', [0.05 0.82 0.9 0.025], ...
              'HorizontalAlignment', 'left', 'FontSize', 10, ...
              'BackgroundColor', [0.94 0.94 0.94]);

    data.popup_file = uicontrol('Parent', panel_ctrl, 'Style', 'popupmenu', ...
                                'String', data.fault_file_map(:, 3), ...
                                'Units', 'normalized', 'Position', [0.05 0.785 0.9 0.035], ...
                                'FontSize', 9, 'Value', 1);

    % Selector de convertidor
    uicontrol('Parent', panel_ctrl, 'Style', 'text', ...
              'String', 'Convertidor:', ...
              'Units', 'normalized', 'Position', [0.05 0.745 0.45 0.025], ...
              'HorizontalAlignment', 'left', 'FontSize', 10, ...
              'BackgroundColor', [0.94 0.94 0.94]);

    data.popup_converter = uicontrol('Parent', panel_ctrl, 'Style', 'popupmenu', ...
                                     'String', strcat(data.converters', ' (', data.converter_types', ')'), ...
                                     'Units', 'normalized', 'Position', [0.05 0.71 0.9 0.035], ...
                                     'FontSize', 9, 'Value', 1);

    % Botón cargar datos
    data.btn_load = uicontrol('Parent', panel_ctrl, 'Style', 'pushbutton', ...
                              'String', 'Cargar Datos', ...
                              'Units', 'normalized', 'Position', [0.05 0.66 0.9 0.04], ...
                              'FontSize', 10, 'FontWeight', 'bold', ...
                              'BackgroundColor', [0.3 0.5 0.8], ...
                              'ForegroundColor', 'white', ...
                              'Callback', @(~,~) load_data(fig));

    % Separador
    uicontrol('Parent', panel_ctrl, 'Style', 'text', ...
              'String', '─────────────────────────', ...
              'Units', 'normalized', 'Position', [0.05 0.63 0.9 0.015], ...
              'BackgroundColor', [0.94 0.94 0.94]);

    % Selector de ventana
    uicontrol('Parent', panel_ctrl, 'Style', 'text', ...
              'String', 'Tamaño de Ventana:', ...
              'Units', 'normalized', 'Position', [0.05 0.60 0.9 0.025], ...
              'HorizontalAlignment', 'left', 'FontSize', 10, ...
              'BackgroundColor', [0.94 0.94 0.94]);

    data.popup_window = uicontrol('Parent', panel_ctrl, 'Style', 'popupmenu', ...
                                  'String', data.window_names, ...
                                  'Units', 'normalized', 'Position', [0.05 0.565 0.9 0.035], ...
                                  'FontSize', 10, 'Value', 4, ...
                                  'Callback', @(~,~) on_window_change(fig));

    % Stride
    uicontrol('Parent', panel_ctrl, 'Style', 'text', ...
              'String', 'Solapamiento (%):', ...
              'Units', 'normalized', 'Position', [0.05 0.525 0.5 0.025], ...
              'HorizontalAlignment', 'left', 'FontSize', 10, ...
              'BackgroundColor', [0.94 0.94 0.94]);

    data.edit_overlap = uicontrol('Parent', panel_ctrl, 'Style', 'edit', ...
                                  'String', '75', ...
                                  'Units', 'normalized', 'Position', [0.6 0.525 0.35 0.03], ...
                                  'FontSize', 10);

    % Separador
    uicontrol('Parent', panel_ctrl, 'Style', 'text', ...
              'String', '─────────────────────────', ...
              'Units', 'normalized', 'Position', [0.05 0.49 0.9 0.015], ...
              'BackgroundColor', [0.94 0.94 0.94]);

    % Botón ejecutar análisis
    data.btn_analyze = uicontrol('Parent', panel_ctrl, 'Style', 'pushbutton', ...
                                 'String', '▶ EJECUTAR ANÁLISIS', ...
                                 'Units', 'normalized', 'Position', [0.05 0.44 0.9 0.045], ...
                                 'FontSize', 11, 'FontWeight', 'bold', ...
                                 'BackgroundColor', [0.2 0.6 0.2], ...
                                 'ForegroundColor', 'white', ...
                                 'Enable', 'off', ...
                                 'Callback', @(~,~) run_analysis(fig));

    % Separador
    uicontrol('Parent', panel_ctrl, 'Style', 'text', ...
              'String', '═════════════════════════', ...
              'Units', 'normalized', 'Position', [0.05 0.41 0.9 0.015], ...
              'BackgroundColor', [0.94 0.94 0.94]);

    %% Controles de navegación de ventana
    uicontrol('Parent', panel_ctrl, 'Style', 'text', ...
              'String', 'Navegación de Ventana:', ...
              'Units', 'normalized', 'Position', [0.05 0.38 0.9 0.025], ...
              'HorizontalAlignment', 'left', 'FontSize', 10, 'FontWeight', 'bold', ...
              'BackgroundColor', [0.94 0.94 0.94]);

    % Slider de posición
    data.slider_window = uicontrol('Parent', panel_ctrl, 'Style', 'slider', ...
                                   'Min', 1, 'Max', 100, 'Value', 1, ...
                                   'Units', 'normalized', 'Position', [0.05 0.34 0.9 0.03], ...
                                   'Enable', 'off', ...
                                   'Callback', @(src,~) on_slider_change(fig, src));

    % Botones de navegación
    data.btn_prev = uicontrol('Parent', panel_ctrl, 'Style', 'pushbutton', ...
                              'String', '◄◄', ...
                              'Units', 'normalized', 'Position', [0.05 0.295 0.2 0.035], ...
                              'FontSize', 10, 'Enable', 'off', ...
                              'Callback', @(~,~) navigate_window(fig, -10));

    data.btn_prev1 = uicontrol('Parent', panel_ctrl, 'Style', 'pushbutton', ...
                               'String', '◄', ...
                               'Units', 'normalized', 'Position', [0.27 0.295 0.15 0.035], ...
                               'FontSize', 10, 'Enable', 'off', ...
                               'Callback', @(~,~) navigate_window(fig, -1));

    data.btn_play = uicontrol('Parent', panel_ctrl, 'Style', 'pushbutton', ...
                              'String', '▶', ...
                              'Units', 'normalized', 'Position', [0.44 0.295 0.12 0.035], ...
                              'FontSize', 10, 'Enable', 'off', ...
                              'Callback', @(~,~) toggle_play(fig));

    data.btn_next1 = uicontrol('Parent', panel_ctrl, 'Style', 'pushbutton', ...
                               'String', '►', ...
                               'Units', 'normalized', 'Position', [0.58 0.295 0.15 0.035], ...
                               'FontSize', 10, 'Enable', 'off', ...
                               'Callback', @(~,~) navigate_window(fig, 1));

    data.btn_next = uicontrol('Parent', panel_ctrl, 'Style', 'pushbutton', ...
                              'String', '►►', ...
                              'Units', 'normalized', 'Position', [0.75 0.295 0.2 0.035], ...
                              'FontSize', 10, 'Enable', 'off', ...
                              'Callback', @(~,~) navigate_window(fig, 10));

    % Info de ventana actual
    data.text_window_info = uicontrol('Parent', panel_ctrl, 'Style', 'text', ...
                                      'String', 'Ventana: -/- | t = - ms', ...
                                      'Units', 'normalized', 'Position', [0.05 0.255 0.9 0.03], ...
                                      'HorizontalAlignment', 'center', 'FontSize', 9, ...
                                      'BackgroundColor', [0.94 0.94 0.94]);

    % Separador
    uicontrol('Parent', panel_ctrl, 'Style', 'text', ...
              'String', '═════════════════════════', ...
              'Units', 'normalized', 'Position', [0.05 0.23 0.9 0.015], ...
              'BackgroundColor', [0.94 0.94 0.94]);

    %% Panel de resultados
    panel_results = uipanel('Parent', panel_ctrl, ...
                            'Title', 'Resultados', ...
                            'FontSize', 10, ...
                            'FontWeight', 'bold', ...
                            'Position', [0.03 0.01 0.94 0.22], ...
                            'BackgroundColor', [1 1 1]);

    data.text_results = uicontrol('Parent', panel_results, 'Style', 'text', ...
                                  'String', 'Cargue datos y ejecute el análisis', ...
                                  'Units', 'normalized', 'Position', [0.02 0.02 0.96 0.94], ...
                                  'HorizontalAlignment', 'left', 'FontSize', 9, ...
                                  'FontName', 'Consolas', ...
                                  'BackgroundColor', [1 1 1]);

    %% Panel de visualización (derecha)
    panel_viz = uipanel('Parent', fig, ...
                        'Title', 'Visualización', ...
                        'FontSize', 11, ...
                        'FontWeight', 'bold', ...
                        'Position', [0.19 0.01 0.80 0.98], ...
                        'BackgroundColor', [1 1 1]);

    % Crear ejes para las gráficas
    % Fila 1: Señal completa con región de falla y timeline de clasificaciones
    data.ax_signal = axes('Parent', panel_viz, 'Position', [0.05 0.78 0.55 0.18]);
    data.ax_timeline = axes('Parent', panel_viz, 'Position', [0.65 0.78 0.32 0.18]);

    % Fila 2: Curva 3D y proyección 2D
    data.ax_3d = axes('Parent', panel_viz, 'Position', [0.05 0.42 0.42 0.32]);
    data.ax_2d = axes('Parent', panel_viz, 'Position', [0.55 0.42 0.42 0.32]);

    % Fila 3: Parámetros y métricas
    data.ax_params = axes('Parent', panel_viz, 'Position', [0.05 0.05 0.42 0.32]);
    data.ax_accuracy = axes('Parent', panel_viz, 'Position', [0.55 0.05 0.42 0.32]);

    %% Guardar datos en figura
    guidata(fig, data);

    %% Añadir rutas a utilidades y clasificadores
    addpath(fullfile('..', 'utils'));
    addpath(fullfile('..', 'classifiers'));
end

%% ========================================================================
%  CALLBACKS
%% ========================================================================

function on_method_change(fig)
    data = guidata(fig);
    idx = get(data.popup_method, 'Value');
    set(data.text_method_desc, 'String', data.method_descriptions{idx});

    % Si ya hay datos cargados, actualizar visualización
    if ~isempty(data.V_abc_full) && ~isempty(data.classifications)
        update_current_window_display(fig);
    end
end

function on_window_change(fig)
    data = guidata(fig);
    % Si hay datos cargados, sugerir re-ejecutar análisis
    if ~isempty(data.V_abc_full)
        set(data.text_results, 'String', 'Tamaño de ventana cambiado. Ejecute el análisis nuevamente.');
    end
end

function load_data(fig)
    data = guidata(fig);

    set(fig, 'Pointer', 'watch');
    drawnow limitrate;

    try
        % Obtener selecciones
        file_idx = get(data.popup_file, 'Value');
        conv_idx = get(data.popup_converter, 'Value');

        filename = data.fault_file_map{file_idx, 1};
        expected_fault = data.fault_file_map{file_idx, 2};
        converter_name = data.converters{conv_idx};

        % Cargar archivo
        filepath = fullfile('data', filename);
        if ~exist(filepath, 'file')
            error('Archivo no encontrado: %s', filepath);
        end

        file_data = load(filepath);

        % Extraer datos del convertidor
        [V_abc, time_vec] = extract_converter_data(file_data, converter_name);

        % Detectar automáticamente la región de falla basándose en la duración
        t_duration = time_vec(end) - time_vec(1);

        if t_duration >= 4.5
            % Archivos largos (5 s): usar región 3.0 - 3.16 s
            data.fault_start = 3.0;
            data.fault_end = 3.16;
        elseif t_duration >= 2.0
            % Archivos medianos (2.5 s): usar región 1.5 - 1.66 s
            data.fault_start = 1.5;
            data.fault_end = 1.66;
        else
            % Archivos cortos: usar el centro del archivo
            t_center = (time_vec(1) + time_vec(end)) / 2;
            fault_duration = 0.16;  % 160 ms (8 ciclos a 50 Hz)
            data.fault_start = t_center - fault_duration / 2;
            data.fault_end = t_center + fault_duration / 2;
        end

        % Guardar datos
        data.V_abc_full = V_abc;
        data.time_full = time_vec;
        data.current_file = filename;
        data.expected_fault = expected_fault;
        data.classifications = {};
        data.time_centers = [];

        % Habilitar botón de análisis
        set(data.btn_analyze, 'Enable', 'on');

        % Actualizar texto de resultados
        set(data.text_results, 'String', sprintf(['Datos cargados:\n' ...
            '  Archivo: %s\n' ...
            '  Falla esperada: %s\n' ...
            '  Convertidor: %s\n' ...
            '  Muestras: %d\n' ...
            '  Duración: %.2f s\n' ...
            '  Región de falla: %.2f - %.2f s\n\n' ...
            'Ejecute el análisis para procesar.'], ...
            filename, expected_fault, converter_name, ...
            size(V_abc, 1), time_vec(end) - time_vec(1), ...
            data.fault_start, data.fault_end));

        % Visualizar señal completa
        plot_full_signal(fig);

        guidata(fig, data);

    catch ME
        set(data.text_results, 'String', sprintf('Error cargando datos:\n%s', ME.message));
    end

    set(fig, 'Pointer', 'arrow');
end

function run_analysis(fig)
    data = guidata(fig);

    if isempty(data.V_abc_full)
        return;
    end

    set(fig, 'Pointer', 'watch');
    set(data.btn_analyze, 'Enable', 'off');
    drawnow limitrate;

    try
        % Obtener parámetros
        method_idx = get(data.popup_method, 'Value');
        window_idx = get(data.popup_window, 'Value');
        overlap = str2double(get(data.edit_overlap, 'String')) / 100;

        method_name = data.methods{method_idx};
        window_cycles = data.window_cycles_options(window_idx);
        window_size = round(window_cycles * data.samples_per_cycle);
        stride = max(1, round(window_size * (1 - overlap)));

        % Función clasificadora según método
        switch method_name
            case 'GA'
                classifier_func = @(V) ga_fault_classifier(V, data.fs, data.f);
            case 'Clarke'
                classifier_func = @(V) clarke_fault_classifier(V);
            case 'PE3D'
                classifier_func = @(V) pe3d_fault_classifier(V);
        end

        % Región de interés (centrada en la falla)
        margin = 0.1;  % 100 ms antes y después
        t_start = max(data.time_full(1), data.fault_start - margin);
        t_end = min(data.time_full(end), data.fault_end + margin);

        idx_start = find(data.time_full >= t_start, 1, 'first');
        idx_end = find(data.time_full <= t_end, 1, 'last');

        % Ventana móvil
        classifications = {};
        time_centers = [];
        params_array = {};
        correct_array = [];

        start_idx = idx_start;
        n_windows = 0;

        while start_idx + window_size - 1 <= idx_end
            end_idx = start_idx + window_size - 1;

            % Extraer ventana
            V_window = data.V_abc_full(start_idx:end_idx, :);

            % Clasificar
            try
                [fault_type, params] = classifier_func(V_window);
            catch
                fault_type = 'Error';
                params = struct();
            end

            % Tiempo central
            center_idx = round((start_idx + end_idx) / 2);
            time_center = data.time_full(center_idx);

            % Verificar si es correcto
            is_correct = check_classification_match(fault_type, data.expected_fault);

            % Guardar
            classifications{end+1} = fault_type;
            time_centers(end+1) = time_center;
            params_array{end+1} = params;
            correct_array(end+1) = is_correct;

            n_windows = n_windows + 1;
            start_idx = start_idx + stride;
        end

        % Verificar si se generaron ventanas
        if n_windows == 0
            % No hay ventanas suficientes
            data_duration = data.time_full(end) - data.time_full(1);
            region_duration = t_end - t_start;

            set(data.text_results, 'String', sprintf(['⚠ No se generaron ventanas\n\n' ...
                'Región de análisis:\n' ...
                '  t_start: %.3f s\n' ...
                '  t_end: %.3f s\n' ...
                '  Duración: %.3f s\n\n' ...
                'Datos disponibles:\n' ...
                '  Inicio: %.3f s\n' ...
                '  Fin: %.3f s\n' ...
                '  Duración: %.3f s\n\n' ...
                'Posibles causas:\n' ...
                '  - Falla fuera del rango\n' ...
                '  - Datos insuficientes\n' ...
                '  - Ventana muy grande'], ...
                t_start, t_end, region_duration, ...
                data.time_full(1), data.time_full(end), data_duration));

            set(data.btn_analyze, 'Enable', 'on');
            set(fig, 'Pointer', 'arrow');
            return;
        end

        % Guardar resultados
        data.classifications = classifications;
        data.time_centers = time_centers;
        data.params_array = params_array;
        data.correct_array = correct_array;
        data.window_size = window_size;
        data.idx_start_region = idx_start;
        data.idx_end_region = idx_end;

        % Calcular métricas
        accuracy = sum(correct_array) / length(correct_array) * 100;

        % En zona de falla
        fault_mask = time_centers >= data.fault_start & time_centers <= data.fault_end;
        if any(fault_mask)
            fault_accuracy = sum(correct_array(fault_mask)) / sum(fault_mask) * 100;
        else
            fault_accuracy = 0;
        end

        % Actualizar resultados
        set(data.text_results, 'String', sprintf(['Análisis completado:\n' ...
            '  Método: %s\n' ...
            '  Ventanas: %d\n' ...
            '  Ventana: %d muestras\n' ...
            '  Stride: %d muestras\n\n' ...
            'Precisión total: %.1f%%\n' ...
            'Precisión en falla: %.1f%%\n' ...
            'Falla esperada: %s'], ...
            method_name, n_windows, window_size, stride, ...
            accuracy, fault_accuracy, data.expected_fault));

        % Configurar controles de navegación
        set(data.slider_window, 'Min', 1, 'Max', max(2, n_windows), 'Value', 1, 'Enable', 'on');

        % Calcular SliderStep con protección contra divisiones por cero
        if n_windows > 1
            slider_step_small = 1 / (n_windows - 1);
            slider_step_large = min(10 / (n_windows - 1), 1);  % Asegurar que no exceda 1
        else
            slider_step_small = 0.1;  % Valor por defecto
            slider_step_large = 0.5;
        end
        set(data.slider_window, 'SliderStep', [slider_step_small, slider_step_large]);

        set([data.btn_prev, data.btn_prev1, data.btn_play, data.btn_next1, data.btn_next], 'Enable', 'on');

        data.current_window_idx = 1;
        guidata(fig, data);

        % Actualizar visualizaciones
        plot_full_signal(fig);
        plot_timeline(fig);
        update_current_window_display(fig);

    catch ME
        set(data.text_results, 'String', sprintf('Error en análisis:\n%s', ME.message));
    end

    set(data.btn_analyze, 'Enable', 'on');
    set(fig, 'Pointer', 'arrow');
end

function on_slider_change(fig, src)
    data = guidata(fig);
    if isempty(data.classifications)
        return;
    end

    idx = round(get(src, 'Value'));
    idx = max(1, min(idx, length(data.classifications)));
    data.current_window_idx = idx;
    guidata(fig, data);

    update_current_window_display(fig);
end

function navigate_window(fig, delta)
    data = guidata(fig);
    if isempty(data.classifications)
        return;
    end

    new_idx = data.current_window_idx + delta;
    new_idx = max(1, min(new_idx, length(data.classifications)));

    data.current_window_idx = new_idx;
    set(data.slider_window, 'Value', new_idx);
    guidata(fig, data);

    update_current_window_display(fig);
end

function toggle_play(fig)
    data = guidata(fig);

    if data.is_playing
        % Detener
        data.is_playing = false;
        if ~isempty(data.play_timer) && isvalid(data.play_timer)
            stop(data.play_timer);
            delete(data.play_timer);
        end
        data.play_timer = [];
        set(data.btn_play, 'String', '▶');
    else
        % Iniciar
        data.is_playing = true;
        set(data.btn_play, 'String', '❚❚');

        data.play_timer = timer('ExecutionMode', 'fixedRate', ...
                                'Period', 0.1, ...
                                'TimerFcn', @(~,~) play_step(fig));
        start(data.play_timer);
    end

    guidata(fig, data);
end

function play_step(fig)
    if ~ishandle(fig)
        return;
    end

    data = guidata(fig);

    if ~data.is_playing
        return;
    end

    new_idx = data.current_window_idx + 1;
    if new_idx > length(data.classifications)
        % Llegó al final, detener
        toggle_play(fig);
        return;
    end

    data.current_window_idx = new_idx;
    set(data.slider_window, 'Value', new_idx);
    guidata(fig, data);

    update_current_window_display(fig);
end

function close_gui(fig, ~)
    data = guidata(fig);

    % Detener timer si existe
    if isfield(data, 'play_timer') && ~isempty(data.play_timer) && isvalid(data.play_timer)
        stop(data.play_timer);
        delete(data.play_timer);
    end

    delete(fig);
end

%% ========================================================================
%  FUNCIONES DE VISUALIZACIÓN
%% ========================================================================

function plot_full_signal(fig)
    data = guidata(fig);

    if isempty(data.V_abc_full)
        return;
    end

    axes(data.ax_signal);
    cla;

    % Definir región de visualización (con margen)
    margin = 0.15;  % 150 ms antes y después
    t_start_view = max(data.time_full(1), data.fault_start - margin);
    t_end_view = min(data.time_full(end), data.fault_end + margin);

    % Extraer solo la región de interés para graficar
    idx_region = find(data.time_full >= t_start_view & data.time_full <= t_end_view);
    t = data.time_full(idx_region);
    V = data.V_abc_full(idx_region, :);

    % Reducir puntos SOLO si la región tiene demasiados puntos
    % (más inteligente: reducir solo si hay más de 5000 puntos en la región visible)
    N = length(t);
    if N > 5000
        step = ceil(N / 5000);
        t = t(1:step:end);
        V = V(1:step:end, :);
    end

    % Dibujar señales
    plot(t * 1000, V(:,1), 'r-', 'LineWidth', 0.8); hold on;
    plot(t * 1000, V(:,2), 'g-', 'LineWidth', 0.8);
    plot(t * 1000, V(:,3), 'b-', 'LineWidth', 0.8);

    % Marcar región de falla
    yl = ylim;
    patch([data.fault_start, data.fault_end, data.fault_end, data.fault_start] * 1000, ...
          [yl(1), yl(1), yl(2), yl(2)], ...
          [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

    % Marcar ventana actual si hay análisis
    if ~isempty(data.classifications) && isfield(data, 'current_window_idx')
        idx = data.current_window_idx;
        if idx <= length(data.time_centers)
            t_center = data.time_centers(idx);
            window_half = data.window_size / 2 / data.fs;
            t_start_win = t_center - window_half;
            t_end_win = t_center + window_half;

            patch([t_start_win, t_end_win, t_end_win, t_start_win] * 1000, ...
                  [yl(1), yl(1), yl(2), yl(2)], ...
                  [0.8 0.8 1], 'EdgeColor', 'b', 'LineWidth', 2, 'FaceAlpha', 0.3);
        end
    end

    xlabel('Tiempo (ms)');
    ylabel('Voltaje');
    title(sprintf('Señal Completa - %s', data.current_file), 'FontWeight', 'bold', 'Interpreter', 'none');
    legend('Va', 'Vb', 'Vc', 'Location', 'best');
    grid on;

    % Ajustar límites del eje X
    xlim([t_start_view * 1000, t_end_view * 1000]);
end

function plot_timeline(fig)
    data = guidata(fig);

    if isempty(data.classifications)
        return;
    end

    axes(data.ax_timeline);
    cla;

    t = data.time_centers * 1000;  % ms
    correct = data.correct_array;

    % Colorear por correcto/incorrecto
    scatter(t(correct == 1), ones(sum(correct == 1), 1), 40, [0 0.6 0], 'filled'); hold on;
    scatter(t(correct == 0), ones(sum(correct == 0), 1), 40, [0.8 0 0], 'filled');

    % Marcar región de falla
    yl = [0.5, 1.5];
    patch([data.fault_start, data.fault_end, data.fault_end, data.fault_start] * 1000, ...
          [yl(1), yl(1), yl(2), yl(2)], ...
          [1 0.9 0.9], 'EdgeColor', 'r', 'LineWidth', 1, 'FaceAlpha', 0.3);

    % Marcar ventana actual
    if isfield(data, 'current_window_idx') && data.current_window_idx <= length(t)
        xline(t(data.current_window_idx), 'b-', 'LineWidth', 2);
    end

    xlabel('Tiempo (ms)');
    ylabel('');
    title('Timeline de Clasificaciones', 'FontWeight', 'bold');
    ylim(yl);
    xlim([t(1) - 5, t(end) + 5]);
    set(gca, 'YTick', []);
    grid on;

    % Leyenda
    legend('Correcto', 'Incorrecto', 'Zona falla', 'Location', 'best');
end

function update_current_window_display(fig)
    data = guidata(fig);

    if isempty(data.classifications) || ~isfield(data, 'current_window_idx')
        return;
    end

    idx = data.current_window_idx;
    if idx > length(data.classifications)
        return;
    end

    % Actualizar info de ventana
    t_center = data.time_centers(idx);
    classification = data.classifications{idx};
    is_correct = data.correct_array(idx);
    params = data.params_array{idx};

    correct_str = 'CORRECTO';
    if ~is_correct
        correct_str = 'INCORRECTO';
    end

    set(data.text_window_info, 'String', ...
        sprintf('Ventana: %d/%d | t = %.1f ms | %s', ...
        idx, length(data.classifications), t_center * 1000, correct_str));

    % Extraer datos de la ventana actual
    window_half = data.window_size / 2;
    center_sample_idx = find(data.time_full >= t_center, 1, 'first');
    start_sample = max(1, center_sample_idx - window_half);
    end_sample = min(size(data.V_abc_full, 1), center_sample_idx + window_half - 1);

    V_window = data.V_abc_full(start_sample:end_sample, :);

    % Actualizar gráfico de señal (marcar ventana actual)
    plot_full_signal(fig);
    plot_timeline(fig);

    % Actualizar visualización 3D
    plot_3d_curve(fig, V_window, classification, is_correct);

    % Actualizar visualización 2D según método
    method_idx = get(data.popup_method, 'Value');
    plot_2d_projection(fig, V_window, params, data.methods{method_idx}, classification, is_correct);

    % Actualizar parámetros
    plot_parameters(fig, params, data.methods{method_idx}, classification);

    % Actualizar gráfico de precisión acumulada
    plot_accuracy_evolution(fig);

    drawnow limitrate;
end

function plot_3d_curve(fig, V_window, classification, is_correct)
    data = guidata(fig);

    axes(data.ax_3d);
    cla;

    % Reducir puntos si necesario
    N = size(V_window, 1);
    if N > 300
        step = ceil(N / 300);
        V_plot = V_window(1:step:end, :);
    else
        V_plot = V_window;
    end

    % Dibujar trayectoria 3D
    plot3(V_plot(:,1), V_plot(:,2), V_plot(:,3), 'b-', 'LineWidth', 1.5);
    hold on;

    % Marcar inicio y fin
    plot3(V_window(1,1), V_window(1,2), V_window(1,3), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
    plot3(V_window(end,1), V_window(end,2), V_window(end,3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

    grid on;
    xlabel('V_a'); ylabel('V_b'); zlabel('V_c');
    view(45, 30);
    axis equal;

    if is_correct
        title_color = [0 0.5 0];
    else
        title_color = [0.8 0 0];
    end

    title(sprintf('Trayectoria 3D: %s', classification), 'FontWeight', 'bold', 'Color', title_color);
end

function plot_2d_projection(fig, V_window, params, method_name, classification, is_correct)
    data = guidata(fig);

    axes(data.ax_2d);
    cla;

    switch method_name
        case 'Clarke'
            % Transformada Clarke
            if isfield(params, 'V_alpha') && isfield(params, 'V_beta')
                Va = params.V_alpha;
                Vb = params.V_beta;

                % Reducir puntos
                N = length(Va);
                if N > 200
                    step = ceil(N / 200);
                    Va = Va(1:step:end);
                    Vb = Vb(1:step:end);
                end

                plot(Va, Vb, 'b.', 'MarkerSize', 5);
                hold on;

                % Dibujar elipse ajustada
                if isfield(params, 'semi_major') && isfield(params, 'semi_minor') && isfield(params, 'angle_rad')
                    draw_ellipse_on_axes(0, 0, params.semi_major, params.semi_minor, params.angle_rad, 'r-', 2);
                end

                xlabel('V_\alpha'); ylabel('V_\beta');
            end

        case 'GA'
            % Datos rotados GA
            if isfield(params, 'x_T') && ~isempty(params.x_T)
                x_T = params.x_T;

                % Reducir puntos
                N = size(x_T, 1);
                if N > 200
                    step = ceil(N / 200);
                    x_T = x_T(1:step:end, :);
                end

                plot(x_T(:,1), x_T(:,2), 'b.', 'MarkerSize', 5);
                hold on;

                % Dibujar elipse
                if isfield(params, 'semi_major') && isfield(params, 'semi_minor') && isfield(params, 'angle_deg')
                    angle_rad = deg2rad(params.angle_deg);
                    cx = 0; cy = 0;
                    if isfield(params, 'center')
                        cx = params.center(1);
                        cy = params.center(2);
                    end
                    draw_ellipse_on_axes(cx, cy, params.semi_major, params.semi_minor, angle_rad, 'r-', 2);
                end

                xlabel('x_1 (rotado)'); ylabel('x_2 (rotado)');
            end

        case 'PE3D'
            % Para PE3D, mostrar proyección simple en plano Va-Vb
            N = size(V_window, 1);
            if N > 200
                step = ceil(N / 200);
                V_plot = V_window(1:step:end, :);
            else
                V_plot = V_window;
            end

            plot(V_plot(:,1), V_plot(:,2), 'b.', 'MarkerSize', 5);
            xlabel('V_a'); ylabel('V_b');
    end

    axis equal;
    grid on;

    if is_correct
        title_color = [0 0.5 0];
    else
        title_color = [0.8 0 0];
    end

    title(sprintf('Proyección 2D (%s): %s', method_name, classification), ...
          'FontWeight', 'bold', 'Color', title_color);
end

function plot_parameters(fig, params, method_name, classification)
    data = guidata(fig);

    axes(data.ax_params);
    cla;
    axis off;

    % Construir string de parámetros
    str = sprintf('%s - Parámetros\n', method_name);
    str = [str, '════════════════════════════\n'];
    str = [str, sprintf('Clasificación: %s\n\n', classification)];

    switch method_name
        case 'Clarke'
            if isfield(params, 'semi_major')
                str = [str, sprintf('Semieje mayor: %.4f\n', params.semi_major)];
            end
            if isfield(params, 'semi_minor')
                str = [str, sprintf('Semieje menor: %.4f\n', params.semi_minor)];
            end
            if isfield(params, 'roundness')
                str = [str, sprintf('Redondez: %.4f\n', params.roundness)];
            end
            if isfield(params, 'V0_normalized')
                str = [str, sprintf('V0 normalizada: %.4f\n', params.V0_normalized)];
            end
            if isfield(params, 'angle_deg')
                str = [str, sprintf('Ángulo: %.2f°\n', params.angle_deg)];
            end

        case 'GA'
            if isfield(params, 'bivector')
                str = [str, sprintf('Bivector:\n')];
                str = [str, sprintf('  σ12: %.4f\n', params.bivector(1))];
                str = [str, sprintf('  σ23: %.4f\n', params.bivector(2))];
                str = [str, sprintf('  σ31: %.4f\n', params.bivector(3))];
            end
            if isfield(params, 'bivector_norm')
                str = [str, sprintf('\nNorma del bivector: %.4f\n', params.bivector_norm)];
            end
            if isfield(params, 'angle_deg')
                str = [str, sprintf('Ángulo de inclinación: %.2f°\n', params.angle_deg)];
            end
            if isfield(params, 'semi_major')
                str = [str, sprintf('\nSemieje mayor: %.4f\n', params.semi_major)];
            end
            if isfield(params, 'semi_minor')
                str = [str, sprintf('Semieje menor: %.4f\n', params.semi_minor)];
            end
            if isfield(params, 'roundness')
                str = [str, sprintf('Redondez: %.4f\n', params.roundness)];
            end
            if isfield(params, 'is_line')
                str = [str, sprintf('Es línea: %s\n', mat2str(params.is_line))];
            end

        case 'PE3D'
            if isfield(params, 'Ax')
                str = [str, sprintf('Ax: %.4f\n', params.Ax)];
            end
            if isfield(params, 'Ay')
                str = [str, sprintf('Ay: %.4f\n', params.Ay)];
            end
            if isfield(params, 'Ay_Ax_ratio')
                str = [str, sprintf('Ratio Ay/Ax: %.4f\n', params.Ay_Ax_ratio)];
            end
            if isfield(params, 'theta')
                str = [str, sprintf('θ (azimuth): %.2f°\n', params.theta)];
            end
            if isfield(params, 'phi')
                str = [str, sprintf('φ (elevación): %.2f°\n', params.phi)];
            end
            if isfield(params, 'psi')
                str = [str, sprintf('ψ (orientación): %.2f°\n', params.psi)];
            end
    end

    text(0.05, 0.95, sprintf(str), ...
         'Units', 'normalized', ...
         'VerticalAlignment', 'top', ...
         'FontSize', 10, ...
         'FontName', 'Consolas', ...
         'BackgroundColor', [0.95 0.95 1], ...
         'EdgeColor', [0.5 0.5 0.8], ...
         'Margin', 8);

    title('Parámetros del Método', 'FontWeight', 'bold');
end

function plot_accuracy_evolution(fig)
    data = guidata(fig);

    if isempty(data.correct_array)
        return;
    end

    axes(data.ax_accuracy);
    cla;

    % Calcular precisión acumulada
    n = length(data.correct_array);
    cumulative_accuracy = zeros(n, 1);
    for i = 1:n
        cumulative_accuracy(i) = sum(data.correct_array(1:i)) / i * 100;
    end

    % Precisión por ventana (suavizada)
    window_smooth = min(10, n);
    smoothed_accuracy = movmean(data.correct_array * 100, window_smooth);

    t = data.time_centers * 1000;

    plot(t, cumulative_accuracy, 'b-', 'LineWidth', 2); hold on;
    plot(t, smoothed_accuracy, 'g-', 'LineWidth', 1.5);

    % Marcar posición actual
    if isfield(data, 'current_window_idx')
        idx = data.current_window_idx;
        xline(t(idx), 'r-', 'LineWidth', 2);
        plot(t(idx), cumulative_accuracy(idx), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    end

    % Marcar zona de falla
    yl = ylim;
    patch([data.fault_start, data.fault_end, data.fault_end, data.fault_start] * 1000, ...
          [yl(1), yl(1), yl(2), yl(2)], ...
          [1 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

    xlabel('Tiempo (ms)');
    ylabel('Precisión (%)');
    title('Evolución de Precisión', 'FontWeight', 'bold');
    legend('Acumulada', 'Suavizada', 'Location', 'best');
    grid on;
    ylim([0 105]);
end

%% ========================================================================
%  FUNCIONES AUXILIARES
%% ========================================================================

function is_match = check_classification_match(classified, expected)
    % Normalizar nombres
    classified = upper(strtrim(classified));
    expected = upper(strtrim(expected));

    % Remover guiones
    classified = strrep(classified, '-', '');
    expected = strrep(expected, '-', '');

    % Mapeos de equivalencias
    equivalences = {
        {'NORMAL', 'ABC', 'BALANCED'};
        {'ABCG', '3PHASE', 'THREEPHASE', 'ABCG'};
        {'AB', 'LLAB'};
        {'BC', 'LLBC'};
        {'CA', 'AC', 'LLCA', 'LLAC'};
        {'AG', 'LGA'};
        {'BG', 'LGB'};
        {'CG', 'LGC'};
        {'ABG', 'LLGAB'};
        {'BCG', 'LLGBC'};
        {'CAG', 'ACG', 'LLGCA', 'LLGAC'};
    };

    % Buscar en equivalencias
    for i = 1:length(equivalences)
        group = equivalences{i};
        if any(strcmp(classified, group)) && any(strcmp(expected, group))
            is_match = true;
            return;
        end
    end

    % Comparación directa
    is_match = strcmp(classified, expected);
end

function draw_ellipse_on_axes(cx, cy, a, b, angle, style, linewidth)
    theta = linspace(0, 2*pi, 100);
    x = a * cos(theta);
    y = b * sin(theta);

    R = [cos(angle), -sin(angle); sin(angle), cos(angle)];
    xy = R * [x; y];

    plot(xy(1,:) + cx, xy(2,:) + cy, style, 'LineWidth', linewidth);
end

function [V_abc, time_vector] = extract_converter_data(data, converter_name)
    % Función auxiliar para extraer datos (duplicada de utils para independencia)

    if ~isfield(data, 'GridConvs')
        error('El archivo .mat no contiene el campo GridConvs');
    end

    GridConvs = data.GridConvs;

    if ~isfield(GridConvs, converter_name)
        error('Convertidor %s no encontrado en GridConvs', converter_name);
    end

    converter_data = GridConvs.(converter_name);

    switch converter_name
        case {'DG1', 'DG2'}
            if ~isfield(converter_data, 'Vac')
                error('Convertidor %s no tiene campo Vac', converter_name);
            end
            V_timeseries = converter_data.Vac;
            V_abc = V_timeseries.Data;
            time_vector = V_timeseries.Time;

        case 'Pex'
            if ~isfield(converter_data, 'Vabc')
                error('Convertidor Pex no tiene campo Vabc');
            end
            V_timeseries = converter_data.Vabc;
            V_abc = V_timeseries.Data;
            time_vector = V_timeseries.Time;

        otherwise
            error('Convertidor desconocido: %s', converter_name);
    end

    if size(time_vector, 2) > size(time_vector, 1)
        time_vector = time_vector';
    end
end
