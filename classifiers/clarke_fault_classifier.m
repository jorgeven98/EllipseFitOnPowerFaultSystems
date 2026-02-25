function [fault_type, params] = clarke_fault_classifier(V_abc)
% CLARKE_ELLIPSE_FAULT_CLASSIFIER Clasificación de faltas mediante Clarke + Elipse
%
% Input:
%   V_abc: Matriz [N x 3] donde cada columna es una fase (Va, Vb, Vc)
%
% Output:
%   fault_type: String con el tipo de falta
%   params: Estructura con parámetros de diagnóstico

%% 1. TRANSFORMADA DE CLARKE
K = sqrt(2/3);
T_clarke = K * [1,        -1/2,           -1/2;
                0,  sqrt(3)/2,  -sqrt(3)/2;
                1/sqrt(2), 1/sqrt(2), 1/sqrt(2)];

N = size(V_abc, 1);
V_alpha = zeros(N, 1);
V_beta = zeros(N, 1);
V_0 = zeros(N, 1);

for i = 1:N
    V_ab0 = T_clarke * V_abc(i, :)';
    V_alpha(i) = V_ab0(1);
    V_beta(i) = V_ab0(2);
    V_0(i) = V_ab0(3);
end

%% 2. AJUSTE DE ELIPSE CON CENTRO EN ORIGEN
ellipse_params = fit_ellipse_centered_origin(V_alpha, V_beta);

if isempty(ellipse_params)
    warning('Ellipse fitting failed - using fallback analysis');
    [fault_type, params] = fallback_analysis(V_abc, V_alpha, V_beta, V_0);
    return;
end

a = ellipse_params.semi_major;
b = ellipse_params.semi_minor;
theta = ellipse_params.angle_rad;

%% 3. MÉTRICAS GEOMÉTRICAS
roundness = b / a;
eccentricity = sqrt(1 - (b/a)^2);

V0_rms = rms(V_0);
V_total_rms = rms(sqrt(V_alpha.^2 + V_beta.^2));
V0_normalized = V0_rms / (V_total_rms + 1e-10);

theta_deg = mod(rad2deg(theta), 180);

%% 4. ANÁLISIS DE FASES AFECTADAS Y MAGNITUD
% Calcular RMS de cada fase
Va_rms = rms(V_abc(:,1));
Vb_rms = rms(V_abc(:,2));
Vc_rms = rms(V_abc(:,3));
V_rms_array = [Va_rms, Vb_rms, Vc_rms];
V_mean = mean(V_rms_array);

% Magnitud normalizada (asumiendo nominal = 1.0 pu RMS)
V_nominal = 1.0;  % Valor nominal esperado en pu
magnitude_ratio = V_mean / V_nominal;

% Identificar fases con voltaje reducido (< 80% de la media)
threshold_fault = 0.80;
fault_phases = V_rms_array < (threshold_fault * V_mean);
num_fault_phases = sum(fault_phases);

%% 5. CLASIFICACIÓN BASADA EN NÚMERO DE FASES AFECTADAS
% MODIFICACIÓN: Umbrales ajustados según criterios de GA para coherencia entre métodos
TH_BALANCED = 0.99;    % Roundness > 0.99 → Sistema circular (Normal o ABCG)
TH_V0_GROUND = 0.15;   % V0 normalizada > 0.15 → Componente de tierra presente
TH_DEGENERATE = 0.05;  % MODIFICADO: 0.15 → 0.05 (coherente con GA, umbral para líneas)
TH_ELLIPTICAL = 0.60;  % Roundness intermedio para faltas L-L
TH_MAGNITUDE = 0.90;   % MODIFICADO: 0.85 → 0.90 (coherente con GA, 90% del nominal)

phase_names = {'A', 'B', 'C'};

% CASO ESPECIAL: Sistema balanceado (circulo) pero magnitud reducida
% MODIFICACIÓN: Criterio mejorado para distinguir Normal vs ABCG usando magnitud
% Para distinguir Normal vs ABCG, usar la magnitud del semieje mayor
% Valor nominal esperado en Clarke: ~sqrt(2) * V_rms para sistema equilibrado
V_clarke_nominal = sqrt(2) * V_nominal;  % Amplitud pico esperada
magnitude_ratio_clarke = a / V_clarke_nominal;

if roundness > TH_BALANCED && V0_normalized < 0.10
    % MODIFICACIÓN: Umbral de magnitud ajustado a 0.90 (coherente con GA)
    if magnitude_ratio < TH_MAGNITUDE || magnitude_ratio_clarke < TH_MAGNITUDE
        % Círculo pequeño (<90% nominal) → Falta trifásica simétrica ABCG
        fault_type = 'A-B-C-G';
    else
        % Círculo normal (≥90% nominal) → Sistema balanceado sin falta
        fault_type = 'A-B-C';
    end
    
elseif V0_normalized > TH_V0_GROUND
    % Falta con componente de tierra
    
    if num_fault_phases == 1
        % UNA sola fase afectada → L-G
        faulted_idx = find(fault_phases);
        fault_type = [phase_names{faulted_idx} '-G'];

    elseif num_fault_phases == 2
        % DOS fases afectadas → L-L-G
        faulted_idx = find(fault_phases);
        llg_map = containers.Map({'AB','BC','AC'}, {'A-B-G','B-C-G','C-A-G'});
        pair_key = [phase_names{faulted_idx(1)} phase_names{faulted_idx(2)}];
        fault_type = llg_map(pair_key);

    elseif num_fault_phases == 3
        % Todas las fases afectadas con V0 → Probablemente A-B-C-G o L-L-G severa
        if roundness > 0.85
            fault_type = 'A-B-C-G';
        else
            fault_type = 'L-L-G';
        end

    else
        % No se detectan fases afectadas claramente
        if roundness > TH_ELLIPTICAL
            [~, min_phase] = min(V_rms_array);
            fault_type = [phase_names{min_phase} '-G'];
        else
            llg_angle_map = containers.Map({'AB','BC','CA'}, {'A-B-G','B-C-G','C-A-G'});
            fault_type = llg_angle_map(identify_phases_by_angle(theta_deg));
        end
    end
    
else
    % Faltas sin componente de tierra (V0 despreciable)
    
    if roundness < TH_DEGENERATE
        % Elipse muy elongada o línea
        if num_fault_phases >= 2
            if num_fault_phases == 3 || roundness < 0.05
                fault_type = 'A-B-C-G';
            else
                fault_type = classify_LL_by_angle(theta_deg);
            end
        else
            fault_type = classify_LL_by_angle(theta_deg);
        end

    elseif roundness < TH_ELLIPTICAL
        % Elipse moderadamente elongada → L-L
        fault_type = classify_LL_by_angle(theta_deg);

    else
        % Roundness alto, sin V0
        if magnitude_ratio < TH_MAGNITUDE
            % Magnitud reducida → A-B-C-G
            fault_type = 'A-B-C-G';
        elseif roundness > 0.85
            fault_type = 'A-B-C';
        else
            fault_type = 'Unclassified';
        end
    end
end

%% 6. RESULTADOS
params.semi_major = a;
params.semi_minor = b;
params.angle_deg = theta_deg;
params.angle_rad = theta;
params.roundness = roundness;
params.V0_normalized = V0_normalized;
params.eccentricity = eccentricity;
params.center = [0, 0];
params.V_alpha = V_alpha;
params.V_beta = V_beta;
params.num_fault_phases = num_fault_phases;
params.phase_rms = V_rms_array;
params.magnitude_ratio = magnitude_ratio;

end

%% ========================================================================
%  FUNCIONES AUXILIARES
%% ========================================================================

function ellipse = fit_ellipse_centered_origin(x, y)
% Ajuste de elipse centrada en el origen usando PCA

n = length(x);
if n < 6
    ellipse = [];
    return;
end

% Remover valores inválidos
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);

if length(x) < 6
    ellipse = [];
    return;
end

% Los datos ya están centrados en origen
data_matrix = [x, y];

% Matriz de covarianza
cov_matrix = (data_matrix' * data_matrix) / (n - 1);

% Eigenvalores y eigenvectores (PCA)
[eigenvectors, eigenvalues] = eig(cov_matrix);
eigenvalues = diag(eigenvalues);

% Ordenar por eigenvalor (mayor primero)
[eigenvalues, idx] = sort(eigenvalues, 'descend');
eigenvectors = eigenvectors(:, idx);

% Proyectar datos sobre los ejes principales
axis1 = eigenvectors(:, 1);
axis2 = eigenvectors(:, 2);

proj1 = data_matrix * axis1;
proj2 = data_matrix * axis2;

% Los semiejes son la mitad del rango total
semi_a = (max(proj1) - min(proj1)) / 2;
semi_b = (max(proj2) - min(proj2)) / 2;

% Verificar que semi_a >= semi_b
if semi_b > semi_a
    temp = semi_a;
    semi_a = semi_b;
    semi_b = temp;
    temp_axis = axis1;
    axis1 = axis2;
    axis2 = temp_axis;
end

% El ángulo es dado por el eje mayor
theta = atan2(axis1(2), axis1(1));

% Normalizar ángulo a [0, π]
theta = mod(theta, pi);

% Para círculos casi perfectos, forzar ángulo a 0
if abs(semi_a - semi_b) / semi_a < 0.05
    theta = 0;
end

% Verificar validez
if ~isfinite(semi_a) || ~isfinite(semi_b) || semi_a <= 0 || semi_b <= 0
    ellipse = [];
    return;
end

ellipse.semi_major = semi_a;
ellipse.semi_minor = semi_b;
ellipse.angle_rad = theta;
ellipse.center_x = 0;
ellipse.center_y = 0;

end

function phases_str = identify_phases_by_angle(theta_deg)
% Identificar fases basándose en el ángulo de la elipse

theta_norm = mod(theta_deg, 180);

if theta_norm >= 15 && theta_norm <= 75
    phases_str = 'AB';
elseif theta_norm >= 75 && theta_norm <= 135
    phases_str = 'CA';
else
    phases_str = 'BC';
end

end

function [fault_type, params] = fallback_analysis(V_abc, V_alpha, V_beta, V_0)
% Análisis de respaldo cuando falla el ajuste de elipse

V0_rms = rms(V_0);
Valpha_rms = rms(V_alpha);
Vbeta_rms = rms(V_beta);
V_total_rms = sqrt(Valpha_rms^2 + Vbeta_rms^2);

Va_rms = rms(V_abc(:,1));
Vb_rms = rms(V_abc(:,2));
Vc_rms = rms(V_abc(:,3));

unbalance = std([Va_rms, Vb_rms, Vc_rms]) / mean([Va_rms, Vb_rms, Vc_rms]);
V0_normalized = V0_rms / (V_total_rms + 1e-10);

if unbalance < 0.05
    fault_type = 'A-B-C';
elseif V0_normalized > 0.15
    [~, min_phase] = min([Va_rms, Vb_rms, Vc_rms]);
    phase_names = {'A', 'B', 'C'};
    fault_type = [phase_names{min_phase} '-G'];
else
    fault_type = 'Unclassified';
end

params.semi_major = (max(V_alpha) - min(V_alpha)) / 2;
params.semi_minor = (max(V_beta) - min(V_beta)) / 2;
params.angle_deg = 0;
params.roundness = params.semi_minor / (params.semi_major + 1e-10);
params.V0_normalized = V0_normalized;
params.eccentricity = sqrt(1 - params.roundness^2);
params.center = [0, 0];
params.V_alpha = V_alpha;
params.V_beta = V_beta;
params.num_fault_phases = 0;
params.phase_rms = [Va_rms, Vb_rms, Vc_rms];
params.magnitude_ratio = mean([Va_rms, Vb_rms, Vc_rms]);

end

function fault_type = classify_LL_by_angle(theta_deg)
% CLASSIFY_LL_BY_ANGLE Clasificación L-L por ángulo de orientación de la elipse
% MODIFICACIÓN: Margen ampliado a ±30° (coherente con GA) para mayor robustez
% Ángulos teóricos: 45° (A-B), 105° (C-A), 165° (B-C)

theta_norm = mod(theta_deg, 180);
MARGIN = 30;  % MODIFICADO: Margen de ±30° (coherente con GA)

% Rangos de clasificación (coherentes con GA):
% A-B:  15° - 75°  (centrado en 45°)
% C-A:  75° - 135° (centrado en 105°)
% B-C:  135° - 180° y 0° - 15° (centrado en 165°/0°)

if theta_norm >= (45 - MARGIN) && theta_norm <= (45 + MARGIN)
    fault_type = 'A-B';  % Falta entre fases A y B
elseif theta_norm >= (105 - MARGIN) && theta_norm <= (105 + MARGIN)
    fault_type = 'C-A';  % Falta entre fases C y A
elseif theta_norm >= (165 - MARGIN) || theta_norm <= 15
    fault_type = 'B-C';  % Falta entre fases B y C (envuelve 180°/0°)
else
    fault_type = 'L-L';  % Falta L-L genérica (ángulo ambiguo)
end

end