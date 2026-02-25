function [fault_type, params] = ga_fault_classifier(V_abc, fs, f)
% GA_ELLIPSE_FAULT_CLASSIFIER Clasificación según metodología del paper
% Usa rotación GA (Montoya Transform) + ajuste GAC + clasificación por bivectores
% Requiere: GA-FuL MATLAB Toolbox
%
% Input:
%   V_abc: Matriz [N x 3] donde cada columna es una fase (Va, Vb, Vc)
%   fs: Frecuencia de muestreo (Hz) [opcional, default: detectar automáticamente]
%   f: Frecuencia del sistema (Hz) [opcional, default: 60 Hz]
%
% Output:
%   fault_type: String con el tipo de falta
%   params: Estructura con parámetros de diagnóstico

persistent vga;
if isempty(vga)
    gafulInit % VGA processor
end

%% 0. PARÁMETROS DE ENTRADA
N = size(V_abc, 1);

% Frecuencia del sistema por defecto
if nargin < 3 || isempty(f)
    f = 50; % Hz
end

% Detectar frecuencia de muestreo si no se proporciona
if nargin < 2 || isempty(fs)
    % Estimar asumiendo que se tienen datos de al menos 1 ciclo
    fs = estimate_sampling_frequency(V_abc, f);
end

%% 1. CONSTRUIR VECTORES 3D EN BASE ORTOGONAL σ
% Vectores en R^3 según ec. (1) del paper
x_vectors = V_abc;

%% 2. CALCULAR BIVECTOR DEL PLANO (ec. 2 del paper)
% B = x_t1 ∧ x_t2 (outer product)
% IMPORTANTE: Los vectores deben estar separados por 90° (cuarto de ciclo)

% Calcular índice de cuarto de ciclo
quarter_cycle_samples = round(fs / (4 * f));

% Validar que tenemos suficientes muestras
if quarter_cycle_samples < 1
    quarter_cycle_samples = 1;
    % warning('GA:InsufficientSamples', ...
    %     'Muy pocas muestras por ciclo (%d). Resultados pueden ser imprecisos.', ...
    %     round(fs/f));
end

if N < quarter_cycle_samples 
    warning('GA:ShortWindow', ...
        'Ventana de datos muy corta. Se necesitan al menos %d muestras (cuarto de ciclo).', ...
        quarter_cycle_samples + 1);
    % Usar primer y último punto disponibles
    idx1 = 1;
    idx2 = N;
else
    % Usar puntos separados por cuarto de ciclo desde el centro de la ventana
    center_idx = round(N / 2);
    idx1 = max(1, center_idx - quarter_cycle_samples);
    idx2 = min(N, idx1 + quarter_cycle_samples);
    idx3 = round(abs(idx1 - idx2)/2);
end

% Construir bivector con vectores separados correctamente
x1_vec = vga.Vector(x_vectors(idx1, :)).DivideByNorm;
x2_vec = vga.Vector(x_vectors(idx2, :)).DivideByNorm;
x3_vec = vga.Vector(x_vectors(idx3, :)).DivideByNorm;

% Bivector mediante outer product
B = x1_vec.Op(x2_vec);

% Normalizar bivector
B_norm = B.DivideByNorm;
K_ref = 1/sqrt(3);

n = vga.Vector([1 1 1]/sqrt(3));
PinB = abs([x1_vec.Sp(n), x2_vec.Sp(n), x3_vec.Sp(n)]);

if all(PinB < 1e-2)
    is_kirchoff = true;
    B_norm = vga.Bivector([1 -1 1]*K_ref);
else
    is_kirchoff = false;
end


B_array = double(B_norm.BivectorToArray1D(3));

% Mapeo de componentes del bivector en 3D:
% Posición 1: <0,1> = σ1∧σ2 = σ12
% Posición 2: <0,2> = σ1∧σ3 = σ13 = -σ31
% Posición 3: <1,2> = σ2∧σ3 = σ23

B12 = B_array(1);   % σ12
B13 = B_array(2);   % σ13
B23 = B_array(3);   % σ23

% Convertir a convención del paper (σ12, σ23, σ31)
% σ31 = -σ13
B31 = -B13;

%% 3. CLASIFICACIÓN PRIMARIA POR COMPONENTES DE BIVECTOR
% MODIFICACIÓN: Detección mejorada de faltas con tierra mediante análisis de Kirchhoff
% Si el producto interno de los vectores con el normal al plano (1,1,1) es cercano a cero,
% el sistema cumple la ley de Kirchhoff (Va + Vb + Vc ≈ 0), indicando ausencia de componente homopolar
tolerance = 0.02;
bivector_components = [B12, B23, B31];

if ~is_kirchoff
    % MODIFICACIÓN: Sistema NO cumple Kirchhoff → Existe componente de tierra (faltas L-G o L-L-G)
    ground_fault_type = classify_by_bivector(bivector_components);
    is_ground_fault = true;
else
    % Sistema cumple Kirchhoff → No hay componente de tierra (faltas L-L o balanceadas)
    ground_fault_type = "";
    is_ground_fault = false;
end
%% 4. CALCULAR ROTOR PARA ALINEAR PLANO B CON σ12
% Rotor R que alinea bivector B con bivector σ12 (plano xy)
% Según ec. (24) del paper: R = (1 + Â B̂†) / ||1 + Â B̂†||

% Bivector objetivo: σ12 (solo componente <0,1>)
A = vga.Bivector([1, 0, 0]);% [σ12, σ13, σ23]

% Reverso de B_hat: B̂†
B_rev = B_norm.Reverse();

% Producto geométrico: 1 + Â B̂†
unity = vga.Scalar(1);
R_unnorm = unity.Add(A.Gp(B_rev));

% Normalizar rotor
R = R_unnorm.DivideByNorm();


%% 5. ROTAR VECTORES AL PLANO σ12 (ec. 3 y 10 del paper)
% x_T = R x R†
% OPTIMIZACIÓN: Pre-calcular el reverso del rotor
x_T = zeros(N, 2);
R_rev = R.Reverse();

% Pre-alocar para vectores GA (optimización de memoria)
for i = 1:N
    % Crear vector en GA y aplicar sandwich product: R x R†
    x = vga.Vector(x_vectors(i, :));
    x_rotated = R.Gp(x).Gp(R_rev);
    

    % Extraer componentes σ1 y σ2 directamente
    %x_array = double(x_rotated.GetVectorPart().VectorToArray1D());
    x_T(i, 1) = x_rotated.VectorPartToVector2D.Scalar1.ScalarValue; %x_array(1);  % σ1
    x_T(i, 2) = x_rotated.VectorPartToVector2D.Scalar2.ScalarValue; %x_array(2);  % σ2
end
%% 6. AJUSTE DE ELIPSE USANDO GAC
is_line = false;  % Flag para indicar si la curva es una línea

try
    if B.Norm > 1e-3
        [center, axes, angle] = fit_ellipse_gac_improved(x_T);

        a = axes(1);  % Semieje mayor
        b = axes(2);  % Semieje menor
        theta = angle;  % Ángulo


        % Verificar centro cerca del origen
        if sqrt(center(1)^2 + center(2)^2) > 0.1 * max(a, b)
            warning('Centro no está en origen: [%.4f, %.4f]', center(1), center(2));
        end

        % Detectar si la elipse es realmente una línea (b muy pequeño)
        if b < 0.01 && b/a < 0.3  % Umbral: b < 1% de a
            is_line = true;
            theta = fit_line(x_T);  % Re-calcular ángulo con ajuste lineal
            b = 0;  % Forzar semieje menor a cero
        end
    else
        % Bivector muy pequeño: curva degenerada a línea
        is_line = true;
        theta = fit_line(x_T);
        a = sqrt(3)/sqrt(2);
        b = 0;
        center = [0, 0];
    end

catch ME
    warning(ME.identifier,'Error en ajuste GAC: %s', ME.message);
    [fault_type, params] = fallback_analysis(V_abc);
    return;
end

%% 7. NORMALIZAR PARÁMETROS DE ELIPSE
% Asegurar que a >= b (semieje mayor >= semieje menor)
if b > a
    [a, b] = deal(b, a);  % Swap optimizado
    theta = theta + pi/2;
end

% Normalizar ángulo a rango [0, 180)
theta_deg = mod(rad2deg(theta), 180);

% Calcular redondez (roundness)
roundness = b / a;

%% 8. CLASIFICACIÓN FINAL
if is_ground_fault
    % Falta con tierra detectada por bivector (Kirchhoff no cumplido)
    fault_type = ground_fault_type;

elseif is_line || roundness < 0.05
    % MODIFICACIÓN: Umbral de roundness ajustado a 0.05 para mejor detección de líneas
    % Curva degenerada a línea: clasificar como falta L-L por ángulo
    % Esto ocurre cuando dos fases están en cortocircuito (elipse colapsa a línea)
    fault_type = classify_by_angle(theta_deg);
    if ~is_line && roundness < 0.05
        % Marcar que fue detectado por bajo roundness
        is_line = true;
    end

elseif roundness > 0.99
    % MODIFICACIÓN: Criterio mejorado para distinguir Normal vs ABCG
    % Sistema circular (roundness > 0.99): usar MAGNITUD del semieje mayor para discriminar
    % Valor nominal esperado: sqrt(3)/sqrt(2) ~ 1.2247 (para 1 pu)
    V_nominal_expected = sqrt(3)/sqrt(2);
    magnitude_ratio = a / V_nominal_expected;

    % MODIFICACIÓN: Umbral de magnitud optimizado a 0.90 (90% del nominal)
    % Basado en análisis empírico para mejor separación entre condiciones
    MAGNITUDE_THRESHOLD = 0.90;

    if magnitude_ratio >= MAGNITUDE_THRESHOLD
        % Magnitud normal → sistema balanceado sin falla (operación normal)
        fault_type = 'A-B-C';
    else
        % MODIFICACIÓN: Magnitud reducida (<90%) → falta trifásica simétrica ABCG
        % El círculo pequeño indica colapso simétrico de voltaje
        fault_type = 'A-B-C-G';
    end

else
    % Faltas L-L: elipse con roundness intermedio (0.05 < r < 0.99)
    % Clasificar por ángulo de orientación de la elipse
    fault_type = classify_by_angle(theta_deg);
end

%% 9. EMPAQUETAR RESULTADOS
params.bivector = [B12, B23, B31];
params.bivector_norm = B.Norm;  % Norma del bivector
params.semi_major = a;
params.semi_minor = b;
params.angle_deg = theta_deg;
params.roundness = roundness;
params.is_line = is_line;  % Indica si la curva es una línea recta
params.x_T = x_T;
params.center = center;

end


%% ========================================================================
%  FUNCIONES AUXILIARES
%% ========================================================================

function fault_type = classify_by_bivector(biv)
% CLASSIFY_BY_BIVECTOR Clasificar falta por componentes de bivector según Fig. 4 del paper
% MODIFICACIÓN: Lógica mejorada para discriminar L-G vs L-L-G basada en diferencias entre componentes
%
% Input:
%   biv: Vector [B12, B23, B31] - componentes del bivector normalizado
%
% Output:
%   fault_type: String con el tipo de falta ('A-G', 'B-G', 'C-G', 'A-B-G', 'B-C-G', 'C-A-G')
%
% Lógica de clasificación:
%   - Ordena componentes por magnitud (mayor a menor)
%   - Calcula diferencias entre componentes consecutivos
%   - Si dif1 < dif2 O dif2 > 0.3 → Falta L-L-G (dos componentes dominantes)
%   - Si no → Falta L-G (una componente dominante)

[biv_sort, idx] = sort(abs(biv), "descend");

dif1 = abs(biv_sort(1) - biv_sort(2));  % Diferencia entre 1er y 2do componente
dif2 = abs(biv_sort(2) - biv_sort(3));  % Diferencia entre 2do y 3er componente

% MODIFICACIÓN: Criterio mejorado para L-L-G
% Si dif1 < dif2: las dos primeras componentes son similares → L-L-G
% Si dif2 > 0.3: separación grande entre 2do y 3er componente → L-L-G
if dif1 < dif2 || dif2 > 0.3
    % Falta L-L-G: Dos líneas a tierra
    % La componente MÍNIMA (idx(3)) indica qué par de fases NO está en falta
    switch idx(3)
        case 1  % σ12 mínimo → Fases B y C no afectadas → Falta A-B-G
            fault_type = 'A-B-G';
        case 2  % σ23 mínimo → Fases C y A no afectadas → Falta B-C-G
            fault_type = 'B-C-G';
        case 3  % σ31 mínimo → Fases A y B no afectadas → Falta C-A-G
            fault_type = 'C-A-G';
    end
else
    % Falta L-G: Una línea a tierra
    % La componente MÁXIMA (idx(1)) indica qué fase NO está en falta
    switch idx(1)
        case 1  % σ12 máximo → Fase C en falta → C-G
            fault_type = 'C-G';
        case 2  % σ23 máximo → Fase A en falta → A-G
            fault_type = 'A-G';
        case 3  % σ31 máximo → Fase B en falta → B-G
            fault_type = 'B-G';
    end

end

end

function fault_type = classify_by_angle(theta_deg)
% CLASSIFY_BY_ANGLE Clasificar falta L-L por ángulo de orientación de la elipse
% Según Fig. 6 del paper - Ángulos teóricos: 45° (A-B), 105° (C-A), 165° (B-C)
%
% MODIFICACIÓN: Margen ampliado a ±30° para mayor robustez ante ruido
% Los rangos se solapan ligeramente para cubrir casos límite

MARGIN = 30;  % MODIFICACIÓN: Margen de ±30° (antes típicamente ±15°)
theta_norm = mod(theta_deg, 180);

% Rangos de clasificación con margen ampliado:
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

function [fault_type, params] = fallback_analysis(V_abc)
% Análisis simple de respaldo

Va_rms = rms(V_abc(:,1));
Vb_rms = rms(V_abc(:,2));
Vc_rms = rms(V_abc(:,3));

unbalance = std([Va_rms, Vb_rms, Vc_rms]) / mean([Va_rms, Vb_rms, Vc_rms]);

if unbalance < 0.05
    fault_type = 'A-B-C';
else
    fault_type = 'Unclassified';
end

params.bivector = [0, 0, 0];
params.bivector_norm = 0;  % Norma del bivector
params.semi_major = 0;
params.semi_minor = 0;
params.angle_deg = 0;
params.roundness = 0;
params.is_line = false;
params.x_T = [];
params.center = [0, 0];

end

function angle = fit_line(points)
    % FIT_LINE Ajuste de línea usando mínimos cuadrados
    %
    % Input:
    %   points - Matriz Nx2 donde cada fila es [x, y]
    %
    % Output:
    %   angle - Ángulo de la línea en radianes [0, pi)

    % Si todos los puntos tienen el mismo x, la línea es vertical → 90°
    if std(points(:, 1)) < 1e-10
        angle = pi/2;
        return;
    end

    % Ajuste lineal: y = mx + b
    coef = polyfit(points(:, 1), points(:, 2), 1);

    % Calcular ángulo y normalizarlo a [0, pi)
    angle = mod(atan(coef(1)), pi);

    % Protección por si polyfit devuelve NaN/Inf (datos casi degenerados)
    if ~isfinite(angle)
        angle = pi/2;   % Fallback: línea vertical
    end
end

function fs_est = estimate_sampling_frequency(V_abc, f)
    % ESTIMATE_SAMPLING_FREQUENCY Estima la frecuencia de muestreo
    % a partir de los datos de voltaje usando detección de cruces por cero
    %
    % Input:
    %   V_abc: Matriz [N x 3] de voltajes
    %   f: Frecuencia nominal del sistema (Hz)
    %
    % Output:
    %   fs_est: Frecuencia de muestreo estimada (Hz)

    N = size(V_abc, 1);

    % Usar la fase A para detección de cruces por cero
    Va = V_abc(:, 1);

    % Detectar cruces por cero (cambios de signo)
    zero_crossings = find(diff(sign(Va)) ~= 0);

    if length(zero_crossings) < 2
        % No hay suficientes cruces por cero, usar estimación por defecto
        warning('GA:NoZeroCrossings', ...
            'No se detectaron cruces por cero. Usando frecuencia de muestreo por defecto.');
        fs_est = 1000; % Default: 1 kHz
        return;
    end

    % Calcular período promedio entre cruces consecutivos
    % Un ciclo completo tiene 2 cruces por cero
    intervals = diff(zero_crossings);
    avg_half_period = mean(intervals);

    % Frecuencia de muestreo estimada
    % fs = muestras_por_ciclo * f = (2 * muestras_por_medio_ciclo) * f
    fs_est = 2 * avg_half_period * f;

    % Validación: frecuencia de muestreo razonable
    if fs_est < 10 * f || fs_est > 1e6
        warning('GA:UnrealisticSamplingRate', ...
            'Frecuencia de muestreo estimada (%.1f Hz) fuera de rango esperado.', fs_est);
        fs_est = max(10 * f, min(fs_est, 1e6));
    end
end

