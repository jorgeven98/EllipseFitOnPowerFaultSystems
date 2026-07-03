function [fault_type, params] = ga_fault_classifier_pure(V_abc, fs, f)
% GA_FAULT_CLASSIFIER_PURE Port MATLAB-puro de ga_fault_classifier (sin GA-FuL).
% Matemáticamente equivalente: el producto exterior x1^x2 se computa como
% producto vectorial (dualidad bivector<->normal en 3D) y el sandwich del
% rotor R x R† como la rotación de Rodrigues que lleva la normal del plano
% al eje z. Paridad verificada etiqueta a etiqueta contra el original
% (ver tests/test_pure_vs_gaful_parity.m).
%
% Input/Output: idénticos a ga_fault_classifier.

%% 0. PARÁMETROS DE ENTRADA
N = size(V_abc, 1);
if nargin < 3 || isempty(f), f = 50; end
if nargin < 2 || isempty(fs)
    fs = estimate_sampling_frequency(V_abc, f);
end

%% 1-2. VECTORES 3D Y BIVECTOR DEL PLANO (via producto vectorial)
x_vectors = V_abc;
quarter_cycle_samples = round(fs / (4 * f));
if quarter_cycle_samples < 1, quarter_cycle_samples = 1; end

if N < quarter_cycle_samples
    warning('GA:ShortWindow', ...
        'Ventana de datos muy corta. Se necesitan al menos %d muestras (cuarto de ciclo).', ...
        quarter_cycle_samples + 1);
    idx1 = 1; idx2 = N; idx3 = 1;
else
    center_idx = round(N / 2);
    idx1 = max(1, center_idx - quarter_cycle_samples);
    idx2 = min(N, idx1 + quarter_cycle_samples);
    idx3 = round(abs(idx1 - idx2)/2);
end

x1 = x_vectors(idx1, :); x1 = x1 / norm(x1);
x2 = x_vectors(idx2, :); x2 = x2 / norm(x2);
x3 = x_vectors(idx3, :); x3 = x3 / norm(x3);

% Bivector B = x1 ^ x2  <=>  normal n = x1 x x2 (dualidad):
% n = [B23, B31, B12]
nvec = cross(x1, x2);
B_norm_val = norm(nvec);                 % ||B|| = ||x1 x x2||
K_ref = 1/sqrt(3);

nK = [1 1 1]/sqrt(3);
PinB = abs([dot(x1, nK), dot(x2, nK), dot(x3, nK)]);

if all(PinB < 1e-2)
    is_kirchoff = true;
    n_unit = [1 1 1]/sqrt(3);            % B_norm = K_ref*[1 -1 1] en (s12,s13,s23) => n = (1,1,1)/sqrt(3)
else
    is_kirchoff = false;
    n_unit = nvec / B_norm_val;
end

B23 = n_unit(1); B31 = n_unit(2); B12 = n_unit(3);

%% 3. CLASIFICACIÓN PRIMARIA POR COMPONENTES DE BIVECTOR
bivector_components = [B12, B23, B31];
if ~is_kirchoff
    ground_fault_type = classify_by_bivector(bivector_components);
    is_ground_fault = true;
else
    ground_fault_type = "";
    is_ground_fault = false;
end

%% 4-5. ROTAR AL PLANO sigma12: Rodrigues n_unit -> z (acción vectorial del rotor)
z = [0 0 1];
c = dot(n_unit, z);
if c > 1 - 1e-12
    R = eye(3);
elseif c < -1 + 1e-12
    R = diag([1 -1 -1]);                 % 180 grados alrededor de x
else
    a = cross(n_unit, z); a = a / norm(a);
    s = sqrt(1 - c^2);
    K = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
    R = eye(3) + s*K + (1-c)*(K*K);
end
x_rot = (R * x_vectors')';
x_T = x_rot(:, 1:2);

%% 6. AJUSTE DE ELIPSE USANDO GAC
is_line = false;
try
    if B_norm_val > 1e-3
        [center, axes, angle] = fit_ellipse_gac_improved(x_T);
        a_ = axes(1); b_ = axes(2); theta = angle;
        if sqrt(center(1)^2 + center(2)^2) > 0.1 * max(a_, b_)
            warning('Centro no está en origen: [%.4f, %.4f]', center(1), center(2));
        end
        if b_ < 0.01 && b_/a_ < 0.3
            is_line = true;
            theta = fit_line(x_T);
            b_ = 0;
        end
    else
        is_line = true;
        theta = fit_line(x_T);
        a_ = sqrt(3)/sqrt(2);
        b_ = 0;
        center = [0, 0];
    end
catch ME
    warning(ME.identifier, 'Error en ajuste GAC: %s', ME.message);
    [fault_type, params] = fallback_analysis(V_abc);
    return;
end

%% 7. NORMALIZAR PARÁMETROS DE ELIPSE
if b_ > a_
    [a_, b_] = deal(b_, a_);
    theta = theta + pi/2;
end
theta_deg = mod(rad2deg(theta), 180);
roundness = b_ / a_;

%% 8. CLASIFICACIÓN FINAL (idéntica al original)
if is_ground_fault
    fault_type = ground_fault_type;
elseif is_line || roundness < 0.05
    fault_type = classify_by_angle(theta_deg);
    if ~is_line && roundness < 0.05
        is_line = true;
    end
elseif roundness > 0.99
    V_nominal_expected = sqrt(3)/sqrt(2);
    magnitude_ratio = a_ / V_nominal_expected;
    MAGNITUDE_THRESHOLD = 0.90;
    if magnitude_ratio >= MAGNITUDE_THRESHOLD
        fault_type = 'A-B-C';
    else
        fault_type = 'A-B-C-G';
    end
else
    fault_type = classify_by_angle(theta_deg);
end

%% 9. EMPAQUETAR RESULTADOS
params.bivector = [B12, B23, B31];
params.bivector_norm = B_norm_val;
params.semi_major = a_;
params.semi_minor = b_;
params.angle_deg = theta_deg;
params.roundness = roundness;
params.is_line = is_line;
params.x_T = x_T;
params.center = center;

end

%% ========================================================================
%  FUNCIONES AUXILIARES (idénticas al original)
%% ========================================================================

function fault_type = classify_by_bivector(biv)
[biv_sort, idx] = sort(abs(biv), "descend");
dif1 = abs(biv_sort(1) - biv_sort(2));
dif2 = abs(biv_sort(2) - biv_sort(3));
if dif1 < dif2 || dif2 > 0.3
    switch idx(3)
        case 1, fault_type = 'A-B-G';
        case 2, fault_type = 'B-C-G';
        case 3, fault_type = 'C-A-G';
    end
else
    switch idx(1)
        case 1, fault_type = 'C-G';
        case 2, fault_type = 'A-G';
        case 3, fault_type = 'B-G';
    end
end
end

function fault_type = classify_by_angle(theta_deg)
MARGIN = 30;
theta_norm = mod(theta_deg, 180);
if theta_norm >= (45 - MARGIN) && theta_norm <= (45 + MARGIN)
    fault_type = 'A-B';
elseif theta_norm >= (105 - MARGIN) && theta_norm <= (105 + MARGIN)
    fault_type = 'C-A';
elseif theta_norm >= (165 - MARGIN) || theta_norm <= 15
    fault_type = 'B-C';
else
    fault_type = 'L-L';
end
end

function [fault_type, params] = fallback_analysis(V_abc)
Va_rms = rms(V_abc(:,1)); Vb_rms = rms(V_abc(:,2)); Vc_rms = rms(V_abc(:,3));
unbalance = std([Va_rms, Vb_rms, Vc_rms]) / mean([Va_rms, Vb_rms, Vc_rms]);
if unbalance < 0.05
    fault_type = 'A-B-C';
else
    fault_type = 'Unclassified';
end
params.bivector = [0, 0, 0];
params.bivector_norm = 0;
params.semi_major = 0;
params.semi_minor = 0;
params.angle_deg = 0;
params.roundness = 0;
params.is_line = false;
params.x_T = [];
params.center = [0, 0];
end

function angle = fit_line(points)
if std(points(:, 1)) < 1e-10
    angle = pi/2;
    return;
end
coef = polyfit(points(:, 1), points(:, 2), 1);
angle = mod(atan(coef(1)), pi);
if ~isfinite(angle)
    angle = pi/2;
end
end

function fs_est = estimate_sampling_frequency(V_abc, f)
Va = V_abc(:, 1);
zero_crossings = find(diff(sign(Va)) ~= 0);
if length(zero_crossings) < 2
    warning('GA:NoZeroCrossings', ...
        'No se detectaron cruces por cero. Usando frecuencia de muestreo por defecto.');
    fs_est = 1000;
    return;
end
intervals = diff(zero_crossings);
avg_half_period = mean(intervals);
fs_est = 2 * avg_half_period * f;
if fs_est < 10 * f || fs_est > 1e6
    warning('GA:UnrealisticSamplingRate', ...
        'Frecuencia de muestreo estimada (%.1f Hz) fuera de rango esperado.', fs_est);
    fs_est = max(10 * f, min(fs_est, 1e6));
end
end
