function [fault_type, params] = r2_ga_classifier_param(V_abc, fs, f, P)
% R2_GA_CLASSIFIER_PARAM  Copia de ga_fault_classifier_pure con los umbrales
% internos expuestos como parametros (segunda revision TPWRD, tarea T03).
%
% Con P vacio o por defecto devuelve EXACTAMENTE lo mismo que
% src/fault_classifiers/ga_fault_classifier_pure.m (misma secuencia de
% operaciones en coma flotante; la identidad se comprueba ventana a ventana
% en experiments/r2_t03_paridad.m). Sirve para barrer umbrales (puerta de
% Kirchhoff, redondez, rama degenerada) sin tocar el clasificador publicado.
%
% Uso:
%   [etq, p] = r2_ga_classifier_param(Vw, fs, f)          % = pure
%   [etq, p] = r2_ga_classifier_param(Vw, fs, f, P)       % umbrales de P
%   P = r2_ga_classifier_param('defaults')                % struct por defecto
%
% Campos de P (valor por defecto = valor codificado en pure):
%   kirchhoff_thr  0.01   puerta de secuencia cero: |x^.k^| < thr en las tres
%                         muestras idx1, idx2, idx3 (linea 48 de pure)
%   bnorm_min      1e-3   ||x1 x x2|| minimo; por debajo, rama degenerada (l. 87)
%   line_b_abs     0.01   semieje menor ABSOLUTO (unidades de la entrada) por
%                         debajo del cual la elipse se trata como linea (l. 93)
%   line_ratio     0.3    y ademas b/a < line_ratio (l. 93)
%   round_line     0.05   redondez por debajo de la cual se marca linea (l. 122)
%   round_circle   0.99   redondez por encima de la cual es circulo (l. 127)
%   mag_thr        0.90   a/(sqrt(3)/sqrt(2)) >= mag_thr -> 'A-B-C' (l. 130)
%   v_nominal      sqrt(3)/sqrt(2)  semieje nominal en las unidades de entrada
%   sector_margin  30     semiancho de los sectores L-L en grados (l. 177)
%   biv_dif_thr    0.3    umbral dif2 de classify_by_bivector (l. 161)
% Variantes de la puerta de secuencia cero (T52; por defecto, la publicada):
%   v0_mode        'samples3'  estadistico que se compara con kirchhoff_thr:
%                         'samples3' = max |x^.k^| en tres muestras normalizadas
%                         (pure); 'rms_win', 'rms_nom', 'rms_max', 'ph_win',
%                         'ph_nom', 'tilt', 'tilt2' = campo homonimo de
%                         r2_v0gate_stats (p. ej. 'rms_win' = sqrt(3) rms(v0) /
%                         rms(||x||) de toda la ventana)
%   v0_nominal     NaN    pico de fase nominal, solo para 'rms_nom' y 'ph_nom'
%   gnd_mode       '2s'   etiqueta de la rama de tierra: '2s' = normal de dos
%                         muestras (pure); 'svd' = normal por SVD de toda la
%                         ventana; 'svdline' = como 'svd' y, si la ventana es una
%                         recta (sv2/sv1 < gnd_line_ratio), por su direccion principal
%   gnd_line_ratio 0.05   ver gnd_mode
%
% Salida params: los campos de pure y, ademas,
%   pin            [1x3] |x^.k^| en idx1, idx2, idx3
%   pinmax         max(pin)
%   v0stat         estadistico comparado con kirchhoff_thr (= pinmax si v0_mode = 'samples3')
%   is_kirchhoff   true si las tres pasan la puerta (rama sin tierra)
%   deg_bnorm      true si ||x1 x x2|| <= bnorm_min (rama degenerada forzada)
%   deg_bsmall     true si la linea se detecto por b < line_b_abs y b/a < line_ratio
%   fallback       true si salto fallback_analysis
%
% OJO con las unidades: line_b_abs, mag_thr y v_nominal son absolutos. Los
% datos de simulacion del paper se clasificaron en voltios (pico ~223 V) y
% los de laboratorio en pu; ver ronda2/informes/T03.md.

if ischar(V_abc) || isstring(V_abc)
    fault_type = defaults_(); params = []; return
end
if nargin < 4 || isempty(P), P = defaults_(); else, P = fill_(P); end

%% 0. PARAMETROS DE ENTRADA
N = size(V_abc, 1);
if nargin < 3 || isempty(f), f = 50; end
if nargin < 2 || isempty(fs)
    error('r2_ga_classifier_param:fs', 'Indica fs explicitamente.');
end

%% 1-2. VECTORES 3D Y BIVECTOR DEL PLANO
x_vectors = V_abc;
quarter_cycle_samples = round(fs / (4 * f));
if quarter_cycle_samples < 1, quarter_cycle_samples = 1; end

if N < quarter_cycle_samples
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

nvec = cross(x1, x2);
B_norm_val = norm(nvec);

nK = [1 1 1]/sqrt(3);
PinB = abs([dot(x1, nK), dot(x2, nK), dot(x3, nK)]);

if strcmp(P.v0_mode, 'samples3')
    gate_pass = all(PinB < P.kirchhoff_thr); v0stat = max(PinB); st0 = [];
else
    st0 = r2_v0gate_stats(V_abc, fs, f, P.v0_nominal);
    v0stat = st0.(P.v0_mode); gate_pass = v0stat < P.kirchhoff_thr;
end
if gate_pass
    is_kirchoff = true;
    n_unit = [1 1 1]/sqrt(3);
else
    is_kirchoff = false;
    n_unit = nvec / B_norm_val;
end

B23 = n_unit(1); B31 = n_unit(2); B12 = n_unit(3);

%% 3. CLASIFICACION PRIMARIA POR COMPONENTES DE BIVECTOR
bivector_components = [B12, B23, B31];
if ~is_kirchoff
    ground_fault_type = classify_by_bivector(bivector_components, P.biv_dif_thr);
    is_ground_fault = true;
else
    ground_fault_type = "";
    is_ground_fault = false;
end

%% 4-5. ROTAR AL PLANO sigma12
z = [0 0 1];
c = dot(n_unit, z);
if c > 1 - 1e-12
    R = eye(3);
elseif c < -1 + 1e-12
    R = diag([1 -1 -1]);
else
    a = cross(n_unit, z); a = a / norm(a);
    s = sqrt(1 - c^2);
    K = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
    R = eye(3) + s*K + (1-c)*(K*K);
end
x_rot = (R * x_vectors')';
x_T = x_rot(:, 1:2);

%% 6. AJUSTE DE ELIPSE USANDO GAC
is_line = false; deg_bnorm = false; deg_bsmall = false;
try
    if B_norm_val > P.bnorm_min
        [center, axes, angle] = fit_ellipse_gac_improved(x_T);
        a_ = axes(1); b_ = axes(2); theta = angle;
        if b_ < P.line_b_abs && b_/a_ < P.line_ratio
            is_line = true; deg_bsmall = true;
            theta = fit_line(x_T);
            b_ = 0;
        end
    else
        is_line = true; deg_bnorm = true;
        theta = fit_line(x_T);
        a_ = P.v_nominal;
        b_ = 0;
        center = [0, 0];
    end
catch
    [fault_type, params] = fallback_analysis(V_abc);
    params.pin = PinB; params.pinmax = max(PinB); params.is_kirchhoff = is_kirchoff;
    params.deg_bnorm = false; params.deg_bsmall = false; params.fallback = true; params.v0stat = v0stat;
    return;
end

%% 7. NORMALIZAR PARAMETROS DE ELIPSE
if b_ > a_
    [a_, b_] = deal(b_, a_);
    theta = theta + pi/2;
end
theta_deg = mod(rad2deg(theta), 180);
roundness = b_ / a_;

%% 8. CLASIFICACION FINAL
if is_ground_fault
    fault_type = ground_fault_type;
elseif is_line || roundness < P.round_line
    fault_type = classify_by_angle(theta_deg, P.sector_margin);
    if ~is_line && roundness < P.round_line
        is_line = true;
    end
elseif roundness > P.round_circle
    magnitude_ratio = a_ / P.v_nominal;
    if magnitude_ratio >= P.mag_thr
        fault_type = 'A-B-C';
    else
        fault_type = 'A-B-C-G';
    end
else
    fault_type = classify_by_angle(theta_deg, P.sector_margin);
end

% Variante de la etiqueta de tierra (T52): solo cambia la etiqueta
if is_ground_fault && ~strcmp(P.gnd_mode, '2s')
    if isempty(st0), st0 = r2_v0gate_stats(V_abc, fs, f, P.v0_nominal); end
    if ~isempty(st0.gnd_svd), fault_type = st0.gnd_svd; end
    if strcmp(P.gnd_mode, 'svdline') && st0.sv(2)/st0.sv(1) < P.gnd_line_ratio && ~isempty(st0.gnd_line)
        fault_type = st0.gnd_line;
    end
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
params.pin = PinB;
params.pinmax = max(PinB);
params.is_kirchhoff = is_kirchoff;
params.deg_bnorm = deg_bnorm;
params.deg_bsmall = deg_bsmall;
params.fallback = false;
params.v0stat = v0stat;
end

%% ========================================================================
function P = defaults_()
P = struct('kirchhoff_thr', 1e-2, 'bnorm_min', 1e-3, 'line_b_abs', 0.01, ...
    'line_ratio', 0.3, 'round_line', 0.05, 'round_circle', 0.99, ...
    'mag_thr', 0.90, 'v_nominal', sqrt(3)/sqrt(2), 'sector_margin', 30, ...
    'biv_dif_thr', 0.3, 'v0_mode', 'samples3', 'v0_nominal', NaN, 'gnd_mode', '2s', ...
    'gnd_line_ratio', 0.05);
end

function P = fill_(Pin)
P = defaults_();
fn = fieldnames(Pin);
for k = 1:numel(fn)
    if ~isfield(P, fn{k})
        error('r2_ga_classifier_param:campo', 'Parametro desconocido: %s', fn{k});
    end
    P.(fn{k}) = Pin.(fn{k});
end
end

function fault_type = classify_by_bivector(biv, dif_thr)
[biv_sort, idx] = sort(abs(biv), "descend");
dif1 = abs(biv_sort(1) - biv_sort(2));
dif2 = abs(biv_sort(2) - biv_sort(3));
if dif1 < dif2 || dif2 > dif_thr
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

function fault_type = classify_by_angle(theta_deg, MARGIN)
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
