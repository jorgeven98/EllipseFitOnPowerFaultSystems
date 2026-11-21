function out = r2_decision_chain(V_abc, t, varargin)
% R2_DECISION_CHAIN  Cadena de decision de referencia (TPWRD ronda 2, T03).
%
% Implementacion unica de la decision completa del paper:
%   ventana deslizante -> ga_fault_classifier_pure (puerta de Kirchhoff
%   incluida) -> puerta adaptativa de redondez -> persistencia de M etiquetas
%   identicas -> declaracion, con todas las marcas temporales al FINAL de la
%   ventana y metricas por evento.
%
% Documentacion completa, con ejemplo: _revision_TPWRD/ronda2/informes/T03_API.md
%
% USO
%   out = r2_decision_chain(V_abc, t)                 % valores por defecto
%   out = r2_decision_chain(V_abc, t, 'truth','A-B','t_on',3.0,'t_off',3.16)
%   out = r2_decision_chain(V_abc, t, opts)           % opts = struct
%   out = r2_decision_chain(scan, [], ...)            % reutiliza out.scan de
%                                                     % una llamada anterior
%                                                     % (no reclasifica)
%   o   = r2_decision_chain('defaults')               % struct de opciones
%
% ENTRADAS
%   V_abc  [N x 3] tensiones fase-tierra, en las MISMAS unidades con que se
%          clasifico en el paper (voltios en simulacion, pu en laboratorio).
%          Los errores de medida se aplican a la senal completa ANTES de
%          llamar a la cadena, conservando la prefalta.
%   t      [N x 1] tiempo en segundos (uniforme).
%
% OPCIONES (nombre, valor por defecto)
%  Ventanas
%   fs            []     Hz; por defecto 1/mediana(diff(t))
%   f             50     Hz
%   W             []     muestras por ventana; por defecto round(fs/(4 f))
%   stride        []     paso; por defecto max(1, floor(W/4)) (solape 75 %)
%   grid_offset   0      la primera ventana empieza en la muestra 1+grid_offset
%   t_range       []     [t0 t1]: solo ventanas con inicio >= t0 y fin <= t1
%                        (la rejilla sigue anclada a la muestra 1+grid_offset)
%   scale         1      la ventana se divide por scale antes de clasificar
%  Clasificador
%   classifier    @ga_fault_classifier_pure   [etq, p] = fn(Vw, fs, f)
%   kirchhoff_thr 0.01   umbral de la puerta de secuencia cero |x^.k^|
%   cls_params    []     struct con otros umbrales internos (ver
%                        r2_ga_classifier_param). Si kirchhoff_thr ~= 0.01 o
%                        cls_params no esta vacio, se usa r2_ga_classifier_param
%   v0_mode       'samples3'  estadistico de la puerta de secuencia cero (T52):
%                        'samples3' (publicada: tres muestras normalizadas) |
%                        'rms_win' | 'rms_nom' | 'rms_max' | 'ph_win' | 'ph_nom' |
%                        'tilt' | 'tilt2' (ver r2_v0gate_stats). El umbral es
%                        kirchhoff_thr. Con otro valor que 'samples3' se usa
%                        r2_ga_classifier_param y out.win.pinmax es ese estadistico
%   v0_nominal    []     pico de fase nominal en las unidades de V_abc (solo
%                        para 'rms_nom' y 'ph_nom')
%   gnd_mode      '2s'   etiqueta de la rama de tierra: '2s' (publicada) | 'svd' |
%                        'svdline' (ver r2_ga_classifier_param)
%  Puerta de redondez
%   gate          'adaptive'  'adaptive' | 'none' | 'fixed'
%   Delta         0.01   margen de la puerta adaptativa
%   K             80     ventanas sanas del nivel ambiente
%   thr_cap       0.99   techo del umbral: thr = min(r_amb - Delta, thr_cap)
%   gate_rc       0.99   umbral si gate = 'fixed'
%   amb_mode      'causal'  'causal' | 'fixed' | 'oracle'
%   r_amb         []     nivel ambiente si amb_mode = 'fixed'
%   amb_init      'assume_healthy'  o un numero (nivel inicial conocido)
%   amb_feed      'healthy_label'   'healthy_label' | 'all'
%   amb_holdoff   0      ventanas sanas consecutivas previas exigidas para
%                        volver a alimentar el nivel ambiente tras una
%                        etiqueta no sana (congelacion)
%  Persistencia
%   M             10     etiquetas identicas consecutivas para declarar
%   healthy_label 'A-B-C'
%  Evento (opcionales; sin t_on solo se devuelve out.summary)
%   truth         ''     clase verdadera: 'A-B','B-C','C-A','A-G','B-G','C-G',
%                        'A-B-G','B-C-G','C-A-G' (se aceptan 'AB','AC','ACG'...)
%   t_on          []     instante de inicio de la falta (s)
%   t_off         []     instante de despeje (s); por defecto fin del registro
%
% SALIDA (struct)
%   out.win      tabla por ventana (ver T03_API.md)
%   out.decl     tabla de declaraciones aseguradas (una fila por declaracion)
%   out.event    struct con las metricas por evento (si se dio t_on)
%   out.summary  struct con tasas sobre todo el tramo procesado
%   out.scan     struct con la clasificacion por ventana, reutilizable
%   out.opts     opciones resueltas
%
% No modifica ni sustituye a ningun script publicado. Autor: T03, nov-2026.

if (ischar(V_abc) || isstring(V_abc)) && strcmpi(V_abc, 'defaults')
    out = defaults_(); return
end
o = parse_(varargin{:});

%% 1. Clasificacion por ventana (o reutilizacion de un escaneo previo)
if isstruct(V_abc)
    S = V_abc;
    req = {'t_start','t_end','label_raw','roundness','angle_deg','semi_major'};
    for k = 1:numel(req)
        if ~isfield(S, req{k})
            error('r2_decision_chain:scan', 'Al escaneo le falta el campo %s', req{k});
        end
    end
    S = complete_scan_(S);
    if isfield(S, 'fs') && ~isempty(S.fs)
        if isempty(o.fs), o.fs = S.fs; end
    elseif ~isempty(o.fs)
        S.fs = o.fs;
    else
        error('r2_decision_chain:scan', 'El escaneo no trae fs: indica la opcion fs.');
    end
else
    t = t(:);
    if size(V_abc, 2) ~= 3 || size(V_abc, 1) ~= numel(t)
        error('r2_decision_chain:dim', 'V_abc debe ser [N x 3] y t [N x 1].');
    end
    if isempty(o.fs), o.fs = 1/median(diff(t)); end
    if isempty(o.W), o.W = round(o.fs/(4*o.f)); end
    if isempty(o.stride), o.stride = max(1, floor(o.W/4)); end
    S = scan_(V_abc, t, o);
end

%% 2. Puerta de redondez, persistencia y declaraciones
D = decide_(S, o);

%% 3. Salidas
out.win = table(S.i_start, S.t_start, S.t_end, S.label_raw, D.label, ...
    S.roundness, S.angle_deg, S.semi_major, S.semi_minor, S.bivector_norm, ...
    S.pinmax, S.is_kirchhoff, S.is_line, S.deg_bnorm, S.deg_bsmall, D.deg_decides, ...
    D.r_amb, D.thr, D.armed, D.fed, D.run, D.fire, ...
    'VariableNames', {'i_start','t_start','t_end','label_raw','label', ...
    'roundness','angle_deg','semi_major','semi_minor','bivector_norm', ...
    'pinmax','is_kirchhoff','is_line','deg_bnorm','deg_bsmall','deg_decides', ...
    'r_amb','thr','armed','fed','run','fire'});
jd = find(D.run == o.M);
out.decl = table(jd, S.t_end(jd), D.label(jd), 'VariableNames', {'j','t_decl','class'});
out.summary = summary_(S, D, true(numel(S.t_end),1), o, '');
out.event = [];
if ~isempty(o.t_on)
    out.event = event_(S, D, out.decl, o);
end
out.scan = S;
o.classifier = func2str(o.classifier);
out.opts = o;
end

%% ========================================================================
function o = defaults_()
o = struct('fs', [], 'f', 50, 'W', [], 'stride', [], 'grid_offset', 0, ...
    't_range', [], 'scale', 1, ...
    'classifier', @ga_fault_classifier_pure, 'kirchhoff_thr', 0.01, 'cls_params', [], ...
    'v0_mode', 'samples3', 'v0_nominal', [], 'gnd_mode', '2s', ...
    'gate', 'adaptive', 'Delta', 0.01, 'K', 80, 'thr_cap', 0.99, 'gate_rc', 0.99, ...
    'amb_mode', 'causal', 'r_amb', [], 'amb_init', 'assume_healthy', ...
    'amb_feed', 'healthy_label', 'amb_holdoff', 0, ...
    'M', 10, 'healthy_label', 'A-B-C', ...
    'truth', '', 't_on', [], 't_off', []);
end

function o = parse_(varargin)
o = defaults_();
if isempty(varargin), args = {};
elseif isstruct(varargin{1})
    s = varargin{1}; fn = fieldnames(s); args = cell(1, 2*numel(fn));
    for k = 1:numel(fn), args{2*k-1} = fn{k}; args{2*k} = s.(fn{k}); end
    args = [args, varargin(2:end)];
else
    args = varargin;
end
if mod(numel(args), 2) ~= 0
    error('r2_decision_chain:opts', 'Las opciones van por pares nombre-valor.');
end
for k = 1:2:numel(args)
    name = char(args{k});
    if ~isfield(o, name)
        error('r2_decision_chain:opts', 'Opcion desconocida: %s', name);
    end
    o.(name) = args{k+1};
end
if ischar(o.classifier) || isstring(o.classifier), o.classifier = str2func(char(o.classifier)); end
o.truth = norm_truth_(o.truth);
% el clasificador de referencia debe estar en la ruta
if exist('ga_fault_classifier_pure', 'file') ~= 2
    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..', 'src', 'fault_classifiers'));
end
end

function tr = norm_truth_(tr)
if isempty(tr), tr = ''; return, end
k = upper(strrep(strrep(char(tr), '-', ''), ' ', ''));
switch k
    case 'AB', tr = 'A-B';
    case {'BC','CB'}, tr = 'B-C';
    case {'CA','AC'}, tr = 'C-A';
    case 'AG', tr = 'A-G';
    case 'BG', tr = 'B-G';
    case 'CG', tr = 'C-G';
    case {'ABG','BAG'}, tr = 'A-B-G';
    case {'BCG','CBG'}, tr = 'B-C-G';
    case {'CAG','ACG'}, tr = 'C-A-G';
    case 'ABC', tr = 'A-B-C';
    case 'ABCG', tr = 'A-B-C-G';
    otherwise, error('r2_decision_chain:truth', 'Clase no reconocida: %s', tr);
end
end

%% ========================================================================
function S = scan_(V, t, o)
% Una pasada del clasificador por la rejilla de ventanas.
N = size(V, 1); W = o.W;
starts = ((1 + o.grid_offset):o.stride:(N - W + 1))';
tol = 0.25/o.fs;
if ~isempty(o.t_range)
    keep = t(starts) >= o.t_range(1) - tol & t(starts + W - 1) <= o.t_range(2) + tol;
    starts = starts(keep);
end
m = numel(starts);
v0_variant = ~strcmp(o.v0_mode, 'samples3'); gnd_variant = ~strcmp(o.gnd_mode, '2s');
use_param = (o.kirchhoff_thr ~= 0.01) || ~isempty(o.cls_params) || v0_variant || gnd_variant;
if use_param
    P = o.cls_params; if isempty(P), P = struct(); end
    P.kirchhoff_thr = o.kirchhoff_thr;
    if v0_variant
        P.v0_mode = o.v0_mode;
        if ~isempty(o.v0_nominal), P.v0_nominal = o.v0_nominal/o.scale; end
    end
    if gnd_variant, P.gnd_mode = o.gnd_mode; end
    cls = @(Vw) r2_ga_classifier_param(Vw, o.fs, o.f, P);
    Pfull = r2_ga_classifier_param('defaults'); fnp = fieldnames(P);
    for k = 1:numel(fnp), Pfull.(fnp{k}) = P.(fnp{k}); end
    bnorm_min = Pfull.bnorm_min;
else
    cls = @(Vw) o.classifier(Vw, o.fs, o.f);
    bnorm_min = 1e-3; v0_variant = false;
end
S.fs = o.fs; S.W = W; S.stride = o.stride;
S.i_start = starts; S.t_start = t(starts); S.t_end = t(starts + W - 1);
S.label_raw = cell(m, 1);
S.roundness = nan(m, 1); S.angle_deg = nan(m, 1);
S.semi_major = nan(m, 1); S.semi_minor = nan(m, 1); S.bivector_norm = nan(m, 1);
S.pinmax = nan(m, 1); S.is_kirchhoff = false(m, 1); S.is_line = false(m, 1);
S.deg_bnorm = false(m, 1); S.deg_bsmall = false(m, 1); S.fallback = false(m, 1);
% indices de las tres muestras de la puerta de Kirchhoff (igual que pure)
q = max(1, round(o.fs/(4*o.f)));
if W < q, i1 = 1; i2 = W; i3 = 1;
else, c = round(W/2); i1 = max(1, c - q); i2 = min(W, i1 + q); i3 = round(abs(i1 - i2)/2);
end
nK = [1 1 1]/sqrt(3);
ws = warning('off', 'all'); cleanup = onCleanup(@() warning(ws));
for j = 1:m
    Vw = V(starts(j):starts(j) + W - 1, :) / o.scale;
    xs = Vw([i1 i2 i3], :); xs = xs ./ sqrt(sum(xs.^2, 2));
    pin = abs(xs * nK');
    S.pinmax(j) = max(pin);
    S.is_kirchhoff(j) = all(pin < o.kirchhoff_thr);
    try
        [ft, p] = cls(Vw);
        S.label_raw{j} = char(ft);
        S.roundness(j) = p.roundness; S.angle_deg(j) = p.angle_deg;
        S.semi_major(j) = p.semi_major; S.semi_minor(j) = p.semi_minor;
        S.bivector_norm(j) = p.bivector_norm; S.is_line(j) = logical(p.is_line);
        S.fallback(j) = isempty(p.x_T);
        if isfield(p, 'is_kirchhoff'), S.is_kirchhoff(j) = p.is_kirchhoff; end
        if v0_variant && isfield(p, 'v0stat'), S.pinmax(j) = p.v0stat; end
        if isfield(p, 'deg_bnorm')
            S.deg_bnorm(j) = p.deg_bnorm; S.deg_bsmall(j) = p.deg_bsmall;
        elseif ~S.fallback(j)
            S.deg_bnorm(j) = p.bivector_norm <= bnorm_min;
            S.deg_bsmall(j) = p.is_line && ~S.deg_bnorm(j) && p.semi_minor == 0;
        end
    catch
        S.label_raw{j} = 'Error';
    end
end
end

function S = complete_scan_(S)
% Rellena campos opcionales de un escaneo externo (p. ej. d1_scan_cache.mat).
m = numel(S.label_raw);
S.label_raw = S.label_raw(:); S.t_start = S.t_start(:); S.t_end = S.t_end(:);
S.roundness = S.roundness(:); S.angle_deg = S.angle_deg(:); S.semi_major = S.semi_major(:);
opt_nan = {'semi_minor','bivector_norm','pinmax'};
for k = 1:numel(opt_nan)
    if ~isfield(S, opt_nan{k}), S.(opt_nan{k}) = nan(m, 1); else, S.(opt_nan{k}) = S.(opt_nan{k})(:); end
end
opt_f = {'is_kirchhoff','is_line','deg_bnorm','deg_bsmall','fallback'};
for k = 1:numel(opt_f)
    if ~isfield(S, opt_f{k}), S.(opt_f{k}) = false(m, 1); else, S.(opt_f{k}) = logical(S.(opt_f{k})(:)); end
end
if ~isfield(S, 'i_start'), S.i_start = nan(m, 1); else, S.i_start = S.i_start(:); end
end

%% ========================================================================
function D = decide_(S, o)
GROUND = {'A-G','B-G','C-G','A-B-G','B-C-G','C-A-G'};
m = numel(S.label_raw);
is_ground = ismember(S.label_raw, GROUND);
invalid = ismember(S.label_raw, {'Error','Unclassified'}) | S.fallback | ~isfinite(S.roundness);
D.label = S.label_raw; D.r_amb = nan(m, 1); D.thr = nan(m, 1);
D.armed = true(m, 1); D.fed = false(m, 1);
vnom = sqrt(3)/sqrt(2); mag_thr = 0.90; sector_margin = 30;
if isstruct(o.cls_params)
    if isfield(o.cls_params, 'v_nominal'), vnom = o.cls_params.v_nominal; end
    if isfield(o.cls_params, 'mag_thr'), mag_thr = o.cls_params.mag_thr; end
    if isfield(o.cls_params, 'sector_margin'), sector_margin = o.cls_params.sector_margin; end
end

    function lab = regate(j, thr)
        % Tierra e invalidas no se tocan; redondez < thr -> L-L por sector;
        % si no, rama redonda del clasificador (A-B-C o A-B-C-G por magnitud).
        if is_ground(j) || invalid(j)
            lab = S.label_raw{j};
        elseif S.roundness(j) < thr
            lab = sector_(S.angle_deg(j), sector_margin);
        elseif S.semi_major(j)/vnom >= mag_thr
            lab = 'A-B-C';
        else
            lab = 'A-B-C-G';
        end
    end

switch lower(o.gate)
    case 'none'
        % etiquetas crudas del clasificador, sin puerta adicional
    case 'fixed'
        D.thr(:) = o.gate_rc;
        for j = 1:m, D.label{j} = regate(j, o.gate_rc); end
    case 'adaptive'
        switch lower(o.amb_mode)
            case 'fixed'
                if isempty(o.r_amb), error('r2_decision_chain:r_amb', 'amb_mode fixed exige r_amb.'); end
                D.r_amb(:) = o.r_amb; D.thr(:) = min(o.r_amb - o.Delta, o.thr_cap);
                for j = 1:m, D.label{j} = regate(j, D.thr(j)); end
            case 'oracle'
                % como la parte D de d1_roundness_gate_design: mediana de las
                % K ultimas ventanas que terminan antes de t_on (no causal en
                % el sentido de que usa el instante de falta conocido)
                if isempty(o.t_on), error('r2_decision_chain:t_on', 'amb_mode oracle exige t_on.'); end
                pre = find(S.t_end <= o.t_on + 0.25/S.fs & isfinite(S.roundness));
                pre = pre(max(1, end - o.K + 1):end);
                ra = median(S.roundness(pre));
                D.r_amb(:) = ra; D.thr(:) = min(ra - o.Delta, o.thr_cap);
                for j = 1:m, D.label{j} = regate(j, D.thr(j)); end
            case 'causal'
                buf = nan(o.K, 1); nb = 0; ib = 0;      % bufer circular
                if isnumeric(o.amb_init)
                    buf(:) = o.amb_init; nb = o.K; armed = true;
                else
                    armed = false;
                end
                run_h = 0;
                for j = 1:m
                    if ~armed
                        % INICIALIZACION: se supone regimen sano; todas las
                        % ventanas con redondez valida alimentan el nivel
                        % ambiente y no se permite declarar
                        D.armed(j) = false; D.thr(j) = o.thr_cap;
                        D.label{j} = regate(j, o.thr_cap);
                        if isfinite(S.roundness(j))
                            ib = mod(ib, o.K) + 1; buf(ib) = S.roundness(j); nb = min(nb + 1, o.K);
                            D.fed(j) = true;
                        end
                        if nb >= o.K, armed = true; end
                    else
                        ra = median(buf(1:nb));
                        D.r_amb(j) = ra; D.thr(j) = min(ra - o.Delta, o.thr_cap);
                        D.label{j} = regate(j, D.thr(j));
                        if strcmp(D.label{j}, o.healthy_label), run_h = run_h + 1; else, run_h = 0; end
                        switch lower(o.amb_feed)
                            case 'healthy_label', feed = run_h > o.amb_holdoff;
                            case 'all', feed = true;
                            otherwise, error('r2_decision_chain:amb_feed', 'amb_feed desconocido.');
                        end
                        if feed && isfinite(S.roundness(j))
                            ib = mod(ib, o.K) + 1; buf(ib) = S.roundness(j);
                            D.fed(j) = true;
                        end
                    end
                end
            otherwise
                error('r2_decision_chain:amb_mode', 'amb_mode desconocido: %s', o.amb_mode);
        end
    otherwise
        error('r2_decision_chain:gate', 'gate desconocido: %s', o.gate);
end

% Persistencia: M etiquetas de falta identicas y consecutivas
healthy = strcmp(D.label, o.healthy_label);
bad = ismember(D.label, {'Error','Unclassified','L-L'});
declarable = D.armed & ~healthy & ~bad;
D.run = zeros(m, 1);
for j = 1:m
    if declarable(j)
        if j > 1 && declarable(j-1) && strcmp(D.label{j}, D.label{j-1})
            D.run(j) = D.run(j-1) + 1;
        else
            D.run(j) = 1;
        end
    end
end
D.fire = D.run >= o.M;
D.healthy = healthy;
D.is_ground = ismember(D.label, GROUND);
D.deg_decides = S.is_line & ~is_ground & ~invalid;   % la rama degenerada fija la etiqueta
end

function s = sector_(theta_deg, MARGIN)
% Identica a classify_by_angle de ga_fault_classifier_pure
th = mod(theta_deg, 180);
if th >= (45 - MARGIN) && th <= (45 + MARGIN), s = 'A-B';
elseif th >= (105 - MARGIN) && th <= (105 + MARGIN), s = 'C-A';
elseif th >= (165 - MARGIN) || th <= 15, s = 'B-C';
else, s = 'L-L';
end
end

%% ========================================================================
function R = summary_(S, D, mask, o, pre)
% Tasas sobre un conjunto de ventanas (secundarias; lo principal es out.event)
n = nnz(mask);
R.([pre 'n_win']) = n;
R.([pre 'n_init']) = nnz(mask & ~D.armed);
R.([pre 'raw_nonhealthy_pct']) = 100*mean(~strcmp(S.label_raw(mask), o.healthy_label));
R.([pre 'nonhealthy_pct']) = 100*mean(~D.healthy(mask));
R.([pre 'ground_pct']) = 100*mean(D.is_ground(mask));
R.([pre 'kirchhoff_fail_pct']) = 100*mean(~S.is_kirchhoff(mask));
R.([pre 'degenerate_pct']) = 100*mean(D.deg_decides(mask));
R.([pre 'deg_bnorm_pct']) = 100*mean(S.deg_bnorm(mask));
R.([pre 'secured_pct']) = 100*mean(D.fire(mask));
R.([pre 'n_decl']) = nnz(mask & D.run == o.M);
jd = find(mask & D.run == o.M);
R.([pre 'decl_classes']) = strjoin(unique(D.label(jd))', '|');
R.([pre 'roundness_med']) = median(S.roundness(mask), 'omitnan');
R.([pre 'pinmax_med']) = median(S.pinmax(mask), 'omitnan');
R.([pre 'pinmax_max']) = max([S.pinmax(mask); NaN]);
end

function E = event_(S, D, decl, o)
tol = 0.25/S.fs; t_on = o.t_on; t_off = o.t_off;
if isempty(t_off), t_off = S.t_end(end); end
truth = o.truth;
pre  = S.t_end <= t_on + tol;                          % solo muestras de prefalta
post = S.t_end > t_on + tol & S.t_end <= t_off + tol;  % contiene falta, antes del despeje
full = S.t_start >= t_on - tol & S.t_end <= t_off + tol;  % integra en la falta
aft  = S.t_end > t_off + tol;
ok   = strcmp(D.label, truth);
okraw = strcmp(S.label_raw, truth);

E.truth = truth; E.t_on = t_on; E.t_off = t_off; E.M = o.M; E.gate = o.gate;
E.n_post = nnz(post); E.n_full = nnz(full);

% primera etiqueta correcta (fin de ventana)
j1 = find(post & ok, 1);
E.first_label_ms = ms_(S.t_end, j1, t_on);
j1r = find(post & okraw, 1);
E.first_label_raw_ms = ms_(S.t_end, j1r, t_on);

% declaraciones aseguradas durante la falta
dp = decl(decl.t_decl > t_on + tol & decl.t_decl <= t_off + tol, :);
E.declared = ~isempty(dp);
if E.declared
    E.decl_class = dp.class{1};
    E.decl_ms = 1e3*(dp.t_decl(1) - t_on);
    if strcmp(dp.class{1}, truth), E.outcome = 'correcta'; else, E.outcome = 'erronea'; end
else
    E.decl_class = ''; E.decl_ms = NaN; E.outcome = 'sin_declaracion';
end
ic = find(strcmp(dp.class, truth), 1);
if isempty(ic), E.correct_decl_ms = NaN; else, E.correct_decl_ms = 1e3*(dp.t_decl(ic) - t_on); end
E.n_decl_correct = nnz(strcmp(dp.class, truth));
E.n_decl_wrong = nnz(~strcmp(dp.class, truth));
E.wrong_classes = strjoin(unique(dp.class(~strcmp(dp.class, truth)))', '|');
if isempty(ic), E.n_wrong_before_correct = E.n_decl_wrong;
else, E.n_wrong_before_correct = nnz(~strcmp(dp.class(1:ic-1), truth)); end

% declaraciones fuera de la falta
dpre = decl(decl.t_decl <= t_on + tol, :);
E.n_decl_pre = height(dpre); E.pre_classes = strjoin(unique(dpre.class)', '|');
daft = decl(decl.t_decl > t_off + tol, :);
E.n_decl_after = height(daft); E.after_classes = strjoin(unique(daft.class)', '|');

% racha correcta maxima durante la falta
r = 0; rmax = 0;
for j = find(post)'
    if ok(j), r = r + 1; rmax = max(rmax, r); else, r = 0; end
end
E.max_run_correct = rmax;

% secundarias, sobre las ventanas integras en la falta
E.retention_pct = 100*mean(ok(full));
E.retention_raw_pct = 100*mean(okraw(full));
E.full_healthy_pct = 100*mean(D.healthy(full));
E.full_ground_pct = 100*mean(D.is_ground(full));
E.full_other_pct = 100*mean(~ok(full) & ~D.healthy(full) & ~D.is_ground(full));
E.full_kirchhoff_fail_pct = 100*mean(~S.is_kirchhoff(full));
E.full_degenerate_pct = 100*mean(D.deg_decides(full));
E.full_secured_correct = nnz(full & D.fire & ok);
E.full_secured_wrong = nnz(full & D.fire & ~ok);
lf = D.label(full);
if isempty(lf), E.full_modal_label = '';
else, [u, ~, ju] = unique(lf); cnt = accumarray(ju, 1); [~, im] = max(cnt); E.full_modal_label = u{im};
end
E.full_roundness_med = median(S.roundness(full), 'omitnan');
E.full_angle_med = median(S.angle_deg(full), 'omitnan');
E.full_pinmax_med = median(S.pinmax(full), 'omitnan');
E.full_pinmax_max = max([S.pinmax(full); NaN]);
E.full_thr_med = median(D.thr(full), 'omitnan');

% prefalta y posfalta (seguridad)
P = summary_(S, D, pre, o, 'pre_');
fn = {'pre_n_win','pre_n_init','pre_raw_nonhealthy_pct','pre_nonhealthy_pct','pre_secured_pct','pre_degenerate_pct'};
for k = 1:numel(fn), E.(fn{k}) = P.(fn{k}); end
E.aft_n_win = nnz(aft);
end

function v = ms_(tend, j, t_on)
if isempty(j), v = NaN; else, v = 1e3*(tend(j) - t_on); end
end
