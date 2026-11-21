function r2_degenerate_security(opts, tag, blocks, NPOOL)
% R2_DEGENERATE_SECURITY  Seguridad de la rama degenerada con persistencia (T41, R1.4 y E3).
%
%   r2_degenerate_security()                    % cadena de referencia, todos los bloques
%   r2_degenerate_security(opts)                % opts = struct con opciones de r2_decision_chain
%   r2_degenerate_security(opts, 'G2')          % salidas en resultados/T41/G2
%   r2_degenerate_security(opts, tag, {'A','C'})% solo algunos bloques
%   r2_degenerate_security(opts, tag, blocks, NPOOL) % ventanas del bloque B (por defecto 2e6;
%                                                    % las cifras de T41 usan 2e7: unos 40 min)
%
% TODO lo que decide (clasificador, puerta de secuencia cero, puerta de
% redondez, persistencia) sale de r2_decision_chain con las opciones de
% `opts`. Para repetir el estudio con otra variante del clasificador basta
% una llamada, p. ej.
%   r2_degenerate_security(struct('kirchhoff_thr', 0.02), 'kthr002')
%   r2_degenerate_security(struct('classifier', @mi_clasificador), 'G1')
%
% Bloques
%   A  Escenarios con exposicion de al menos 10 s por escenario y convertidor
%      (Pex y DG1), en voltios (como se publico) y en pu: sano con ruido de
%      60 a 0 dB, linea sin tension (ruido, cero exacto, ruido cuantizado,
%      tension inducida residual), hueco trifasico equilibrado profundo con
%      ruido, perdida de tension de medida en una y dos fases, y control
%      positivo con faltas francas sinteticas.
%   B  Persistencia frente a etiquetas aleatorias: secuencias largas de
%      ventanas (linea sin tension, sano a 30 dB, perdida de dos fases, sano
%      a 40 dB sin puerta) con probabilidades de clase, coincidencia entre
%      ventanas consecutivas, histograma de rachas, recuento directo de
%      rachas de M y extrapolacion (independencia y cadena de Markov).
%   C  Suite de 320 eventos de seguridad del proyecto DGC (8 familias),
%      generada por r2_security_events_gen.
%   D  Longitud real de la trayectoria en las 24 faltas del paper (el
%      clasificador fuerza el semieje mayor en la rama degenerada).
%
% Salidas (resultados/T41[/tag]): escenarios.csv, declaraciones.csv,
% persistencia.csv, persistencia_rachas.csv, eventos_seguridad.csv,
% eventos_seguridad_resumen.csv, longitud_linea_faltas.csv, log.txt.
% Semillas fijas. No modifica ningun archivo existente.
if nargin < 1 || isempty(opts), opts = struct(); end
if nargin < 2 || isempty(tag), tag = ''; end
if nargin < 3 || isempty(blocks), blocks = {'A','B','C','D'}; end
if nargin < 4 || isempty(NPOOL), NPOOL = 2e6; end
here = fileparts(mfilename('fullpath')); RAIZ = fullfile(here, '..');
addpath(here); addpath(fullfile(RAIZ, 'src', 'fault_classifiers')); addpath(fullfile(RAIZ, 'src', 'utils'));
OUT = fullfile(RAIZ, '_revision_TPWRD', 'ronda2', 'resultados', 'T41', tag);
if ~exist(OUT, 'dir'), mkdir(OUT); end
logf = fullfile(OUT, sprintf('log_%s.txt', strjoin(blocks, '')));
if exist(logf, 'file'), delete(logf); end
diary(logf); cobj = onCleanup(@() diary('off'));
fprintf('R2_DEGENERATE_SECURITY  %s  MATLAB %s\nbloques: %s   tag: "%s"\nopciones de la cadena distintas de las de referencia:\n', ...
    datestr(now, 31), version, strjoin(blocks, ','), tag); %#ok<TNOW1,DATST>
disp(opts);
warning('off', 'all');

C.opts = opts; C.fs = 4000; C.f = 50; C.W = 20; C.STRIDE = 5;
d0 = r2_decision_chain('defaults'); C.M = d0.M; if isfield(opts, 'M'), C.M = opts.M; end
C.GR = {'A-G','B-G','C-G','A-B-G','B-C-G','C-A-G'};
C.LABS = {'A-B-C','A-B','B-C','C-A','A-G','B-G','C-G','A-B-G','B-C-G','C-A-G','A-B-C-G'};
C.OUT = OUT;
[Tc, info] = r2_sim_cases(); C.info = info; C.Tc = Tc;
base = load(fullfile(info.dir, 'DataGridRL_AB.mat'));
CONV = info.converters;
for c = 1:numel(CONV)
    [V, t] = extract_converter_data(base, CONV{c});
    C.Npk(c) = max(abs(V(t >= 2.4 & t < 2.8, :)), [], 'all');
    seg = t >= 0.4 & t < 2.95;                   % tramo sano real (2.55 s), sin el arranque
    C.Vh{c} = V(seg, :); C.th{c} = t(seg);
end
C.CONV = CONV; C.T_EV = 0.65;                    % inicio del evento: 0.25 s (200 ventanas) de sano previo

if ismember('A', blocks), block_A(C); end
if ismember('B', blocks), block_B(C, NPOOL); end
if ismember('C', blocks), block_C(C); end
if ismember('D', blocks), block_D(C); end
fprintf('\nFIN r2_degenerate_security (%s).\n', datestr(now, 31)); %#ok<TNOW1,DATST>
end

%% =====================================================================
function block_A(C)
fprintf('\n================ BLOQUE A: escenarios (exposicion >= 10 s por escenario y convertidor) ================\n');
NREP = 6;                                       % 6 x 1.90 s de exposicion tras el inicio = 11.4 s (14.7 s en los escenarios sin evento)
SC = {};                                        % nombre, familia, tipo, parametros, n_rep
for snr = [60 50 40 30 20 10 0]
    SC(end+1, :) = {sprintf('sano + ruido %d dB', snr), 'sano con ruido', 'noise', struct('snr', snr), NREP}; %#ok<AGROW>
end
for sg = [1e-4 1e-3 1e-2 5e-2]
    SC(end+1, :) = {sprintf('linea sin tension, ruido %g del pico', sg), 'linea sin tension', 'dead', struct('sig', sg, 'lsb', 0), NREP}; %#ok<AGROW>
end
SC(end+1, :) = {'linea sin tension, cero exacto', 'linea sin tension', 'dead', struct('sig', 0, 'lsb', 0), 1};
LSB = 1/2048;                                   % conversor de 12 bits con fondo de escala igual al pico nominal
SC(end+1, :) = {'linea sin tension, ruido de 0.3 LSB cuantizado (12 bits)', 'linea sin tension', 'dead', struct('sig', 0.3*LSB, 'lsb', LSB), NREP};
SC(end+1, :) = {'linea sin tension, ruido de 1 LSB cuantizado (12 bits)', 'linea sin tension', 'dead', struct('sig', 1.0*LSB, 'lsb', LSB), NREP};
SC(end+1, :) = {'linea sin tension, ruido de 3 LSB cuantizado (12 bits)', 'linea sin tension', 'dead', struct('sig', 3.0*LSB, 'lsb', LSB), NREP};
U = {'comun (1,1,1)', [1 1 1]; 'una fase (1,0,0)', [1 0 0]; 'entre fases (1,-1,0)', [1 -1 0]};
for iu = 1:size(U, 1)
    for ep = [0.005 0.02 0.05]
        SC(end+1, :) = {sprintf('linea sin tension, inducida %s de %g %% + ruido 0.1 %%', U{iu, 1}, 100*ep), 'linea sin tension con tension inducida', ...
            'induced', struct('u', U{iu, 2}, 'eps', ep, 'sig', 1e-3), NREP}; %#ok<AGROW>
    end
end
for r = [0.5 0.2 0.1 0.05 0.02 0.01]
    for snr = [60 40 30]
        SC(end+1, :) = {sprintf('hueco trifasico equilibrado, tension residual %g pu, ruido %d dB', r, snr), 'hueco trifasico equilibrado', ...
            'dip3', struct('r', r, 'snr', snr), NREP}; %#ok<AGROW>
    end
end
PH = {'A', 1; 'B', 2; 'C', 3};
for ip = 1:3
    SC(end+1, :) = {sprintf('perdida de tension en %s (ruido residual 0.002)', PH{ip, 1}), 'perdida de tension, una fase', 'lop', struct('ph', PH{ip, 2}, 'sig', 0.002, 'snr', Inf), NREP}; %#ok<AGROW>
end
SC(end+1, :) = {'perdida de tension en A (cero exacto)', 'perdida de tension, una fase', 'lop', struct('ph', 1, 'sig', 0, 'snr', Inf), 1};
SC(end+1, :) = {'perdida de tension en A (ruido residual 0.01)', 'perdida de tension, una fase', 'lop', struct('ph', 1, 'sig', 0.01, 'snr', Inf), NREP};
SC(end+1, :) = {'perdida de tension en A (residual 0.002) + ruido 40 dB en todo', 'perdida de tension, una fase', 'lop', struct('ph', 1, 'sig', 0.002, 'snr', 40), NREP};
PP = {'B y C', [2 3]; 'A y B', [1 2]; 'C y A', [1 3]};
for ip = 1:3
    SC(end+1, :) = {sprintf('perdida de tension en %s (ruido residual 0.002)', PP{ip, 1}), 'perdida de tension, dos fases', 'lop', struct('ph', PP{ip, 2}, 'sig', 0.002, 'snr', Inf), NREP}; %#ok<AGROW>
end
SC(end+1, :) = {'perdida de tension en B y C (cero exacto)', 'perdida de tension, dos fases', 'lop', struct('ph', [2 3], 'sig', 0, 'snr', Inf), 1};
SC(end+1, :) = {'perdida de tension en B y C (ruido residual 0.01)', 'perdida de tension, dos fases', 'lop', struct('ph', [2 3], 'sig', 0.01, 'snr', Inf), NREP};
SC(end+1, :) = {'perdida de tension en B y C (residual 0.002) + ruido 40 dB en todo', 'perdida de tension, dos fases', 'lop', struct('ph', [2 3], 'sig', 0.002, 'snr', 40), NREP};
BO = {'A-B', [1 2 3]; 'B-C', [2 3 1]; 'C-A', [3 1 2]};
for ib = 1:3
    SC(end+1, :) = {sprintf('control positivo: %s franca sintetica', BO{ib, 1}), 'control positivo', 'bolted', struct('perm', BO{ib, 2}, 'truth', BO{ib, 1}), 1}; %#ok<AGROW>
end

rowsS = {}; rowsD = {};
UN = {'voltios', 'pu'};
for c = 1:numel(C.CONV)
    V0 = C.Vh{c}; tl = C.th{c}; Npk = C.Npk(c); on = (tl - tl(1)) >= C.T_EV; t_on = tl(find(on, 1));
    vnom = sqrt(3/2)*Npk;
    for i = 1:size(SC, 1)
        rs = RandStream('mt19937ar', 'Seed', 4100 + 100*c + i);
        A = acc_init(C);
        A2 = A;
        for rep = 1:SC{i, 5}
            [V, has_ev] = gen_(SC{i, 3}, SC{i, 4}, V0, tl, on, Npk, rs, C.f);
            if has_ev, t_ev = t_on; else, t_ev = []; end
            [aw, ~] = r2_window_amplitude(V, (1:C.STRIDE:(size(V, 1) - C.W + 1))', C.W, C.fs, C.f);
            for u = 1:2
                if u == 1, sc_ = 1; else, sc_ = Npk; end
                o = r2_decision_chain(V, tl, r2_opts_nominal(C.opts, Npk), 'scale', sc_);
                [R, dd] = analyze_(o, t_ev, aw/vnom, C);
                if u == 1, A = acc_add(A, R); else, A2 = acc_add(A2, R); end
                for q = 1:size(dd, 1)
                    rowsD(end+1, :) = [{SC{i, 1}, SC{i, 2}, C.CONV{c}, UN{u}, rep}, dd(q, :)]; %#ok<AGROW>
                end
            end
        end
        for u = 1:2
            if u == 1, Au = A; else, Au = A2; end
            rowsS(end+1, :) = [{SC{i, 1}, SC{i, 2}, C.CONV{c}, UN{u}}, acc_row(Au, C)]; %#ok<AGROW>
            fprintf('%-74s %-3s %-7s | %5.1f s | linea %7.3f %% decide %7.3f %% biv.nulo %7.3f %% | no sana %5.1f %% tierra %5.1f %% | racha max %3d | decl %3d en %d/%d reg: %s\n', ...
                SC{i, 1}, C.CONV{c}, UN{u}, rowsS{end, 6}, rowsS{end, 8}, rowsS{end, 9}, rowsS{end, 10}, rowsS{end, 13}, rowsS{end, 14}, rowsS{end, 15}, rowsS{end, 16}, rowsS{end, 17}, rowsS{end, 5}, rowsS{end, 18});
        end
    end
end
vn = [{'escenario','familia','convertidor','unidades'}, acc_names(C)];
TS = cell2table(rowsS, 'VariableNames', vn);
writetable(TS, fullfile(C.OUT, 'escenarios.csv'), 'QuoteStrings', true);
TD = cell2table(cell(0, 12), 'VariableNames', {'escenario','familia','convertidor','unidades','rep','t_decl_ms','clase','via','frac_linea_en_racha','frac_biv_nulo_en_racha','semieje_mayor_rel','redondez_med'});
if ~isempty(rowsD), TD = cell2table(rowsD, 'VariableNames', TD.Properties.VariableNames); end
writetable(TD, fullfile(C.OUT, 'declaraciones.csv'), 'QuoteStrings', true);
fprintf('Bloque A: %d filas de escenario y %d declaraciones.\n', height(TS), height(TD));
end

function [V, has_ev] = gen_(kind, p, V0, tl, on, Npk, rs, f)
% Construye el registro del escenario en voltios. El evento empieza en `on`.
V = V0; has_ev = true; n = nnz(on);
switch kind
    case 'noise'
        V = V0 + (Npk/sqrt(2))*10^(-p.snr/20)*randn(rs, size(V0)); has_ev = false;
    case 'dead'
        X = p.sig*Npk*randn(rs, n, 3);
        if p.lsb > 0, X = round(X/(p.lsb*Npk))*(p.lsb*Npk); end
        V(on, :) = X;
    case 'induced'
        ph = 2*pi*rand(rs);
        V(on, :) = p.eps*Npk*cos(2*pi*f*tl(on) + ph)*p.u + p.sig*Npk*randn(rs, n, 3);
    case 'dip3'
        V(on, :) = p.r*V0(on, :);
        V = V + (Npk/sqrt(2))*10^(-p.snr/20)*randn(rs, size(V0));
    case 'lop'
        V(on, p.ph) = p.sig*Npk*randn(rs, n, numel(p.ph));
        if isfinite(p.snr), V = V + (Npk/sqrt(2))*10^(-p.snr/20)*randn(rs, size(V0)); end
    case 'bolted'
        % falta franca entre las fases perm(1) y perm(2): ambas valen -v3/2
        k1 = p.perm(1); k2 = p.perm(2); k3 = p.perm(3);
        V(on, k1) = -V0(on, k3)/2; V(on, k2) = -V0(on, k3)/2;
    otherwise
        error('tipo de escenario desconocido: %s', kind);
end
end

function [R, dd] = analyze_(o, t_ev, a_rel, C)
% Recuentos de un registro tras el inicio del evento (o en todo el registro armado).
w = o.win; M = o.opts.M;
if isempty(t_ev), m = w.armed; dm = o.decl; t0 = w.t_start(1);
else, m = w.armed & w.t_start >= t_ev - 1e-9; dm = o.decl(o.decl.t_decl > t_ev, :); t0 = t_ev;
end
R.n = nnz(m); R.exp_s = nnz(m)*o.scan.stride/o.scan.fs;
R.is_line = nnz(w.is_line(m)); R.decides = nnz(w.deg_decides(m)); R.bnorm = nnz(w.deg_bnorm(m)); R.bsmall = nnz(w.deg_bsmall(m));
R.nonh = nnz(~strcmp(w.label(m), 'A-B-C')); R.ground = nnz(ismember(w.label(m), C.GR));
R.rawnonh = nnz(~strcmp(w.label_raw(m), 'A-B-C'));
R.fire = nnz(w.fire(m)); R.maxrun = max([w.run(m); 0]);
R.cnt = zeros(1, numel(C.LABS) + 1);
[tf, loc] = ismember(w.label(m), C.LABS);
R.cnt(1:numel(C.LABS)) = accumarray(loc(tf), 1, [numel(C.LABS) 1])';
R.cnt(end) = nnz(~tf);
R.ndecl = height(dm); R.rec_decl = double(height(dm) > 0);
R.classes = dm.class;
if isempty(dm), R.first_ms = NaN; else, R.first_ms = 1e3*(dm.t_decl(1) - t0); end
dd = cell(height(dm), 7); R.paths = cell(height(dm), 1);
for q = 1:height(dm)
    j = dm.j(q); run = (j - M + 1):j;
    fl = mean(w.is_line(run)); fb = mean(w.deg_bnorm(run));
    if ismember(dm.class{q}, C.GR), via = 'tierra';
    elseif strcmp(dm.class{q}, 'A-B-C-G'), via = 'trifasica por magnitud';
    elseif fl == 1, via = 'rama degenerada';
    elseif fl == 0, via = 'puerta de redondez';
    else, via = 'mixta (degenerada y redondez)';
    end
    R.paths{q} = via;
    ia = (w.i_start(run) - 1)/o.scan.stride + 1;              % indice de ventana en la rejilla completa
    dd(q, :) = {1e3*(dm.t_decl(q) - t0), dm.class{q}, via, fl, fb, median(a_rel(ia), 'omitnan'), median(w.roundness(run), 'omitnan')};
end
end

function A = acc_init(C)
A = struct('nrec', 0, 'n', 0, 'exp_s', 0, 'is_line', 0, 'decides', 0, 'bnorm', 0, 'bsmall', 0, 'nonh', 0, 'rawnonh', 0, 'ground', 0, ...
    'fire', 0, 'maxrun', 0, 'cnt', zeros(1, numel(C.LABS) + 1), 'ndecl', 0, 'rec_decl', 0, 'classes', {{}}, 'paths', {{}}, 'first_ms', []);
end
function A = acc_add(A, R)
A.nrec = A.nrec + 1;
fn = {'n','exp_s','is_line','decides','bnorm','bsmall','nonh','rawnonh','ground','fire','ndecl','rec_decl'};
for k = 1:numel(fn), A.(fn{k}) = A.(fn{k}) + R.(fn{k}); end
A.cnt = A.cnt + R.cnt; A.maxrun = max(A.maxrun, R.maxrun);
A.classes = [A.classes; R.classes]; A.paths = [A.paths; R.paths]; A.first_ms(end+1) = R.first_ms;
end
function n = acc_names(C)
lab = strcat('pct_', strrep(C.LABS, '-', ''));
n = [{'n_registros','exposicion_s','n_ventanas','linea_pct','rama_degenerada_decide_pct','biv_nulo_pct','semieje_menor_abs_pct', ...
      'etiqueta_cruda_no_sana_pct','etiqueta_no_sana_pct','etiqueta_tierra_pct','racha_max','n_declaraciones','registros_con_declaracion', ...
      'clases_declaradas','vias_declaracion','ventanas_aseguradas_pct','decl_por_10s','primera_decl_ms_mediana','primera_decl_ms_max'}, lab, {'pct_otra'}];
end
function r = acc_row(A, C) %#ok<INUSD>
p = @(x) 100*x/max(A.n, 1);
r = [{A.nrec, A.exp_s, A.n, p(A.is_line), p(A.decides), p(A.bnorm), p(A.bsmall), p(A.rawnonh), p(A.nonh), p(A.ground), A.maxrun, A.ndecl, A.rec_decl, ...
      cnt_str(A.classes), cnt_str(A.paths), p(A.fire), 10*A.ndecl/max(A.exp_s, eps), median(A.first_ms, 'omitnan'), max([A.first_ms(~isnan(A.first_ms)), NaN], [], 'omitnan')}, ...
      num2cell(100*A.cnt/max(A.n, 1))];
end
function s = cnt_str(c)
if isempty(c), s = ''; return, end
[u, ~, j] = unique(c); n = accumarray(j, 1);
s = strjoin(arrayfun(@(k) sprintf('%s x%d', u{k}, n(k)), 1:numel(u), 'UniformOutput', false), '; ');
end

%% =====================================================================
function block_B(C, NPOOL)
fprintf('\n================ BLOQUE B: persistencia frente a etiquetas aleatorias ================\n');
% Muestras que usa el clasificador para la etiqueta de tierra
q = round(C.fs/(4*C.f)); cc = round(C.W/2); i1 = max(1, cc - q); i2 = min(C.W, i1 + q); i3 = round(abs(i1 - i2)/2);
fprintf('Ventana de %d muestras, paso %d: el bivector usa las muestras %d y %d; la puerta de secuencia cero, %d, %d y %d.\n', C.W, C.STRIDE, i1, i2, i1, i2, i3);
sh = 0; for k = 1:20, sh = sh + numel(intersect([i1 i2] + C.STRIDE*k, [i1 i2])); end
fprintf('Muestras del bivector compartidas con cualquiera de las 20 ventanas siguientes: %d (de las %d que se solapan en ventanas consecutivas).\n', sh, C.W - C.STRIDE);
c = 1; V0 = C.Vh{c}; tl = C.th{c}; Npk = C.Npk(c); nseg = size(V0, 1);
wps = C.fs/C.STRIDE;
% nombre, tipo, parametros, n ventanas, puerta
P = {'linea sin tension (ruido blanco puro)', 'purenoise', struct('sig', 1e-2), NPOOL, 'adaptive'
     'sano Pex + ruido 30 dB', 'noise', struct('snr', 30), 2e5, 'adaptive'
     'sano Pex + ruido 40 dB, cadena completa', 'noise', struct('snr', 40), 2e5, 'adaptive'
     'sano Pex + ruido 40 dB, sin puerta adaptativa', 'noise', struct('snr', 40), 2e5, 'none'
     'perdida de tension en B y C (residual 0.002)', 'lop', struct('ph', [2 3], 'sig', 0.002, 'snr', Inf), 2e5, 'adaptive'
     'hueco trifasico residual 0.02 pu + ruido 40 dB', 'dip3', struct('r', 0.02, 'snr', 40), 2e5, 'adaptive'};
rows = {}; rowsR = {};
for i = 1:size(P, 1)
    rs = RandStream('mt19937ar', 'Seed', 4500 + i);
    lab = {}; chunk = []; isl = []; tic; nc = 0; ndecl = 0; dcls = {};
    while numel(lab) < P{i, 4}
        nc = nc + 1;
        if strcmp(P{i, 2}, 'purenoise')
            n = 400000; tt = (0:n-1)'/C.fs; V = P{i, 3}.sig*Npk*randn(rs, n, 3);
            o = r2_decision_chain(V, tt, r2_opts_nominal(C.opts, Npk), 'gate', P{i, 5}, 'amb_init', 0.99);   % sin tramo sano previo: se arma desde la primera ventana
            keep = true(height(o.win), 1);
        else
            on = true(nseg, 1); on(1:round(C.T_EV*C.fs)) = false;
            [V, has_ev] = gen_(P{i, 2}, P{i, 3}, V0, tl, on, Npk, rs, C.f);
            o = r2_decision_chain(V, tl, r2_opts_nominal(C.opts, Npk), 'gate', P{i, 5});
            if has_ev, keep = o.win.armed & o.win.t_start >= tl(find(on, 1)) - 1e-9; else, keep = o.win.armed; end
        end
        lab = [lab; o.win.label(keep)]; chunk = [chunk; nc*ones(nnz(keep), 1)]; isl = [isl; o.win.is_line(keep)]; %#ok<AGROW>
        kj = find(keep); dj = o.decl(ismember(o.decl.j, kj) & ismember(o.decl.j - C.M + 1, kj), :);
        ndecl = ndecl + height(dj); dcls = [dcls; dj.class]; %#ok<AGROW>
    end
    S = persist_stats(lab, chunk, C);
    fprintf('\n-- %s (puerta %s): %d ventanas (%.0f s, %d tramos, %.0f s de calculo)\n', P{i, 1}, P{i, 5}, S.n, S.n/wps, nc, toc);
    fprintf('   clases: %s\n', S.dist);
    fprintf('   coincidencia con la ventana anterior: %.4f (independencia: %.4f); a 2, 3 y 4 ventanas: %.4f %.4f %.4f\n', S.agree(1), S.sum_p2, S.agree(2), S.agree(3), S.agree(4));
    fprintf('   coincidencia a 1 ventana entre etiquetas de tierra: %.4f (indep. %.4f, n = %d); entre etiquetas sin tierra: %.4f (indep. %.4f, n = %d)\n', S.agree_g, S.indep_g, S.n_g, S.agree_ll, S.indep_ll, S.n_ll);
    fprintf('   rachas de M = %d: observadas %d (cadena: %d declaraciones: %s); esperadas con independencia %.3g; con cadena de Markov %.3g\n', C.M, S.n_runM, ndecl, cnt_str(dcls), S.exp_iid, S.exp_markov);
    fprintf('   tasa por hora: observada %.3g (cota superior al 95 %%: %.3g); independencia %.3g; Markov %.3g\n', S.n_runM/S.n*wps*3600, S.ub95/S.n*wps*3600, S.exp_iid/S.n*wps*3600, S.exp_markov/S.n*wps*3600);
    rows(end+1, :) = {P{i, 1}, P{i, 5}, S.n, S.n/wps, S.dist, S.sum_p2, S.agree(1), S.agree(2), S.agree(3), S.agree(4), S.agree_g, S.indep_g, S.n_g, S.agree_ll, S.indep_ll, S.n_ll, ...
        S.maxrun, S.n_runM, ndecl, cnt_str(dcls), S.ub95, S.exp_iid, S.exp_markov, S.n_runM/S.n, S.ub95/S.n, S.exp_iid/S.n, S.exp_markov/S.n, ...
        S.n_runM/S.n*wps*3600, S.ub95/S.n*wps*3600, S.exp_iid/S.n*wps*3600, S.exp_markov/S.n*wps*3600, 100*mean(isl)}; %#ok<AGROW>
    for L = 1:numel(S.runs_ge)
        rowsR(end+1, :) = {P{i, 1}, P{i, 5}, L, S.runs_eq(L), S.runs_ge(L), S.runs_ge_iid(L), S.runs_ge_markov(L)}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'secuencia','puerta','n_ventanas','duracion_s','distribucion_etiquetas','suma_p2','coinc_lag1','coinc_lag2','coinc_lag3','coinc_lag4', ...
    'coinc_lag1_tierra','indep_tierra','n_pares_tierra','coinc_lag1_sin_tierra','indep_sin_tierra','n_pares_sin_tierra', ...
    'racha_max','n_rachas_M_observadas','n_declaraciones_cadena','clases_declaradas','cota95_rachas_M','esperadas_independencia','esperadas_markov', ...
    'prob_por_ventana_observada','prob_por_ventana_cota95','prob_por_ventana_independencia','prob_por_ventana_markov', ...
    'por_hora_observada','por_hora_cota95','por_hora_independencia','por_hora_markov','linea_pct'});
writetable(T, fullfile(C.OUT, 'persistencia.csv'), 'QuoteStrings', true);
TR = cell2table(rowsR, 'VariableNames', {'secuencia','puerta','L','n_rachas_longitud_L','n_rachas_al_menos_L','esperadas_al_menos_L_independencia','esperadas_al_menos_L_markov'});
writetable(TR, fullfile(C.OUT, 'persistencia_rachas.csv'), 'QuoteStrings', true);
end

function S = persist_stats(lab, chunk, C)
% Estadistica de rachas de etiquetas de falta identicas (misma definicion que la cadena).
n = numel(lab); S.n = n; M = C.M;
[u, ~, id] = unique(lab); K = numel(u); p = accumarray(id, 1)/n;
[ps, od] = sort(p, 'descend');
S.dist = strjoin(arrayfun(@(k) sprintf('%s %.2f%%', u{od(k)}, 100*ps(k)), 1:min(K, 8), 'UniformOutput', false), '; ');
decl = ~ismember(u, {'A-B-C', 'Error', 'Unclassified', 'L-L'});      % clases que pueden declarar
isg = ismember(u, C.GR);
S.sum_p2 = sum(p.^2);
for k = 1:4
    same = chunk(1+k:end) == chunk(1:end-k);
    S.agree(k) = mean(id(1+k:end) == id(1:end-k) & same)/mean(same);
end
same = chunk(2:end) == chunk(1:end-1); a = id(2:end); b = id(1:end-1);
gg = same & isg(a) & isg(b); S.n_g = nnz(gg);
S.agree_g = mean(a(gg) == b(gg)); pg = p(isg)/max(sum(p(isg)), eps); S.indep_g = sum(pg.^2);
ll = same & decl(a) & decl(b) & ~isg(a) & ~isg(b); S.n_ll = nnz(ll);
S.agree_ll = mean(a(ll) == b(ll)); pl = p(decl & ~isg)/max(sum(p(decl & ~isg)), eps); S.indep_ll = sum(pl.^2);
if S.n_g == 0, S.agree_g = NaN; S.indep_g = NaN; end
if S.n_ll == 0, S.agree_ll = NaN; S.indep_ll = NaN; end
% rachas
LMAX = max(12, M + 2); runs_eq = zeros(LMAX, 1); r = 0; maxrun = 0; n_runM = 0;
for j = 1:n
    if decl(id(j))
        if j > 1 && chunk(j) == chunk(j-1) && id(j) == id(j-1), r = r + 1; else
            if r > 0, runs_eq(min(r, LMAX)) = runs_eq(min(r, LMAX)) + 1; end
            r = 1;
        end
        if r == M, n_runM = n_runM + 1; end
        maxrun = max(maxrun, r);
    else
        if r > 0, runs_eq(min(r, LMAX)) = runs_eq(min(r, LMAX)) + 1; end
        r = 0;
    end
end
if r > 0, runs_eq(min(r, LMAX)) = runs_eq(min(r, LMAX)) + 1; end
S.runs_eq = runs_eq; S.runs_ge = flipud(cumsum(flipud(runs_eq))); S.maxrun = maxrun; S.n_runM = n_runM;
% prediccion con independencia: P(racha que alcanza L en la ventana j) = p^L (1 - p)
pd = p(decl);
S.runs_ge_iid = arrayfun(@(L) n*sum(pd.^L.*(1 - pd)), (1:LMAX)');
% cadena de Markov de primer orden estimada de la secuencia
Tm = accumarray([b(same), a(same)], 1, [K K]); rowsum = sum(Tm, 2); Tm = Tm./max(rowsum, 1);
pin = zeros(K, 1);                         % probabilidad de empezar racha de la clase c
for cidx = 1:K
    if decl(cidx), pin(cidx) = sum(p.*Tm(:, cidx)) - p(cidx)*Tm(cidx, cidx); end
end
tcc = diag(Tm);
S.runs_ge_markov = arrayfun(@(L) n*sum(pin(decl).*tcc(decl).^(L-1)), (1:LMAX)');
S.exp_iid = S.runs_ge_iid(M); S.exp_markov = S.runs_ge_markov(M);
S.ub95 = gammaincinv(0.975, n_runM + 1);   % cota superior de Poisson al 97.5 % unilateral
end

%% =====================================================================
function block_C(C)
fprintf('\n================ BLOQUE C: suite de eventos de seguridad del proyecto DGC ================\n');
E = r2_security_events_gen(40);
fprintf('%d eventos generados (rng(7), mismo orden de sorteos que exp_security_events.m)\n', numel(E));
UN = {'pu', 'voltios'}; SCL = [1, 1/C.Npk(1)];          % el generador da pu; en "voltios" se multiplica por el pico de Pex
rows = {}; rowsD = {};
for k = 1:numel(E)
    e = E(k); t = e.t; v = e.v;
    [aw, ~] = r2_window_amplitude(v, (1:C.STRIDE:(size(v, 1) - C.W + 1))', C.W, C.fs, C.f);
    a_rel = aw/sqrt(3/2);
    for u = 1:2
        o = r2_decision_chain(v, t, r2_opts_nominal(C.opts, 1), 'scale', SCL(u));
        w = o.win; pre = w.armed & w.t_end <= e.tE + 1e-9; post = w.t_end > e.tE + 1e-9;
        dpre = o.decl(o.decl.t_decl <= e.tE + 1e-9, :); dpost = o.decl(o.decl.t_decl > e.tE + 1e-9, :);
        vias = cell(height(o.decl), 1);
        for q = 1:height(o.decl)
            j = o.decl.j(q); run = (j - C.M + 1):j; fl = mean(w.is_line(run));
            if ismember(o.decl.class{q}, C.GR), vias{q} = 'tierra';
            elseif strcmp(o.decl.class{q}, 'A-B-C-G'), vias{q} = 'trifasica por magnitud';
            elseif fl == 1, vias{q} = 'rama degenerada';
            elseif fl == 0, vias{q} = 'puerta de redondez';
            else, vias{q} = 'mixta (degenerada y redondez)';
            end
            ia = (w.i_start(run) - 1)/C.STRIDE + 1;
            rowsD(end+1, :) = {e.kind, e.q, UN{u}, 1e3*(o.decl.t_decl(q) - e.tE), o.decl.class{q}, vias{q}, median(a_rel(ia)), min(a_rel(ia)), median(w.roundness(run)), e.par.p1, e.par.p2}; %#ok<AGROW>
        end
        ipost = o.decl.t_decl > e.tE + 1e-9;
        if isempty(dpost), l1 = NaN; else, l1 = 1e3*(dpost.t_decl(1) - e.tE); end
        rows(end+1, :) = {e.kind, e.q, UN{u}, e.par.p1, e.par.p2, e.par.p3, e.par.desc, ...
            height(dpre), cnt_str(dpre.class), height(dpost), cnt_str(dpost.class), cnt_str(vias(ipost)), l1, ...
            100*mean(~strcmp(w.label(pre), 'A-B-C')), 100*mean(~strcmp(w.label(post), 'A-B-C')), 100*mean(ismember(w.label(post), C.GR)), ...
            100*mean(w.deg_decides(post)), 100*mean(w.is_line(post)), max([w.run(post); 0]), median(w.roundness(pre), 'omitnan'), min(w.roundness(post)), ...
            median(w.pinmax(pre)), max(w.pinmax(post)), min(a_rel((w.i_start(post) - 1)/C.STRIDE + 1))}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'familia','semilla','unidades','p1','p2','p3','parametros','n_decl_previas','clases_previas','n_decl_evento','clases_evento','vias_evento','primera_decl_ms', ...
    'pre_no_sana_pct','post_no_sana_pct','post_tierra_pct','post_rama_degenerada_pct','post_linea_pct','post_racha_max','pre_redondez_med','post_redondez_min','pre_pinmax_med','post_pinmax_max','post_semieje_mayor_rel_min'});
writetable(T, fullfile(C.OUT, 'eventos_seguridad.csv'), 'QuoteStrings', true);
TD = cell2table(cell(0, 11), 'VariableNames', {'familia','semilla','unidades','t_decl_ms','clase','via','semieje_mayor_rel_med','semieje_mayor_rel_min','redondez_med','p1','p2'});
if ~isempty(rowsD), TD = cell2table(rowsD, 'VariableNames', TD.Properties.VariableNames); end
writetable(TD, fullfile(C.OUT, 'eventos_seguridad_declaraciones.csv'), 'QuoteStrings', true);
% resumen por familia y unidades
fam = unique(T.familia, 'stable'); rowsR = {};
for u = 1:2
    for i = 1:numel(fam)
        s = strcmp(T.familia, fam{i}) & strcmp(T.unidades, UN{u}); Ts = T(s, :);
        sd = strcmp(TD.familia, fam{i}) & strcmp(TD.unidades, UN{u}) & TD.t_decl_ms > 0;
        rowsR(end+1, :) = {fam{i}, UN{u}, height(Ts), nnz(Ts.n_decl_previas > 0), nnz(Ts.n_decl_evento > 0), sum(Ts.n_decl_evento), ...
            cnt_str(TD.clase(sd)), cnt_str(TD.via(sd)), nnz(Ts.post_rama_degenerada_pct > 0), max(Ts.post_rama_degenerada_pct), ...
            median(Ts.post_no_sana_pct), max(Ts.post_racha_max), median(Ts.primera_decl_ms, 'omitnan')}; %#ok<AGROW>
        fprintf('%-12s %-7s | eventos con declaracion previa %2d, durante o tras el evento %2d de %2d (%3d decl.: %s; vias: %s) | rama degenerada en %d eventos (max %.2f %%)\n', ...
            fam{i}, UN{u}, rowsR{end, 4}, rowsR{end, 5}, rowsR{end, 3}, rowsR{end, 6}, rowsR{end, 7}, rowsR{end, 8}, rowsR{end, 9}, rowsR{end, 10});
    end
    s = strcmp(T.unidades, UN{u});
    fprintf('TOTAL %-7s: %d de %d eventos con alguna declaracion tras el evento; %d con declaracion previa\n', UN{u}, nnz(T.n_decl_evento(s) > 0), nnz(s), nnz(T.n_decl_previas(s) > 0));
end
TR = cell2table(rowsR, 'VariableNames', {'familia','unidades','n_eventos','eventos_con_decl_previa','eventos_con_decl','n_declaraciones','clases','vias','eventos_con_rama_degenerada','rama_degenerada_pct_max','post_no_sana_pct_mediana','racha_max','primera_decl_ms_mediana'});
writetable(TR, fullfile(C.OUT, 'eventos_seguridad_resumen.csv'), 'QuoteStrings', true);
end

%% =====================================================================
function block_D(C)
fprintf('\n================ BLOQUE D: longitud real de la trayectoria en las faltas del paper ================\n');
info = C.info; Tc = C.Tc; rows = {};
for i = 1:height(Tc)
    d = load(fullfile(info.dir, Tc.file{i}));
    for c = 1:numel(C.CONV)
        [V, t] = extract_converter_data(d, C.CONV{c});
        o = r2_decision_chain(V, t, r2_opts_nominal(C.opts, C.Npk(c)), 'truth', Tc.truth{i}, 't_on', info.t_on, 't_off', info.t_off, 't_range', info.t_range);
        w = o.win; full = w.t_start >= info.t_on & w.t_end <= info.t_off; pre = w.armed & w.t_end <= info.t_on;
        [a, b] = r2_window_amplitude(V, w.i_start, C.W, C.fs, C.f);
        a0 = median(a(pre)); e = o.event;
        rows(end+1, :) = {Tc.file{i}, Tc.truth{i}, Tc.family{i}, C.CONV{c}, e.outcome, e.decl_ms, 100*mean(w.is_line(full)), 100*mean(w.deg_decides(full)), 100*mean(w.deg_bnorm(full)), ...
            median(a(full))/a0, min(a(full))/a0, max(a(full))/a0, median(b(full))/a0, median(w.semi_major(full))/a0, median(w.semi_major(full & ~w.deg_bnorm), 'omitnan')/a0, median(w.roundness(full)), a0}; %#ok<AGROW>
        fprintf('%-18s %-5s %-3s | %-15s %6.2f ms | linea %5.1f %% decide %5.1f %% biv.nulo %5.1f %% | semieje mayor real %.3f (min %.3f, max %.3f) del nominal; menor %.3f; del clasificador %.4g\n', ...
            Tc.file{i}, Tc.truth{i}, C.CONV{c}, e.outcome, e.decl_ms, rows{end, 7}, rows{end, 8}, rows{end, 9}, rows{end, 10}, rows{end, 11}, rows{end, 12}, rows{end, 13}, rows{end, 14});
    end
end
T = cell2table(rows, 'VariableNames', {'archivo','clase','familia','convertidor','resultado','decl_ms','linea_pct','rama_degenerada_decide_pct','biv_nulo_pct', ...
    'semieje_mayor_real_rel_med','semieje_mayor_real_rel_min','semieje_mayor_real_rel_max','semieje_menor_real_rel_med','semieje_mayor_clasificador_rel_med','semieje_mayor_clasificador_sin_biv_nulo_rel_med','redondez_med','semieje_nominal_V'});
writetable(T, fullfile(C.OUT, 'longitud_linea_faltas.csv'), 'QuoteStrings', true);
end
