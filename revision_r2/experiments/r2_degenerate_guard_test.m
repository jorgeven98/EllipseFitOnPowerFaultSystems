function r2_degenerate_guard_test(opts, tag)
% R2_DEGENERATE_GUARD_TEST  Prueba de la decision D5 (T41): condiciones de amplitud.
%
%   r2_degenerate_guard_test()             % cadena de referencia
%   r2_degenerate_guard_test(opts, tag)    % otra variante de la cadena; salidas en T41/tag
%
% Compara, con la misma cadena (r2_decision_chain y las opciones de `opts`),
% cinco variantes de r2_ga_classifier_guard:
%   referencia        sin condiciones (clasificador tal cual)
%   semieje_linea     semieje mayor >= 0.5 del nominal en la rama degenerada
%                     (la comprobacion que describia el texto del paper)
%   tension_minima    bloqueo de cualquier etiqueta si el semieje mayor < 0.1 del nominal
%   arranque_subtension  falta solo si alguna fase baja de 0.9 pu
%   semieje+minima    las dos primeras a la vez
% sobre (1) los 320 eventos de seguridad, (2) los escenarios de T41 que dan
% declaraciones falsas mas un caso minimo de linea de baja amplitud en el
% plano de Kirchhoff, y (3) las 24 faltas del paper (coste en fiabilidad).
% Salidas en resultados/T41[/tag]: guard_eventos.csv, guard_escenarios.csv,
% guard_faltas.csv, guard_resumen.csv, log_guard.txt.
if nargin < 1 || isempty(opts), opts = struct(); end
if nargin < 2 || isempty(tag), tag = ''; end
here = fileparts(mfilename('fullpath')); RAIZ = fullfile(here, '..');
addpath(here); addpath(fullfile(RAIZ, 'src', 'fault_classifiers')); addpath(fullfile(RAIZ, 'src', 'utils'));
OUT = fullfile(RAIZ, '_revision_TPWRD', 'ronda2', 'resultados', 'T41', tag);
if ~exist(OUT, 'dir'), mkdir(OUT); end
logf = fullfile(OUT, 'log_guard.txt'); if exist(logf, 'file'), delete(logf); end
diary(logf); cobj = onCleanup(@() diary('off'));
fprintf('R2_DEGENERATE_GUARD_TEST  %s  MATLAB %s\nopciones de la cadena distintas de las de referencia:\n', datestr(now, 31), version); disp(opts); %#ok<TNOW1,DATST>
warning('off', 'all');
fs = 4000; f = 50; W = 20; STRIDE = 5; GR = {'A-G','B-G','C-G','A-B-G','B-C-G','C-A-G'};
GU = {'referencia', struct()
      'semieje_linea', struct('a_min_line', 0.5)
      'tension_minima', struct('v_min', 0.1)
      'arranque_subtension', struct('uv_start', 0.9)
      'semieje+minima', struct('a_min_line', 0.5, 'v_min', 0.1)};
    function o = runchain(V, t, vnom, g, varargin)
        % clasificador base segun opts (la cadena ignora 'classifier' si se cambian umbrales o variantes internas)
        [base, copts] = base_from_opts(opts, vnom/sqrt(3/2));
        G = GU{g, 2}; G.base = base; G.v_nominal = vnom;
        o = r2_decision_chain(V, t, copts, 'classifier', @(Vw, fs_, f_) r2_ga_classifier_guard(Vw, fs_, f_, G), varargin{:});
    end
[Tc, info] = r2_sim_cases(); CONV = info.converters;
d0 = load(fullfile(info.dir, 'DataGridRL_AB.mat'));
for c = 1:2
    [V, t] = extract_converter_data(d0, CONV{c}); Npk(c) = max(abs(V(t >= 2.4 & t < 2.8, :)), [], 'all'); %#ok<AGROW>
    seg = t >= 0.4 & t < 2.95; Vh{c} = V(seg, :); th{c} = t(seg); %#ok<AGROW>
end

%% (1) eventos de seguridad
fprintf('\n== (1) Suite de 320 eventos de seguridad (pu)\n');
E = r2_security_events_gen(40); rows = {};
for k = 1:numel(E)
    e = E(k);
    for g = 1:size(GU, 1)
        o = runchain(e.v, e.t, sqrt(3/2), g);
        dp = o.decl(o.decl.t_decl > e.tE, :); dpre = o.decl(o.decl.t_decl <= e.tE, :);
        post = o.win.t_end > e.tE;
        rows(end+1, :) = {e.kind, e.q, GU{g, 1}, e.par.p1, e.par.p2, height(dpre), height(dp), cstr(dp.class), 100*mean(strcmp(o.win.label_raw(post), 'Unclassified'))}; %#ok<AGROW>
    end
end
TE = cell2table(rows, 'VariableNames', {'familia','semilla','variante','p1','p2','n_decl_previas','n_decl_evento','clases','ventanas_bloqueadas_pct'});
writetable(TE, fullfile(OUT, 'guard_eventos.csv'), 'QuoteStrings', true);
fam = unique(TE.familia, 'stable'); rowsR = {};
for g = 1:size(GU, 1)
    fprintf('%-20s:', GU{g, 1});
    for i = 1:numel(fam)
        s = strcmp(TE.familia, fam{i}) & strcmp(TE.variante, GU{g, 1});
        fprintf(' %s %d/%d', fam{i}, nnz(TE.n_decl_evento(s) > 0), nnz(s));
        rowsR(end+1, :) = {'eventos de seguridad', fam{i}, GU{g, 1}, nnz(s), nnz(TE.n_decl_evento(s) > 0), sum(TE.n_decl_evento(s)), cstr(TE.clases(s & TE.n_decl_evento > 0))}; %#ok<AGROW>
    end
    fprintf('\n');
end

%% (2) escenarios de T41 con declaraciones falsas y caso minimo para D5
fprintf('\n== (2) Escenarios (voltios, como se publico; 5 realizaciones de 2.3 s)\n');
SC = {'linea sin tension, ruido 0.01 del pico', 'dead', struct('sig', 1e-2, 'lsb', 0), 5
      'linea sin tension, cero exacto', 'dead', struct('sig', 0, 'lsb', 0), 1
      'linea sin tension, ruido de 0.3 LSB cuantizado (12 bits)', 'dead', struct('sig', 0.3/2048, 'lsb', 1/2048), 5
      'linea sin tension, inducida una fase (1,0,0) de 5 % + ruido 0.1 %', 'induced', struct('u', [1 0 0], 'eps', 0.05, 'sig', 1e-3), 5
      'linea sin tension, inducida entre fases (1,-1,0) de 5 % + ruido 0.1 %', 'induced', struct('u', [1 -1 0], 'eps', 0.05, 'sig', 1e-3), 5
      'CASO MINIMO D5: inducida entre fases (1,-1,0) de 5 % + ruido 0.001 %', 'induced', struct('u', [1 -1 0], 'eps', 0.05, 'sig', 1e-5), 5
      'CASO MINIMO D5: inducida entre fases (1,-1,0) de 20 % + ruido 0.01 %', 'induced', struct('u', [1 -1 0], 'eps', 0.20, 'sig', 1e-4), 5
      'hueco trifasico equilibrado, tension residual 0.05 pu, ruido 60 dB', 'dip3', struct('r', 0.05, 'snr', 60), 5
      'perdida de tension en A (ruido residual 0.002)', 'lop', struct('ph', 1, 'sig', 0.002), 5
      'perdida de tension en B y C (ruido residual 0.002)', 'lop', struct('ph', [2 3], 'sig', 0.002), 5
      'control positivo: A-B franca sintetica', 'bolted', struct('perm', [1 2 3]), 1
      'control positivo: B-C franca sintetica', 'bolted', struct('perm', [2 3 1]), 1
      'control positivo: C-A franca sintetica', 'bolted', struct('perm', [3 1 2]), 1};
rows = {};
for c = 1:2
    V0 = Vh{c}; tl = th{c}; on = (tl - tl(1)) >= 0.65; t_ev = tl(find(on, 1)); n = nnz(on);
    for i = 1:size(SC, 1)
        p = SC{i, 3};
        for g = 1:size(GU, 1)
            rs = RandStream('mt19937ar', 'Seed', 4700 + 100*c + i);       % mismas realizaciones para todas las variantes
            nd = 0; nrec = 0; cls = {}; nb = 0; nw = 0; ndeg = 0; vias_deg = 0;
            for rep = 1:SC{i, 4}
                V = V0;
                switch SC{i, 2}
                    case 'dead', Xn = p.sig*Npk(c)*randn(rs, n, 3); if p.lsb > 0, Xn = round(Xn/(p.lsb*Npk(c)))*(p.lsb*Npk(c)); end, V(on, :) = Xn;
                    case 'induced', ph = 2*pi*rand(rs); V(on, :) = p.eps*Npk(c)*cos(2*pi*f*tl(on) + ph)*p.u + p.sig*Npk(c)*randn(rs, n, 3);
                    case 'dip3', V(on, :) = p.r*V0(on, :); V = V + (Npk(c)/sqrt(2))*10^(-p.snr/20)*randn(rs, size(V0));
                    case 'lop', V(on, p.ph) = p.sig*Npk(c)*randn(rs, n, numel(p.ph));
                    case 'bolted', V(on, p.perm(1)) = -V0(on, p.perm(3))/2; V(on, p.perm(2)) = -V0(on, p.perm(3))/2;
                end
                o = runchain(V, tl, sqrt(3/2)*Npk(c), g); w = o.win; m = w.t_start >= t_ev - 1e-9;
                dm = o.decl(o.decl.t_decl > t_ev, :);
                nd = nd + height(dm); nrec = nrec + (height(dm) > 0); cls = [cls; dm.class]; %#ok<AGROW>
                nb = nb + nnz(strcmp(w.label_raw(m), 'Unclassified')); nw = nw + nnz(m); ndeg = ndeg + nnz(w.deg_decides(m));
                for q = 1:height(dm)
                    if ~ismember(dm.class{q}, GR) && all(w.is_line((dm.j(q) - 9):dm.j(q))), vias_deg = vias_deg + 1; end
                end
            end
            rows(end+1, :) = {SC{i, 1}, CONV{c}, GU{g, 1}, SC{i, 4}, nw*STRIDE/fs, nd, nrec, cstr(cls), vias_deg, 100*nb/nw, 100*ndeg/nw}; %#ok<AGROW>
            if g == 1 || nd ~= rows{end - g + 1, 6}
                fprintf('%-70s %-3s %-20s | decl %3d en %d/%d reg (%s) | por rama degenerada %d | bloqueadas %5.1f %%\n', SC{i, 1}, CONV{c}, GU{g, 1}, nd, nrec, SC{i, 4}, cstr(cls), vias_deg, 100*nb/nw);
            end
        end
    end
end
TS = cell2table(rows, 'VariableNames', {'escenario','convertidor','variante','n_registros','exposicion_s','n_declaraciones','registros_con_declaracion','clases','n_decl_por_rama_degenerada','ventanas_bloqueadas_pct','rama_degenerada_decide_pct'});
writetable(TS, fullfile(OUT, 'guard_escenarios.csv'), 'QuoteStrings', true);

%% (3) coste en fiabilidad: las 24 faltas del paper
fprintf('\n== (3) Las 24 faltas del paper\n'); rows = {};
for i = 1:height(Tc)
    d = load(fullfile(info.dir, Tc.file{i}));
    for c = 1:2
        [V, t] = extract_converter_data(d, CONV{c});
        for g = 1:size(GU, 1)
            o = runchain(V, t, sqrt(3/2)*Npk(c), g, 'truth', Tc.truth{i}, 't_on', info.t_on, 't_off', info.t_off, 't_range', info.t_range); e = o.event;
            full = o.win.t_start >= info.t_on & o.win.t_end <= info.t_off;
            rows(end+1, :) = {Tc.file{i}, Tc.truth{i}, Tc.family{i}, CONV{c}, GU{g, 1}, e.outcome, e.decl_class, e.decl_ms, e.n_decl_pre, 100*mean(strcmp(o.win.label_raw(full), 'Unclassified'))}; %#ok<AGROW>
        end
    end
end
TF = cell2table(rows, 'VariableNames', {'archivo','clase','familia','convertidor','variante','resultado','clase_declarada','decl_ms','n_decl_prefalta','ventanas_falta_bloqueadas_pct'});
writetable(TF, fullfile(OUT, 'guard_faltas.csv'), 'QuoteStrings', true);
for g = 1:size(GU, 1)
    s = strcmp(TF.variante, GU{g, 1}); ok = strcmp(TF.resultado, 'correcta');
    fprintf('%-20s: %2d/24 correctas (Pex %2d/12, DG1 %2d/12); latencia mediana %.2f ms, max %.2f ms; no declaradas: %s\n', GU{g, 1}, nnz(s & ok), nnz(s & ok & strcmp(TF.convertidor, 'Pex')), nnz(s & ok & strcmp(TF.convertidor, 'DG1')), ...
        median(TF.decl_ms(s & ok)), max(TF.decl_ms(s & ok)), strjoin(unique(strcat(TF.archivo(s & ~ok), '/', TF.convertidor(s & ~ok)))', ' '));
    rowsR(end+1, :) = {'faltas del paper', 'las 24', GU{g, 1}, nnz(s), nnz(s & ok), NaN, sprintf('latencia mediana %.2f ms, max %.2f ms', median(TF.decl_ms(s & ok)), max(TF.decl_ms(s & ok)))}; %#ok<AGROW>
end
for i = 1:size(SC, 1)
    for g = 1:size(GU, 1)
        s = strcmp(TS.escenario, SC{i, 1}) & strcmp(TS.variante, GU{g, 1});
        rowsR(end+1, :) = {'escenarios', SC{i, 1}, GU{g, 1}, sum(TS.n_registros(s)), sum(TS.registros_con_declaracion(s)), sum(TS.n_declaraciones(s)), cstr(TS.clases(s & TS.n_declaraciones > 0))}; %#ok<AGROW>
    end
end
TR = cell2table(rowsR, 'VariableNames', {'bloque','caso','variante','n','n_con_declaracion_o_correctas','n_declaraciones','detalle'});
writetable(TR, fullfile(OUT, 'guard_resumen.csv'), 'QuoteStrings', true);
fprintf('\nFIN r2_degenerate_guard_test.\n');
end

function [base, copts] = base_from_opts(opts, vpk)
% Clasificador base equivalente al que usaria la cadena con `opts`, y las
% opciones restantes. La cadena ignora 'classifier' cuando se cambia
% kirchhoff_thr, cls_params, v0_mode o gnd_mode (variantes de T52), asi que
% aqui se mueven al clasificador base (misma construccion que scan_ de la
% cadena, con scale = 1). vpk: pico de fase nominal de la senal.
opts = r2_opts_nominal(opts, vpk);
copts = opts; P = struct(); usep = false;
if isfield(opts, 'v0_mode') && ~strcmp(opts.v0_mode, 'samples3')
    P.v0_mode = opts.v0_mode; usep = true;
    if isfield(opts, 'v0_nominal') && ~isempty(opts.v0_nominal), P.v0_nominal = opts.v0_nominal; end
end
if isfield(opts, 'gnd_mode') && ~strcmp(opts.gnd_mode, '2s'), P.gnd_mode = opts.gnd_mode; usep = true; end
for fn = {'v0_mode', 'v0_nominal', 'gnd_mode'}
    if isfield(copts, fn{1}), copts = rmfield(copts, fn{1}); end
end
if isfield(opts, 'kirchhoff_thr'), P.kirchhoff_thr = opts.kirchhoff_thr; end
if isfield(opts, 'cls_params') && ~isempty(opts.cls_params), fq = fieldnames(opts.cls_params); for q = 1:numel(fq), P.(fq{q}) = opts.cls_params.(fq{q}); end, usep = true; copts.cls_params = []; end
if isfield(opts, 'kirchhoff_thr') && opts.kirchhoff_thr ~= 0.01, usep = true; end
if isfield(copts, 'kirchhoff_thr'), copts.kirchhoff_thr = 0.01; end
if usep
    base = @(Vw, fs, f) r2_ga_classifier_param(Vw, fs, f, P);
elseif isfield(opts, 'classifier') && ~isempty(opts.classifier)
    base = opts.classifier; if ischar(base) || isstring(base), base = str2func(char(base)); end
else
    base = @ga_fault_classifier_pure;
end
if isfield(copts, 'classifier'), copts = rmfield(copts, 'classifier'); end
end

function s = cstr(c)
if isempty(c), s = ''; return, end
c = c(~cellfun(@isempty, c)); if isempty(c), s = ''; return, end
[u, ~, j] = unique(c); n = accumarray(j, 1);
s = strjoin(arrayfun(@(k) sprintf('%s x%d', u{k}, n(k)), 1:numel(u), 'UniformOutput', false), '; ');
end
