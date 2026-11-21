function r2_gate_sensitivity(opts, tag, blocks)
% R2_GATE_SENSITIVITY  Sensibilidad de la puerta adaptativa de redondez a Delta y K (T61, R2.2).
%
%   r2_gate_sensitivity()                     % cadena de referencia, todos los bloques
%   r2_gate_sensitivity(opts)                 % opts = struct con opciones de r2_decision_chain
%   r2_gate_sensitivity(opts, 'G2')           % salidas en resultados/T61/G2
%   r2_gate_sensitivity(opts, tag, {'S','F'}) % solo algunos bloques
%
% Toda decision sale de r2_decision_chain (nivel ambiente CAUSAL) con las
% opciones de `opts`; aqui solo se barren Delta y K. Cada senal se clasifica
% una vez y su escaneo (out.scan) se reutiliza para todas las parejas.
% Para repetir con otra variante del clasificador:
%   r2_gate_sensitivity(struct('kirchhoff_thr', 0.02), 'kthr002')
%
% Rejilla: Delta en {0.005, 0.01, 0.02, 0.03} x K en {20, 40, 80, 160} ventanas.
% Bloques
%   S  Seguridad en sano con desequilibrio ambiente V2 de 0 a 2 % (convenio de
%      d1_roundness_gate_design.m), con barrido de la fase de V2 y con ruido.
%   F  Faltas del paper (12 archivos x 2 convertidores) con V2 ambiente de
%      0, 0.5, 1 y 2 %: declaracion, clase y latencia.
%   D  Suelo de detectabilidad: faltas A-B sinteticas de severidad k.
%   I  Inicializacion, congelacion y rearme: falta antes de completar K
%      ventanas, escalon legitimo del desequilibrio ambiente, falta de
%      evolucion lenta, segunda falta tras el despeje y recuperacion real.
%   L  Registros de laboratorio (Fault_emulator_Outkontaktor): ambiente real.
%   G  Dispersion de la redondez sana y Delta minimo seguro por traza (guia).
% Salidas en resultados/T61[/tag]: sens_delta_K.csv y los CSV de detalle.
if nargin < 1 || isempty(opts), opts = struct(); end
if nargin < 2 || isempty(tag), tag = ''; end
if nargin < 3 || isempty(blocks), blocks = {'S','F','D','I','L','G'}; end
here = fileparts(mfilename('fullpath')); RAIZ = fullfile(here, '..');
addpath(here); addpath(fullfile(RAIZ, 'src', 'fault_classifiers')); addpath(fullfile(RAIZ, 'src', 'utils'));
OUT = fullfile(RAIZ, '_revision_TPWRD', 'ronda2', 'resultados', 'T61', tag);
if ~exist(OUT, 'dir'), mkdir(OUT); end
logf = fullfile(OUT, sprintf('log_%s.txt', strjoin(blocks, '')));
if exist(logf, 'file'), delete(logf); end
diary(logf); cobj = onCleanup(@() diary('off'));
fprintf('R2_GATE_SENSITIVITY  %s  MATLAB %s\nbloques: %s   tag: "%s"\nopciones de la cadena distintas de las de referencia:\n', ...
    datestr(now, 31), version, strjoin(blocks, ','), tag); %#ok<TNOW1,DATST>
disp(opts); warning('off', 'all');

C.opts = opts; C.OUT = OUT; C.RAIZ = RAIZ;
C.DELTAS = [0.005 0.01 0.02 0.03]; C.KS = [20 40 80 160];
C.DFINE = [0.001:0.001:0.010, 0.0125:0.0025:0.05];       % para el Delta minimo seguro
C.fs = 4000; C.f = 50; C.W = 20; C.STRIDE = 5;
d0 = r2_decision_chain('defaults'); C.M = d0.M; if isfield(opts, 'M'), C.M = opts.M; end
[C.Tc, C.info] = r2_sim_cases(); C.CONV = C.info.converters;
C.GR = {'A-G','B-G','C-G','A-B-G','B-C-G','C-A-G'};
fprintf('Rejilla: Delta = %s; K = %s ventanas.\n', mat2str(C.DELTAS), mat2str(C.KS));
for K = C.KS
    fprintf('  K = %3d ventanas abarcan %d muestras = %.2f ms = %.2f ciclos a 4 kHz (paso %d, W = %d)\n', K, (K-1)*C.STRIDE + C.W, 1e3*((K-1)*C.STRIDE + C.W)/C.fs, ((K-1)*C.STRIDE + C.W)/C.fs*C.f, C.STRIDE, C.W);
end
G = struct();
if ismember('S', blocks), G.S = block_S(C); end
if ismember('F', blocks), G.F = block_F(C); end
if ismember('D', blocks), G.D = block_D(C); end
if ismember('L', blocks), G.L = block_L(C); end
if ismember('I', blocks), block_I(C); end
if all(ismember({'S','L'}, blocks)) && ismember('G', blocks), block_G(C, G); end
if all(ismember({'S','F','D','L'}, blocks)), grid_summary(C, G); end
fprintf('\nFIN r2_gate_sensitivity (%s).\n', datestr(now, 31)); %#ok<TNOW1,DATST>
end

%% ===================================================================== utilidades
function v2 = v2wave(k2, npk, f, t, phi)
% Secuencia inversa sostenida de amplitud k2*npk (convenio de d1_roundness_gate_design.m; phi = 0 alli)
v2 = k2*npk*[cos(2*pi*f*t + phi), cos(2*pi*f*t + phi + 2*pi/3), cos(2*pi*f*t + phi - 2*pi/3)];
end
function S2 = subscan(S, keep)
% Subconjunto de ventanas de un escaneo (equivale a empezar a procesar mas tarde)
m = numel(S.label_raw); S2 = S; fn = fieldnames(S);
for k = 1:numel(fn)
    if size(S.(fn{k}), 1) == m && ~ischar(S.(fn{k})), S2.(fn{k}) = S.(fn{k})(keep, :); end
end
end
function s = cnt_str(c)
if isempty(c), s = ''; return, end
[u, ~, j] = unique(c); n = accumarray(j, 1);
s = strjoin(arrayfun(@(k) sprintf('%s x%d', u{k}, n(k)), 1:numel(u), 'UniformOutput', false), '; ');
end
function v = pct_(x, p)
x = sort(x(isfinite(x))); n = numel(x);
if n == 0, v = NaN; elseif n == 1, v = x; else, v = interp1(((1:n) - 0.5)/n, x, min(max(p/100, 0.5/n), (n - 0.5)/n)); end
end

%% ===================================================================== S: sano con V2 ambiente
function R = block_S(C)
fprintf('\n================ BLOQUE S: seguridad en sano con desequilibrio ambiente ================\n');
base = load(fullfile(C.info.dir, 'DataGridRL_AB.mat'));
V2S = [0 0.005 0.01 0.015 0.02]; PHIS = [0 90 180 270]; SNRS = [Inf 50 40];
rows = {}; R.sc = {}; R.meta = {};
for c = 1:numel(C.CONV)
    [V, t] = extract_converter_data(base, C.CONV{c});
    npk = max(abs(V(t >= 2.4 & t < 2.8, :)), [], 'all');
    seg = t >= 0.4 & t < 2.95; Vh = V(seg, :); th = t(seg);
    for k2 = V2S
        for phi = PHIS
            if k2 == 0 && phi ~= 0, continue, end
            for snr = SNRS
                if isfinite(snr) && phi ~= 0, continue, end
                Vs = Vh + v2wave(k2, npk, C.f, th, deg2rad(phi));
                if isfinite(snr)
                    rs = RandStream('mt19937ar', 'Seed', 6100 + 10*c + round(1000*k2) + snr);
                    Vs = Vs + (npk/sqrt(2))*10^(-snr/20)*randn(rs, size(Vs));
                end
                o = r2_decision_chain(Vs, th, r2_opts_nominal(C.opts, npk));
                name = sprintf('%s V2 = %.1f %% fase %d ruido %g dB', C.CONV{c}, 100*k2, phi, snr);
                R.sc{end+1} = o.scan; R.meta(end+1, :) = {name, C.CONV{c}, 100*k2, phi, snr, 'simulacion'};
                for D = C.DELTAS
                    for K = C.KS
                        a = r2_decision_chain(o.scan, [], C.opts, 'Delta', D, 'K', K);
                        w = a.win; m = w.armed;
                        rows(end+1, :) = {name, C.CONV{c}, 100*k2, phi, snr, D, K, nnz(m), nnz(m)*C.STRIDE/C.fs, height(a.decl), cnt_str(a.decl.class), 100*mean(w.fire(m)), ...
                            100*mean(~strcmp(w.label(m), 'A-B-C')), 100*mean(~strcmp(w.label_raw(m), 'A-B-C')), 100*mean(w.fed(m)), median(w.r_amb(m), 'omitnan'), median(w.thr(m), 'omitnan'), max([w.run(m); 0])}; %#ok<AGROW>
                    end
                end
            end
        end
    end
end
T = cell2table(rows, 'VariableNames', {'traza','convertidor','v2_pct','fase_v2_grados','snr_db','delta','K','n_ventanas','exposicion_s','n_declaraciones','clases','ventanas_aseguradas_pct', ...
    'etiqueta_no_sana_pct','etiqueta_cruda_no_sana_pct','ventanas_que_alimentan_pct','r_amb_med','umbral_med','racha_max'});
writetable(T, fullfile(C.OUT, 'sano_detalle.csv'), 'QuoteStrings', true);
R.T = T;
% resumen por pantalla: trazas publicadas (fase 0, sin ruido)
fprintf('Trazas con alguna declaracion falsa, de las %d sin ruido (2 convertidores, V2 de 0 a 2 %%, cuatro fases):\n', nnz(isinf(cell2mat(R.meta(:, 5)))));
for D = C.DELTAS
    fprintf('  Delta = %.3f:', D);
    for K = C.KS
        s = T.delta == D & T.K == K & isinf(T.snr_db);
        fprintf('  K=%3d: %2d trazas (%3d decl.)', K, nnz(T.n_declaraciones(s) > 0), sum(T.n_declaraciones(s)));
    end
    fprintf('\n');
end
fprintf('Lo mismo con ruido de 50 y 40 dB (fase 0):\n');
for D = C.DELTAS
    fprintf('  Delta = %.3f:', D);
    for K = C.KS
        s = T.delta == D & T.K == K & isfinite(T.snr_db);
        fprintf('  K=%3d: %2d de %2d trazas (%3d decl.)', K, nnz(T.n_declaraciones(s) > 0), nnz(s), sum(T.n_declaraciones(s)));
    end
    fprintf('\n');
end
end

%% ===================================================================== F: faltas del paper
function R = block_F(C)
fprintf('\n================ BLOQUE F: faltas del paper con desequilibrio ambiente ================\n');
info = C.info; Tc = C.Tc; rows = {}; rowsM = {};
V2S = [0 0.005 0.01 0.02];
for i = 1:height(Tc)
    d = load(fullfile(info.dir, Tc.file{i})); isLL = ismember(Tc.family{i}, {'LL', 'LL_RL'});
    for c = 1:numel(C.CONV)
        [V, t] = extract_converter_data(d, C.CONV{c});
        npk = max(abs(V(t >= 2.4 & t < 2.8, :)), [], 'all');
        for k2 = V2S
            if k2 > 0 && isLL, PHIS = [0 90 180 270]; else, PHIS = 0; end
            for phi = PHIS
                Vs = V + v2wave(k2, npk, C.f, t, deg2rad(phi));
                ev = {'truth', Tc.truth{i}, 't_on', info.t_on, 't_off', info.t_off};
                o = r2_decision_chain(Vs, t, r2_opts_nominal(C.opts, npk), ev{:}, 't_range', info.t_range);
                % Delta maximo con el que la falta se sigue declarando bien (y con todas las Delta menores)
                dmax = nan(1, numel(C.KS));
                for ik = 1:numel(C.KS)
                    for D = C.DFINE
                        a = r2_decision_chain(o.scan, [], C.opts, 'Delta', D, 'K', C.KS(ik), ev{:});
                        if strcmp(a.event.outcome, 'correcta'), dmax(ik) = D; else, break, end
                    end
                end
                rowsM(end+1, :) = [{Tc.file{i}, Tc.truth{i}, Tc.family{i}, C.CONV{c}, 100*k2, phi}, num2cell(dmax)]; %#ok<AGROW>
                for D = C.DELTAS
                    for K = C.KS
                        a = r2_decision_chain(o.scan, [], C.opts, 'Delta', D, 'K', K, ev{:}); e = a.event;
                        w = a.win; flt = w.t_start >= info.t_on & w.t_end <= info.t_off;
                        rows(end+1, :) = {Tc.file{i}, Tc.truth{i}, Tc.family{i}, C.CONV{c}, 100*k2, phi, D, K, e.outcome, e.decl_class, e.decl_ms, e.first_label_ms, e.max_run_correct, e.retention_pct, ...
                            e.n_decl_pre, e.pre_classes, e.n_decl_wrong, e.wrong_classes, e.n_decl_after, e.after_classes, e.full_roundness_med, e.full_thr_med, 100*mean(w.fed(flt))}; %#ok<AGROW>
                    end
                end
            end
        end
    end
end
T = cell2table(rows, 'VariableNames', {'archivo','clase','familia','convertidor','v2_pct','fase_v2_grados','delta','K','resultado','clase_declarada','decl_ms','primera_etq_ms','racha_max','retencion_pct', ...
    'n_decl_prefalta','clases_prefalta','n_decl_erroneas','clases_erroneas','n_decl_posfalta','clases_posfalta','redondez_falta_med','umbral_falta_med','ventanas_falta_que_alimentan_pct'});
writetable(T, fullfile(C.OUT, 'faltas_detalle.csv'), 'QuoteStrings', true);
R.T = T;
kn = arrayfun(@(K) sprintf('delta_max_fiable_K%d', K), C.KS, 'UniformOutput', false);
TM = cell2table(rowsM, 'VariableNames', [{'archivo','clase','familia','convertidor','v2_pct','fase_v2_grados'}, kn]);
writetable(TM, fullfile(C.OUT, 'faltas_delta_max.csv'), 'QuoteStrings', true);
pubM = ismember(TM.familia, {'LL_RL'}) & ismember(TM.v2_pct, [0 1]) & TM.fase_v2_grados == 0; nomM = TM.v2_pct == 0;
fprintf('Delta maximo fiable (rejilla fina hasta %.3f; NaN = no se declara ni con el menor) [minimo sobre las 12 publicadas | minimo sobre las 24 nominales]:\n', C.DFINE(end));
for ik = 1:numel(C.KS)
    x = TM.(kn{ik}); fprintf('  K = %3d: %.4f | %.4f\n', C.KS(ik), min(x(pubM), [], 'includenan'), min(x(nomM), [], 'includenan'));
end
ok = strcmp(T.resultado, 'correcta');
pub = ismember(T.familia, {'LL_RL'}) & ismember(T.v2_pct, [0 1]) & T.fase_v2_grados == 0;
nom = T.v2_pct == 0;
fprintf('Faltas declaradas con la clase correcta [publicadas 12 = 3 RL x 2 conv x V2 0 y 1 %% | 24 nominales | todas con V2 > 0]:\n');
for D = C.DELTAS
    for K = C.KS
        s = T.delta == D & T.K == K;
        fprintf('  Delta = %.3f K = %3d: %2d/%2d (lat. %5.2f a %5.2f ms) | %2d/%2d (mediana %5.2f, max %5.2f ms) | %3d/%3d | decl. en prefalta %d\n', D, K, ...
            nnz(ok & s & pub), nnz(s & pub), min(T.decl_ms(ok & s & pub)), max(T.decl_ms(ok & s & pub)), nnz(ok & s & nom), nnz(s & nom), median(T.decl_ms(ok & s & nom)), max(T.decl_ms(ok & s & nom)), ...
            nnz(ok & s & ~nom), nnz(s & ~nom), sum(T.n_decl_prefalta(s)));
    end
end
end

%% ===================================================================== D: suelo de detectabilidad
function R = block_D(C)
fprintf('\n================ BLOQUE D: suelo de detectabilidad (faltas A-B sinteticas) ================\n');
fs = C.fs; f = C.f; T0S = 0.5;
tsyn = (0:4399)'/fs;                                   % 1.1 s; inicio en 0.5 s (como d1)
vbase = [cos(2*pi*f*tsyn), cos(2*pi*f*tsyn - 2*pi/3), cos(2*pi*f*tsyn + 2*pi/3)];
KSEV = [0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.05 0.06 0.08 0.10 0.15 0.20 0.30];
V2S = [0 0.005 0.01 0.02]; rows = {};
a = exp(2j*pi/3);
for k2 = V2S
    if k2 == 0, PHIS = 0; else, PHIS = 0:45:315; end
    for phi = PHIS
        for k = KSEV
            V = vbase; dv = (k/2)*(V(:, 1) - V(:, 2)); msk = tsyn >= T0S;
            V(msk, 1) = V(msk, 1) - dv(msk); V(msk, 2) = V(msk, 2) + dv(msk);
            V = V + v2wave(k2, 1, f, tsyn, deg2rad(phi));
            % desequilibrio V2/V1 antes y durante la falta (fasores por minimos cuadrados, un ciclo)
            u = zeros(1, 2); segs = {tsyn >= 0.3 & tsyn < 0.32, tsyn >= 0.8 & tsyn < 0.82};
            for q = 1:2
                X = [cos(2*pi*f*tsyn(segs{q})), -sin(2*pi*f*tsyn(segs{q}))]; cf = X\V(segs{q}, :); ph = cf(1, :) + 1j*cf(2, :);
                u(q) = abs((ph(1) + a^2*ph(2) + a*ph(3))/(ph(1) + a*ph(2) + a^2*ph(3)));
            end
            ev = {'truth', 'A-B', 't_on', T0S};
            o = r2_decision_chain(V, tsyn, r2_opts_nominal(C.opts, 1), ev{:});
            for D = C.DELTAS
                for K = C.KS
                    b = r2_decision_chain(o.scan, [], C.opts, 'Delta', D, 'K', K, ev{:}); e = b.event;
                    w = b.win; flt = w.t_start >= T0S;
                    rows(end+1, :) = {k, 100*k2, phi, D, K, 100*u(1), 100*u(2), e.outcome, e.decl_class, e.decl_ms, e.retention_pct, e.n_decl_pre, e.full_roundness_med, e.full_thr_med, 100*mean(w.fed(flt)), w.r_amb(end)}; %#ok<AGROW>
                end
            end
        end
    end
end
T = cell2table(rows, 'VariableNames', {'k_severidad','v2_ambiente_pct','fase_v2_grados','delta','K','v2v1_prefalta_pct','v2v1_falta_pct','resultado','clase_declarada','decl_ms','retencion_pct','n_decl_prefalta', ...
    'redondez_falta_med','umbral_falta_med','ventanas_falta_que_alimentan_pct','r_amb_final'});
writetable(T, fullfile(C.OUT, 'suelo_detalle.csv'), 'QuoteStrings', true);
% suelo: menor k a partir del cual TODAS las severidades mayores se declaran bien
rowsF = {};
for D = C.DELTAS
    for K = C.KS
        for k2 = V2S
            if k2 == 0, PHIS = 0; else, PHIS = 0:45:315; end
            kf = nan(size(PHIS));
            for ip = 1:numel(PHIS)
                s = T.delta == D & T.K == K & T.v2_ambiente_pct == 100*k2 & T.fase_v2_grados == PHIS(ip);
                Ts = sortrows(T(s, :), 'k_severidad'); okk = strcmp(Ts.resultado, 'correcta');
                lastbad = find(~okk, 1, 'last');
                if isempty(lastbad), kf(ip) = Ts.k_severidad(1); elseif lastbad < height(Ts), kf(ip) = Ts.k_severidad(lastbad + 1); end
            end
            [kw, iw] = max(kf); if any(isnan(kf)), kw = NaN; iw = find(isnan(kf), 1); end
            rowsF(end+1, :) = {D, K, 100*k2, kf(1), min(kf), kw, PHIS(iw), 100*kf(1)/2, 100*kw/2}; %#ok<AGROW>
        end
    end
end
TF = cell2table(rowsF, 'VariableNames', {'delta','K','v2_ambiente_pct','k_suelo_fase0','k_suelo_mejor_fase','k_suelo_peor_fase','fase_peor_grados','v2_falta_suelo_fase0_pct','v2_falta_suelo_peor_pct'});
writetable(TF, fullfile(C.OUT, 'suelo_resumen.csv'), 'QuoteStrings', true);
R.T = T; R.TF = TF;
fprintf('Suelo de severidad k (la falta sintetica aporta V2/V1 = k/2 aprox.) [fase 0 | peor fase]; NaN = no se declara ni con k = 0.30\n');
for D = C.DELTAS
    for K = C.KS
        s = TF.delta == D & TF.K == K;
        fprintf('  Delta = %.3f K = %3d:', D, K);
        for k2 = V2S, q = s & TF.v2_ambiente_pct == 100*k2; fprintf('  V2 %.1f %%: %.3f | %.3f', 100*k2, TF.k_suelo_fase0(q), TF.k_suelo_peor_fase(q)); end
        fprintf('\n');
    end
end
end

%% ===================================================================== L: laboratorio
function R = block_L(C)
fprintf('\n================ BLOQUE L: registros de laboratorio (ambiente real) ================\n');
labdir = fullfile(C.RAIZ, 'data', 'lab', 'Fault_emulator_Outkontaktor'); files = dir(fullfile(labdir, '*.mat'));
fsL = 1e4; WL = 50; STL = 12; CH = 59:61; PERM = [3 2 1];            % como lab_latency_r3m1.m y T03
emap = containers.Map({'ab','bc','ac','abc'}, {'A-B','B-C','C-A','A-B-C-G'});
for K = C.KS, fprintf('  K = %3d ventanas abarcan %.1f ms = %.2f ciclos a 10 kHz (paso %d, W = %d)\n', K, 1e3*((K-1)*STL + WL)/fsL, ((K-1)*STL + WL)/fsL*C.f, STL, WL); end
rows = {}; R.sc = {}; R.meta = {}; rowsA = {};
for k = 1:numel(files)
    fn = files(k).name; pref = regexp(fn, '^[a-z]+', 'match', 'once'); truth = emap(pref);
    dl = load(fullfile(labdir, fn)); Rr = dl.results; t = Rr(1, :)'; t = t - t(1); V = Rr(CH, :)'; V = V(:, PERM);
    [ki, ke] = r2_lab_inception(V, t, fsL, C.f);
    t_on = t(ki); t_off = t(ke);
    ev = {'truth', truth, 't_on', t_on, 't_off', t_off};
    o = r2_decision_chain(V, t, r2_opts_nominal(C.opts, 1), ev{:}, 'fs', fsL, 'W', WL, 'stride', STL);
    w0 = o.win; pre = w0.t_end <= t_on - 0.06; post = w0.t_start >= t_off + 0.5;     % sano: hasta 3 ciclos antes y desde 0.5 s despues
    R.sc{end+1} = subscan(o.scan, pre); R.meta(end+1, :) = {sprintf('laboratorio %s prefalta', fn(1:end-4)), 'laboratorio', NaN, NaN, NaN, 'laboratorio'};
    rr = w0.roundness;
    rowsA(end+1, :) = {fn(1:end-4), truth, t_on, t_off, nnz(pre), nnz(pre)*STL/fsL, median(rr(pre)), pct_(rr(pre), 1), pct_(rr(pre), 5), pct_(rr(pre), 95), min(rr(pre)), std(rr(pre)), ...
        nnz(post), nnz(post)*STL/fsL, median(rr(post)), pct_(rr(post), 1), min(rr(post)), 100*mean(~strcmp(w0.label_raw(pre), 'A-B-C')), 100*mean(ismember(w0.label_raw(pre), C.GR)), median(w0.pinmax(pre))}; %#ok<AGROW>
    for D = C.DELTAS
        for K = C.KS
            a = r2_decision_chain(o.scan, [], C.opts, 'Delta', D, 'K', K, ev{:}); e = a.event; w = a.win;
            dpost = a.decl(a.decl.t_decl >= t_off + 0.5, :); dmid = a.decl(a.decl.t_decl > t_off & a.decl.t_decl < t_off + 0.5, :);
            dpre = a.decl(a.decl.t_decl <= t_on - 0.06, :);
            rows(end+1, :) = {fn(1:end-4), truth, D, K, e.outcome, e.decl_class, e.decl_ms, e.first_label_ms, e.max_run_correct, e.n_decl_pre, e.pre_classes, height(dpre), ...
                height(dpost), cnt_str(dpost.class), height(dmid), cnt_str(dmid.class), median(w.r_amb(pre), 'omitnan'), median(w.thr(pre), 'omitnan'), e.full_roundness_med, 100*mean(w.fed(pre & w.armed)), ...
                100*mean(w.fed(post)), median(w.r_amb(post), 'omitnan'), nnz(pre & w.armed)*STL/fsL, nnz(post)*STL/fsL}; %#ok<AGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'registro','clase','delta','K','resultado','clase_declarada','decl_ms','primera_etq_ms','racha_max','n_decl_prefalta','clases_prefalta','n_decl_prefalta_hasta_3_ciclos_antes', ...
    'n_decl_tras_despeje','clases_tras_despeje','n_decl_recuperacion_0p5s','clases_recuperacion','r_amb_prefalta_med','umbral_prefalta_med','redondez_falta_med','prefalta_alimenta_pct', ...
    'posfalta_alimenta_pct','r_amb_posfalta_med','exposicion_prefalta_s','exposicion_posfalta_s'});
writetable(T, fullfile(C.OUT, 'lab_detalle.csv'), 'QuoteStrings', true);
TA = cell2table(rowsA, 'VariableNames', {'registro','clase','t_on_s','t_off_s','n_ventanas_prefalta','prefalta_s','redondez_med','redondez_p1','redondez_p5','redondez_p95','redondez_min','redondez_std', ...
    'n_ventanas_posfalta','posfalta_s','redondez_posfalta_med','redondez_posfalta_p1','redondez_posfalta_min','prefalta_cruda_no_sana_pct','prefalta_cruda_tierra_pct','prefalta_pinmax_med'});
writetable(TA, fullfile(C.OUT, 'lab_ambiente.csv'), 'QuoteStrings', true);
R.T = T; R.TA = TA;
isLL = ~strcmp(T.clase, 'A-B-C-G'); ok = strcmp(T.resultado, 'correcta');
fprintf('Ambiente real: redondez sana mediana de %.4f a %.4f; p1 de %.4f a %.4f; minimo %.4f; exposicion sana total %.1f s antes y %.1f s despues de la falta\n', ...
    min(TA.redondez_med), max(TA.redondez_med), min(TA.redondez_p1), max(TA.redondez_p1), min(TA.redondez_min), sum(TA.prefalta_s), sum(TA.posfalta_s));
for D = C.DELTAS
    for K = C.KS
        s = T.delta == D & T.K == K;
        fprintf('  Delta = %.3f K = %3d: L-L %2d/%2d correctas (lat. mediana %5.2f, max %5.2f ms) | trifasicas %d/%d | declaraciones falsas: prefalta %d, recuperacion %d, tras el despeje %d\n', D, K, ...
            nnz(ok & s & isLL), nnz(s & isLL), median(T.decl_ms(ok & s & isLL)), max(T.decl_ms(ok & s & isLL)), nnz(ok & s & ~isLL), nnz(s & ~isLL), sum(T.n_decl_prefalta(s)), sum(T.n_decl_recuperacion_0p5s(s)), sum(T.n_decl_tras_despeje(s)));
    end
end
end

%% ===================================================================== G: dispersion y Delta minimo seguro
function block_G(C, G)
fprintf('\n================ BLOQUE G: dispersion de la redondez sana y Delta minimo seguro ================\n');
SC = [G.S.sc, G.L.sc]; META = [G.S.meta; G.L.meta]; rows = {};
for i = 1:numel(SC)
    S = SC{i}; r = S.roundness; r = r(isfinite(r)); n = numel(r);
    if n < 400, continue, end
    med = median(r); mm = movmax(r, [C.M - 1, 0]); q10 = min(mm(C.M:end));          % minimo del maximo movil de M ventanas
    dmin = nan(1, numel(C.KS));
    for ik = 1:numel(C.KS)
        for D = C.DFINE
            a = r2_decision_chain(S, [], C.opts, 'Delta', D, 'K', C.KS(ik));
            if height(a.decl) == 0, dmin(ik) = D; break, end
        end
    end
    % declaraciones con el Delta de referencia y K = 80, y cuantas son de tierra (no dependen de Delta)
    a = r2_decision_chain(S, [], C.opts, 'Delta', 0.01, 'K', 80); ng = nnz(ismember(a.decl.class, C.GR));
    rows(end+1, :) = [{META{i, 1}, META{i, 6}, META{i, 3}, META{i, 4}, META{i, 5}, n, med, pct_(r, 1), pct_(r, 5), min(r), std(r), 1.4826*median(abs(r - med)), ...
        med - pct_(r, 1), med - pct_(r, 5), med - min(r), med - q10, pct_(r, 95), pct_(r, 99), max(r), q10, pct_(r, 95) - q10, pct_(r, 99) - q10, max(r) - q10}, num2cell(dmin), {max(dmin, [], 'includenan'), ng}]; %#ok<AGROW>
    fprintf('%-52s mediana %.4f | med-p1 %.4f | med-min %.4f | med-q(racha de %d) %.4f | Delta minimo seguro por K: %s | decl. de tierra %d\n', META{i, 1}, med, med - pct_(r, 1), med - min(r), C.M, med - q10, mat2str(dmin), ng);
end
kn = arrayfun(@(K) sprintf('delta_min_seguro_K%d', K), C.KS, 'UniformOutput', false);
T = cell2table(rows, 'VariableNames', [{'traza','origen','v2_pct','fase_v2_grados','snr_db','n_ventanas','redondez_med','redondez_p1','redondez_p5','redondez_min','redondez_std','redondez_sigma_mad', ...
    'med_menos_p1','med_menos_p5','med_menos_min','med_menos_min_del_max_movil_M','redondez_p95','redondez_p99','redondez_max','min_del_max_movil_M', ...
    'p95_menos_min_del_max_movil_M','p99_menos_min_del_max_movil_M','max_menos_min_del_max_movil_M'}, kn, {'delta_min_seguro_todas_K','n_decl_tierra_delta001_K80'}]);
writetable(T, fullfile(C.OUT, 'dispersion_delta_min.csv'), 'QuoteStrings', true);
ok = isfinite(T.delta_min_seguro_todas_K);
fprintf('Trazas: %d. Delta minimo seguro: maximo %.4f. Trazas en que med - p1 >= Delta minimo: %d de %d; med - min >= Delta minimo: %d de %d\n', height(T), max(T.delta_min_seguro_todas_K), ...
    nnz(T.med_menos_p1(ok) >= T.delta_min_seguro_todas_K(ok) - 1e-12), nnz(ok), nnz(T.med_menos_min(ok) >= T.delta_min_seguro_todas_K(ok) - 1e-12), nnz(ok));
cands = {'med_menos_p1','med_menos_min','med_menos_min_del_max_movil_M','p95_menos_min_del_max_movil_M','p99_menos_min_del_max_movil_M','max_menos_min_del_max_movil_M'};
for q = 1:numel(cands)
    x = T.(cands{q}); d = T.delta_min_seguro_todas_K;
    fprintf('  candidato %-34s: cubre el Delta minimo (con la resolucion de 0.001 de la rejilla) en %2d de %2d trazas; exceso mediano %.4f, maximo %.4f, minimo %.4f\n', cands{q}, nnz(x(ok) >= d(ok) - 0.001 - 1e-12), nnz(ok), median(x(ok) - d(ok)), max(x(ok) - d(ok)), min(x(ok) - d(ok)));
end
sim = strcmp(T.origen, 'simulacion');
fprintf('  Delta minimo seguro, simulacion: K=20 max %.4f; K=40 max %.4f; K=80 max %.4f; K=160 max %.4f. Laboratorio: max %.4f\n', max(T.delta_min_seguro_K20(sim)), max(T.delta_min_seguro_K40(sim)), max(T.delta_min_seguro_K80(sim)), max(T.delta_min_seguro_K160(sim)), max(T.delta_min_seguro_todas_K(~sim)));
end

%% ===================================================================== I: inicializacion, congelacion y rearme
function block_I(C)
fprintf('\n================ BLOQUE I: inicializacion, congelacion y rearme ================\n');
info = C.info; Tc = C.Tc; fs = C.fs; f = C.f;
%% I1. La falta llega antes de completar K ventanas
fprintf('\n-- I1. Falta que llega antes de completar las K ventanas de inicializacion\n');
LEADS = [0 25 50 75 100 125 150 250 500]; rows = {};
for i = 1:height(Tc)
    d = load(fullfile(info.dir, Tc.file{i}));
    for c = 1:numel(C.CONV)
        [V, t] = extract_converter_data(d, C.CONV{c});
        ev = {'truth', Tc.truth{i}, 't_on', info.t_on, 't_off', info.t_off};
        o = r2_decision_chain(V, t, r2_opts_nominal(C.opts, max(abs(V(t >= 2.4 & t < 2.8, :)), [], 'all')), ev{:}, 't_range', [2.0 5.0]);
        ra0 = median(o.scan.roundness(o.scan.t_end <= info.t_on & o.scan.t_start >= 2.4));
        for lead = LEADS
            S = subscan(o.scan, o.scan.t_start >= info.t_on - lead/1e3 - 1e-9);
            for K = C.KS
                for mode = 1:2
                    if mode == 1, extra = {}; mname = 'arranque en frio (se supone sano)';
                    else, extra = {'amb_init', ra0}; mname = 'arranque con nivel guardado'; end
                    a = r2_decision_chain(S, [], C.opts, 'K', K, ev{:}, extra{:}); e = a.event; w = a.win;
                    flt = w.t_start >= info.t_on & w.t_end <= info.t_off; ja = find(w.armed, 1);
                    aft = a.decl(a.decl.t_decl > info.t_off + 0.2, :);
                    rows(end+1, :) = {Tc.file{i}, Tc.truth{i}, Tc.family{i}, C.CONV{c}, lead, K, mname, e.outcome, e.decl_class, e.decl_ms, 1e3*(w.t_end(ja) - info.t_on), ...
                        100*mean(w.fed(flt)), w.r_amb(find(flt, 1, 'last')), ra0, e.n_decl_after, e.after_classes, height(aft), cnt_str(aft.class), w.r_amb(end)}; %#ok<AGROW>
                end
            end
        end
    end
end
T1 = cell2table(rows, 'VariableNames', {'archivo','clase','familia','convertidor','sano_previo_ms','K','modo','resultado','clase_declarada','decl_ms','armado_ms_respecto_falta', ...
    'ventanas_falta_que_alimentan_pct','r_amb_fin_falta','r_amb_sano','n_decl_posfalta','clases_posfalta','n_decl_desde_200ms_tras_despeje','clases_desde_200ms','r_amb_final'});
writetable(T1, fullfile(C.OUT, 'init_falta_temprana.csv'), 'QuoteStrings', true);
okk = strcmp(T1.resultado, 'correcta'); frio = startsWith(T1.modo, 'arranque en frio');
for K = C.KS
    fprintf('  K = %3d, arranque en frio: correctas de 24 por sano previo (ms) ', K);
    for lead = LEADS, s = T1.K == K & T1.sano_previo_ms == lead & frio; fprintf('%d:%d ', lead, nnz(okk & s)); end
    s = T1.K == K & frio; fprintf('| latencia max %.1f ms | decl. >200 ms tras el despeje: %d\n', max(T1.decl_ms(okk & s)), sum(T1.n_decl_desde_200ms_tras_despeje(s)));
    fprintf('  K = %3d, nivel guardado:   correctas de 24 por sano previo (ms) ', K);
    for lead = LEADS, s = T1.K == K & T1.sano_previo_ms == lead & ~frio; fprintf('%d:%d ', lead, nnz(okk & s)); end
    s = T1.K == K & ~frio; fprintf('| latencia max %.1f ms | decl. >200 ms tras el despeje: %d\n', max(T1.decl_ms(okk & s)), sum(T1.n_decl_desde_200ms_tras_despeje(s)));
end

%% I2. Recuperacion real tras el despeje (registros completos hasta 5 s)
fprintf('\n-- I2. Rearme tras el despeje en los 24 registros (prefalta desde 2.0 s, registro hasta 5.0 s)\n');
rows = {};
for i = 1:height(Tc)
    d = load(fullfile(info.dir, Tc.file{i}));
    for c = 1:numel(C.CONV)
        [V, t] = extract_converter_data(d, C.CONV{c});
        ev = {'truth', Tc.truth{i}, 't_on', info.t_on, 't_off', info.t_off};
        o = r2_decision_chain(V, t, r2_opts_nominal(C.opts, max(abs(V(t >= 2.4 & t < 2.8, :)), [], 'all')), ev{:}, 't_range', [2.0 5.0]);
        for D = C.DELTAS
            for K = C.KS
                a = r2_decision_chain(o.scan, [], C.opts, 'Delta', D, 'K', K, ev{:}); e = a.event; w = a.win;
                flt = w.t_start >= info.t_on & w.t_end <= info.t_off; post = w.t_start > info.t_off;
                pre = w.armed & w.t_end <= info.t_on; ra_pre = median(w.r_amb(pre), 'omitnan');
                jf = find(post & w.fed, 1); if isempty(jf), t_fed = NaN; else, t_fed = 1e3*(w.t_end(jf) - info.t_off); end
                jr = find(post & abs(w.r_amb - ra_pre) < 0.002 & w.fed, 1); if isempty(jr), t_rec = NaN; else, t_rec = 1e3*(w.t_end(jr) - info.t_off); end
                jh = find(post & ~strcmp(w.label, 'A-B-C'), 1, 'last'); if isempty(jh), t_lastnh = 0; else, t_lastnh = 1e3*(w.t_end(jh) - info.t_off); end
                late = a.decl(a.decl.t_decl > info.t_off + 0.2, :);
                rows(end+1, :) = {Tc.file{i}, Tc.truth{i}, Tc.family{i}, C.CONV{c}, D, K, e.outcome, e.decl_ms, 100*mean(w.fed(flt)), ra_pre, w.r_amb(find(flt, 1, 'last')), ...
                    t_fed, t_rec, t_lastnh, e.n_decl_after, e.after_classes, height(late), cnt_str(late.class), w.r_amb(end), 100*mean(w.fed(post)), 100*mean(~strcmp(w.label(w.t_start > info.t_off + 0.2), 'A-B-C'))}; %#ok<AGROW>
            end
        end
    end
end
T2 = cell2table(rows, 'VariableNames', {'archivo','clase','familia','convertidor','delta','K','resultado','decl_ms','ventanas_falta_que_alimentan_pct','r_amb_prefalta','r_amb_fin_falta', ...
    'primera_actualizacion_tras_despeje_ms','r_amb_recuperado_ms','ultima_etiqueta_no_sana_ms','n_decl_posfalta','clases_posfalta','n_decl_desde_200ms','clases_desde_200ms','r_amb_final','posfalta_alimenta_pct','no_sana_desde_200ms_pct'});
writetable(T2, fullfile(C.OUT, 'rearme_tras_despeje.csv'), 'QuoteStrings', true);
for D = C.DELTAS
    for K = C.KS
        s = T2.delta == D & T2.K == K;
        fprintf('  Delta = %.3f K = %3d: alimenta durante la falta (max) %.1f %% | r_amb recuperado: mediana %.0f ms, max %.0f ms, sin recuperar %d | registros con declaracion >200 ms tras el despeje: %d (%s)\n', D, K, ...
            max(T2.ventanas_falta_que_alimentan_pct(s)), median(T2.r_amb_recuperado_ms(s), 'omitnan'), max(T2.r_amb_recuperado_ms(s)), nnz(isnan(T2.r_amb_recuperado_ms(s))), nnz(T2.n_decl_desde_200ms(s) > 0), cnt_str(T2.clases_desde_200ms(s & T2.n_decl_desde_200ms > 0)));
    end
end

%% I3. Escalon legitimo del desequilibrio ambiente
fprintf('\n-- I3. Escalon del desequilibrio ambiente V2 (no es una falta) en t = 1.2 s, 1.75 s de registro posterior\n');
base = load(fullfile(info.dir, 'DataGridRL_AB.mat')); rows = {};
STEPS = [0 0.0025; 0 0.005; 0 0.01; 0 0.02; 0.01 0.02; 0.005 0.01; 0.02 0.01; 0.01 0];
for c = 1:numel(C.CONV)
    [V, t] = extract_converter_data(base, C.CONV{c}); npk = max(abs(V(t >= 2.4 & t < 2.8, :)), [], 'all');
    seg = t >= 0.4 & t < 2.95; Vh = V(seg, :); th = t(seg); tS = 1.2; af = th >= tS;
    for q = 1:size(STEPS, 1)
        Vs = Vh + v2wave(STEPS(q, 1), npk, f, th, 0).*(~af) + v2wave(STEPS(q, 2), npk, f, th, 0).*af;
        o = r2_decision_chain(Vs, th, r2_opts_nominal(C.opts, npk));
        for D = C.DELTAS
            for K = C.KS
                a = r2_decision_chain(o.scan, [], C.opts, 'Delta', D, 'K', K); w = a.win; post = w.t_start >= tS;
                dp = a.decl(a.decl.t_decl > tS, :); dpre = a.decl(a.decl.t_decl <= tS, :);
                jl = find(post & w.fire, 1, 'last'); if isempty(jl), t_last = NaN; else, t_last = w.t_end(jl) - tS; end
                rows(end+1, :) = {C.CONV{c}, 100*STEPS(q, 1), 100*STEPS(q, 2), D, K, height(dpre), height(dp), cnt_str(dp.class), 100*mean(w.fire(post)), t_last, 100*mean(w.fed(post)), ...
                    median(w.r_amb(w.armed & ~post), 'omitnan'), w.r_amb(end), median(w.roundness(post)), double(~isnan(t_last) && t_last >= (w.t_end(end) - tS) - 0.05)}; %#ok<AGROW>
            end
        end
    end
end
T3 = cell2table(rows, 'VariableNames', {'convertidor','v2_antes_pct','v2_despues_pct','delta','K','n_decl_antes','n_decl_despues','clases','ventanas_aseguradas_despues_pct','ultima_asegurada_s', ...
    'despues_alimenta_pct','r_amb_antes','r_amb_final','redondez_despues_med','declaracion_permanente'});
writetable(T3, fullfile(C.OUT, 'escalon_ambiente.csv'), 'QuoteStrings', true);
for q = 1:size(STEPS, 1)
    fprintf('  V2 %.2f -> %.2f %%:', 100*STEPS(q, 1), 100*STEPS(q, 2));
    for D = C.DELTAS
        s = T3.v2_antes_pct == 100*STEPS(q, 1) & T3.v2_despues_pct == 100*STEPS(q, 2) & T3.delta == D;
        fprintf('  D=%.3f: decl. en %d/%d casos, hasta el final en %d', D, nnz(T3.n_decl_despues(s) > 0), nnz(s), nnz(T3.declaracion_permanente(s) > 0));
    end
    fprintf('\n');
end

%% I4. Falta de evolucion lenta (rampa de severidad) y segunda falta tras el despeje
fprintf('\n-- I4. Falta A-B sintetica con severidad en rampa de 0 a 0.20\n');
rows = {};
for Tr = [0.05 0.2 1 2 5 10]
    n = round((0.5 + Tr + 0.5)*fs); ts = (0:n-1)'/fs;
    vb = [cos(2*pi*f*ts), cos(2*pi*f*ts - 2*pi/3), cos(2*pi*f*ts + 2*pi/3)];
    k = 0.20*min(max((ts - 0.5)/Tr, 0), 1); dv = (k/2).*(vb(:, 1) - vb(:, 2));
    V = vb; V(:, 1) = V(:, 1) - dv; V(:, 2) = V(:, 2) + dv;
    ev = {'truth', 'A-B', 't_on', 0.5};
    o = r2_decision_chain(V, ts, r2_opts_nominal(C.opts, 1), ev{:});
    for D = C.DELTAS
        for K = C.KS
            a = r2_decision_chain(o.scan, [], C.opts, 'Delta', D, 'K', K, ev{:}); e = a.event; w = a.win;
            if e.declared, kd = 0.20*min(e.decl_ms/1e3/Tr, 1); else, kd = NaN; end
            rows(end+1, :) = {Tr, 0.20/Tr, D, K, e.outcome, e.decl_class, e.decl_ms, kd, 100*kd/2, w.r_amb(end), 100*mean(w.fed(w.t_start >= 0.5)), 100*mean(w.fire(w.t_start >= 0.5 + Tr))}; %#ok<AGROW>
        end
    end
end
T4 = cell2table(rows, 'VariableNames', {'rampa_s','pendiente_k_por_s','delta','K','resultado','clase_declarada','decl_ms','k_en_la_declaracion','v2v1_en_la_declaracion_pct','r_amb_final','ventanas_falta_que_alimentan_pct','aseguradas_tras_rampa_pct'});
writetable(T4, fullfile(C.OUT, 'falta_lenta.csv'), 'QuoteStrings', true);
for Tr = unique(T4.rampa_s)'
    fprintf('  rampa de %5.2f s:', Tr);
    for K = C.KS
        s = T4.rampa_s == Tr & T4.K == K & T4.delta == 0.01;
        fprintf('  K=%3d: %s k=%.3f', K, T4.resultado{s}, T4.k_en_la_declaracion(s));
    end
    fprintf('   (Delta = 0.01)\n');
end

fprintf('\n-- I4b. Dos faltas sinteticas sucesivas (A-B de 0.5 a 0.66 s y B-C desde 1.2 s, k = 0.3)\n');
n = round(1.6*fs); ts = (0:n-1)'/fs; vb = [cos(2*pi*f*ts), cos(2*pi*f*ts - 2*pi/3), cos(2*pi*f*ts + 2*pi/3)];
V = vb; m1 = ts >= 0.5 & ts < 0.66; m2 = ts >= 1.2; kk = 0.3;
dv = (kk/2)*(vb(:, 1) - vb(:, 2)); V(m1, 1) = V(m1, 1) - dv(m1); V(m1, 2) = V(m1, 2) + dv(m1);
dv = (kk/2)*(vb(:, 2) - vb(:, 3)); V(m2, 2) = V(m2, 2) - dv(m2); V(m2, 3) = V(m2, 3) + dv(m2);
o = r2_decision_chain(V, ts, r2_opts_nominal(C.opts, 1)); rows = {};
for D = C.DELTAS
    for K = C.KS
        a = r2_decision_chain(o.scan, [], C.opts, 'Delta', D, 'K', K); w = a.win; dc = a.decl;
        d1 = dc(dc.t_decl > 0.5 & dc.t_decl <= 0.66, :); d2 = dc(dc.t_decl > 1.2, :); dmid = dc(dc.t_decl > 0.66 + 0.01 & dc.t_decl <= 1.2, :);
        l1 = NaN; c1 = ''; l2 = NaN; c2 = '';
        if ~isempty(d1), l1 = 1e3*(d1.t_decl(1) - 0.5); c1 = d1.class{1}; end
        if ~isempty(d2), l2 = 1e3*(d2.t_decl(1) - 1.2); c2 = d2.class{1}; end
        jf = find(w.t_start >= 0.66 & w.fed, 1); if isempty(jf), tf_ = NaN; else, tf_ = 1e3*(w.t_end(jf) - 0.66); end
        rows(end+1, :) = {D, K, c1, l1, c2, l2, height(dmid), 100*mean(w.fed(w.t_start >= 0.5 & w.t_end <= 0.66)), tf_, w.r_amb(find(w.t_end <= 1.2, 1, 'last'))}; %#ok<AGROW>
    end
end
T5 = cell2table(rows, 'VariableNames', {'delta','K','clase_falta_1','decl_1_ms','clase_falta_2','decl_2_ms','n_decl_entre_faltas','ventanas_falta_1_que_alimentan_pct','primera_actualizacion_tras_despeje_ms','r_amb_antes_de_la_segunda'});
writetable(T5, fullfile(C.OUT, 'dos_faltas.csv'), 'QuoteStrings', true); disp(T5);
end

%% ===================================================================== resumen de la rejilla
function grid_summary(C, G)
fprintf('\n================ RESUMEN DE LA REJILLA (sens_delta_K.csv) ================\n');
TS = G.S.T; TF = G.F.T; TFl = G.D.TF; TL = G.L.T; rows = {};
okF = strcmp(TF.resultado, 'correcta'); okL = strcmp(TL.resultado, 'correcta'); isLLlab = ~strcmp(TL.clase, 'A-B-C-G');
pub = ismember(TF.familia, {'LL_RL'}) & ismember(TF.v2_pct, [0 1]) & TF.fase_v2_grados == 0; nom = TF.v2_pct == 0;
LLf = ismember(TF.familia, {'LL', 'LL_RL'});
for D = C.DELTAS
    for K = C.KS
        s = TS.delta == D & TS.K == K; cl = isinf(TS.snr_db); pb = cl & TS.fase_v2_grados == 0 & ismember(TS.v2_pct, [0 0.5 1 2]);
        f_ = TF.delta == D & TF.K == K; l = TL.delta == D & TL.K == K;
        fl = @(k2) TFl.k_suelo_peor_fase(TFl.delta == D & TFl.K == K & TFl.v2_ambiente_pct == k2);
        fl0 = @(k2) TFl.k_suelo_fase0(TFl.delta == D & TFl.K == K & TFl.v2_ambiente_pct == k2);
        gfl = strcmp(TF.convertidor, 'Pex'); v2p = TF.v2_pct > 0;
        rows(end+1, :) = {D, K, 1e3*((K-1)*C.STRIDE + C.W)/C.fs, ((K-1)*C.STRIDE + C.W)/C.fs*C.f, ...
            nnz(s & pb), nnz(TS.n_declaraciones(s & pb) > 0), max(TS.ventanas_aseguradas_pct(s & pb)), ...
            nnz(s & cl), nnz(TS.n_declaraciones(s & cl) > 0), sum(TS.n_declaraciones(s & cl)), ...
            nnz(s & ~cl), nnz(TS.n_declaraciones(s & ~cl) > 0), sum(TS.n_declaraciones(s & ~cl)), ...
            nnz(f_ & pub), nnz(okF & f_ & pub), min(TF.decl_ms(okF & f_ & pub)), max(TF.decl_ms(okF & f_ & pub)), ...
            nnz(f_ & nom), nnz(okF & f_ & nom), median(TF.decl_ms(okF & f_ & nom)), max(TF.decl_ms(okF & f_ & nom)), ...
            nnz(f_ & LLf & v2p & gfl), nnz(okF & f_ & LLf & v2p & gfl), nnz(f_ & LLf & v2p & ~gfl), nnz(okF & f_ & LLf & v2p & ~gfl), ...
            nnz(f_ & LLf & TF.v2_pct > 0 & TF.v2_pct <= 1 & ~gfl), nnz(okF & f_ & LLf & TF.v2_pct > 0 & TF.v2_pct <= 1 & ~gfl), ...
            nnz(f_ & ~LLf & v2p), nnz(okF & f_ & ~LLf & v2p), sum(TF.n_decl_prefalta(f_)), sum(TF.n_decl_erroneas(f_)), ...
            fl0(0), fl0(1), fl(1), fl(2), ...
            nnz(l & isLLlab), nnz(okL & l & isLLlab), median(TL.decl_ms(okL & l & isLLlab)), max(TL.decl_ms(okL & l & isLLlab)), nnz(okL & l & ~isLLlab), ...
            sum(TL.n_decl_prefalta(l)), sum(TL.n_decl_tras_despeje(l)), sum(TL.n_decl_recuperacion_0p5s(l))}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'delta','K','K_ms','K_ciclos', ...
    'sano_pub_trazas','sano_pub_trazas_con_decl','sano_pub_aseguradas_pct_max', ...
    'sano_fases_trazas','sano_fases_trazas_con_decl','sano_fases_n_decl', ...
    'sano_ruido_trazas','sano_ruido_trazas_con_decl','sano_ruido_n_decl', ...
    'faltas_pub_n','faltas_pub_correctas','faltas_pub_lat_min_ms','faltas_pub_lat_max_ms', ...
    'faltas_nominales_n','faltas_nominales_correctas','faltas_nominales_lat_med_ms','faltas_nominales_lat_max_ms', ...
    'LL_v2_GFL_n','LL_v2_GFL_correctas','LL_v2_GFM_n','LL_v2_GFM_correctas','LL_v2_hasta1_GFM_n','LL_v2_hasta1_GFM_correctas', ...
    'tierra_v2_n','tierra_v2_correctas','faltas_n_decl_prefalta','faltas_n_decl_erroneas', ...
    'suelo_k_v2_0','suelo_k_v2_1_fase0','suelo_k_v2_1_peor_fase','suelo_k_v2_2_peor_fase', ...
    'lab_LL_n','lab_LL_correctas','lab_LL_lat_med_ms','lab_LL_lat_max_ms','lab_trifasicas_correctas','lab_n_decl_prefalta','lab_n_decl_tras_despeje','lab_n_decl_recuperacion'});
writetable(T, fullfile(C.OUT, 'sens_delta_K.csv'), 'QuoteStrings', true);
disp(T(:, {'delta','K','sano_pub_trazas_con_decl','sano_fases_trazas_con_decl','sano_ruido_trazas_con_decl','faltas_pub_correctas','faltas_nominales_correctas','LL_v2_GFL_correctas','LL_v2_GFM_correctas','LL_v2_hasta1_GFM_correctas','suelo_k_v2_0','suelo_k_v2_1_peor_fase','lab_LL_correctas','lab_n_decl_prefalta','lab_n_decl_tras_despeje'}));
end
