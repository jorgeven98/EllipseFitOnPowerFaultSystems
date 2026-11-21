function r2_final_emt_resumen(OUT)
% R2_FINAL_EMT_RESUMEN  Resumenes, controles, cifras e informe del postproceso
% definitivo de la campana EMT (T90, bloque EMT). No reclasifica (salvo el
% diagnostico ventana a ventana de las diferencias del control a1): lee los
% CSV por evento que escribe r2_final_emt (eventos_todos.csv, S2_ruido_eventos.csv)
% y, como referencia de la puerta antigua, las columnas T03_* de
% ronda2\resultados\T23\S1_preliminar.csv y S2_preliminar.csv y las filas G5_005
% de ronda2\resultados\T52\emt_eventos.csv.
%
% Salidas en OUT (por defecto ronda2\resultados\T90\emt):
%   controles.csv, controles_diferencias.csv, control_a1_diagnostico.csv, cambios_vs_T23.csv
%   S2_umbrales.csv, S2_umbrales_T23.csv, S2_barrido_tabla.csv, S2_recuentos.csv, S2_incidencia.csv,
%   S2_tierra_resistiva.csv, S2_latencias.csv
%   S2_ruido_resumen.csv, S2_ruido_umbrales.csv, S2_ruido_cambios.csv
%   S1_por_factor.csv, S1_totales.csv, S1_angulo.csv, S1_no_correctas.csv
%   y, a traves de r2_final_emt_cifras, cifras_emt.csv y seccion_informe_emt.md
% Los valores "antiguos" se calculan aqui con las mismas funciones de resumen
% sobre las columnas T03_* de T23 (puerta publicada) y la misma poblacion (Pex y
% DG1, medidas validas); no se transcriben de los informes.
here = fileparts(mfilename('fullpath')); RAIZ = fileparts(here);
R2 = fullfile(RAIZ, '_revision_TPWRD', 'ronda2');
if nargin < 1 || isempty(OUT), OUT = fullfile(R2, 'resultados', 'T90', 'emt'); end
addpath(here, fullfile(RAIZ, 'src', 'fault_classifiers'), fullfile(RAIZ, 'src', 'utils'));
if strcmp(get(0, 'Diary'), 'off'), diary(fullfile(OUT, 'r2_final_emt_log.txt')); cu = onCleanup(@() diary('off')); end
fprintf('\n==== R2_FINAL_EMT_RESUMEN  %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
ZBASE = 5.775;
E = r2_final_emt_leer(fullfile(OUT, 'eventos_todos.csv'));
A = antiguo_(E, R2);                      % mismas filas con el resultado de la puerta publicada (T23)
Ev = E(E.valida == 1, :); Av = A(A.valida == 1, :);
fprintf('Eventos (caso y terminal): %d; validos %d; excluidos %d\n', height(E), height(Ev), height(E) - height(Ev));

%% ---- controles
[C, D] = controles_(E, A, R2);
writetable(C, fullfile(OUT, 'controles.csv')); writetable(D, fullfile(OUT, 'controles_diferencias.csv'));
diag_a1_(D(strcmp(D.control, 'a1'), :), R2, OUT);
disp(C);
% cambios respecto a T23 (puerta publicada), solo medidas validas
dif = ~strcmp(Ev.outcome, Av.outcome) | ~strcmp(Ev.decl_class, Av.decl_class) | ~igual_(Ev.decl_ms, Av.decl_ms);
X = Ev(dif, {'id', 'conjunto', 'terminal', 'tipo', 'Rf_ohm', 'Rtierra_ohm', 'ang_inc_pedido_deg', 'V0V1_falta_pct', 'V2V1_falta_pct', 'full_v0stat_med'});
X.T23_outcome = Av.outcome(dif); X.T23_decl_class = Av.decl_class(dif); X.T23_decl_ms = Av.decl_ms(dif);
X.final_outcome = Ev.outcome(dif); X.final_decl_class = Ev.decl_class(dif); X.final_decl_ms = Ev.decl_ms(dif);
X.cambia_resultado = double(~strcmp(X.T23_outcome, X.final_outcome) | ~strcmp(X.T23_decl_class, X.final_decl_class));
writetable(X, fullfile(OUT, 'cambios_vs_T23.csv'));
fprintf('Eventos validos que cambian respecto a T23 (puerta publicada): %d (resultado o clase: %d; solo latencia: %d)\n', height(X), nnz(X.cambia_resultado), nnz(~X.cambia_resultado));

%% ---- S2: umbrales, tabla del barrido, recuentos, incidencia, tierra resistiva, latencias
sw = ismember(Ev.conjunto, {'S2a', 'S2d'});
U = umbrales_(Ev(sw, :), ZBASE); U0 = umbrales_(Av(sw, :), ZBASE);
U = [U, renombra_(U0(:, {'secuencia', 'Rf_ok_contigua_ohm', 'Rf_ok_max_ohm', 'Rf_decl_contigua_ohm', 'Rf_decl_max_ohm', 'Rf_nodecl_desde_ohm', 'clases_erroneas'}), 'T23_')];
writetable(U, fullfile(OUT, 'S2_umbrales.csv')); writetable(U0, fullfile(OUT, 'S2_umbrales_T23.csv'));
writetable(barrido_(Ev(sw, :), Av(sw, :)), fullfile(OUT, 'S2_barrido_tabla.csv'));
RC = recuentos_s2_(Ev, Av); writetable(RC, fullfile(OUT, 'S2_recuentos.csv'));
INC = incidencia_(Ev, Av); writetable(INC, fullfile(OUT, 'S2_incidencia.csv'));
s2c = strcmp(Ev.conjunto, 'S2c');
TR = Ev(s2c, {'id', 'terminal', 'tipo', 'Rf_ohm', 'Rtierra_ohm', 'V0V1_falta_pct', 'V2V1_falta_pct', 'full_v0stat_med', 'outcome', 'decl_class', 'decl_ms', 'first_label_ms', 'correct_decl_ms', 'full_ground_pct', 'full_modal_label'});
TR.T23_outcome = Av.outcome(s2c); TR.T23_decl_class = Av.decl_class(s2c); TR.T23_decl_ms = Av.decl_ms(s2c);
TR = sortrows(TR, {'terminal', 'Rtierra_ohm'}, {'descend', 'ascend'});
writetable(TR, fullfile(OUT, 'S2_tierra_resistiva.csv'));
LAT = latencias_(Ev, Av); writetable(LAT, fullfile(OUT, 'S2_latencias.csv'));

%% ---- S2 con ruido de 40 dB y clase 0.5
fr = fullfile(OUT, 'S2_ruido_eventos.csv'); NR = table(); NU = table();
if exist(fr, 'file')
    N = r2_final_emt_leer(fr); N = N(N.valida == 1, :);
    [NR, NU, NC] = ruido_(N, ZBASE);
    writetable(NR, fullfile(OUT, 'S2_ruido_resumen.csv')); writetable(NU, fullfile(OUT, 'S2_ruido_umbrales.csv')); writetable(NC, fullfile(OUT, 'S2_ruido_cambios.csv'));
    % coherencia: las filas 'base' de la parte de ruido deben coincidir con S2_final
    Nb = N(strcmp(N.condicion, 'base'), :);
    [tf, loc] = ismember(strcat(Nb.id, '|', Nb.terminal), strcat(Ev.id, '|', Ev.terminal));
    nd = nnz(~tf) + nnz(~strcmp(Nb.outcome(tf), Ev.outcome(loc(tf))) | ~strcmp(Nb.decl_class(tf), Ev.decl_class(loc(tf))) | ~igual_(Nb.decl_ms(tf), Ev.decl_ms(loc(tf))));
    fprintf('Coherencia interna: filas base de S2_ruido_eventos frente a S2_final: %d diferencias en %d eventos\n', nd, height(Nb));
end

%% ---- S1
s1 = startsWith(E.conjunto, 'S1');
[SF, ST] = s1_factor_(E(s1, :), A(s1, :)); writetable(SF, fullfile(OUT, 'S1_por_factor.csv')); writetable(ST, fullfile(OUT, 'S1_totales.csv'));
[~, ST0] = s1_factor_(A(s1, :), A(s1, :)); writetable(ST0, fullfile(OUT, 'S1_totales_T23.csv'));   % T90: los mismos totales con la puerta publicada (T23), para r2_final_emt_cifras
SA = s1_angulo_(Ev(startsWith(Ev.conjunto, 'S1'), :), Av(startsWith(Av.conjunto, 'S1'), :)); writetable(SA, fullfile(OUT, 'S1_angulo.csv'));
nd = startsWith(Ev.conjunto, 'S1') & ~strcmp(Ev.outcome, 'correcta');
SN = Ev(nd, {'id', 'conjunto', 'factor', 'terminal', 'tipo', 'Rf_ohm', 'R_bucle_ohm', 'Lf_H', 'red_Z_ohm', 'red_XR', 'red_f_Hz', 'dg1_on', 'dg2_on', 'V2V1_falta_pct', 'theta_ley_deg', 'theta_emt_deg', 'outcome', 'resultado4', 'decl_class', 'decl_ms', 'full_modal_label', 'retention_pct'});
writetable(sortrows(SN, {'terminal', 'resultado4', 'id'}, {'descend', 'ascend', 'ascend'}), fullfile(OUT, 'S1_no_correctas.csv'));
disp(U(:, [1:3 5:17])); disp(RC(:, 1:9)); disp(ST); disp(SA);
if ~isempty(NR), disp(NR(:, [1:5 7:14])); disp(NU); end

%% ---- cifras e informe
if exist('r2_final_emt_cifras', 'file') == 2, r2_final_emt_cifras(OUT); end
fprintf('FIN r2_final_emt_resumen\n');
end

% =========================================================================
function A = antiguo_(E, R2)
% Tabla con las mismas filas que E y los resultados de la puerta publicada
% (columnas T03_* de T23, cadena de T03 con sus valores por defecto).
B1 = r2_final_emt_leer(fullfile(R2, 'resultados', 'T23', 'S1_preliminar.csv')); B2 = r2_final_emt_leer(fullfile(R2, 'resultados', 'T23', 'S2_preliminar.csv'));
cols = {'id', 'convertidor', 'medida_valida', 'V2V1_falta_pct', 'V0V1_falta_pct', 'theta_pred_deg', 'theta_mediana_deg', 'T21_theta_pred_Pex_deg', ...
    'T03_outcome', 'T03_decl_class', 'T03_decl_ms', 'T03_first_label_ms', 'T03_correct_decl_ms', 'T03_n_decl_pre', 'T03_n_decl_correct', 'T03_n_decl_wrong', ...
    'T03_wrong_classes', 'T03_n_decl_after', 'T03_full_ground_pct', 'T03_full_modal_label', 'T03_retention_pct'};
B = [B1(:, cols); B2(:, cols)];
[tf, loc] = ismember(strcat(E.id, '|', E.terminal), strcat(B.id, '|', B.convertidor));
if ~all(tf), error('r2_final_emt_resumen:ref', 'Faltan %d eventos en la linea base de T23.', nnz(~tf)); end
B = B(loc, :); A = E;
A.outcome = B.T03_outcome; A.decl_class = B.T03_decl_class; A.decl_ms = B.T03_decl_ms;
A.first_label_ms = B.T03_first_label_ms; A.correct_decl_ms = B.T03_correct_decl_ms;
A.n_decl_pre = B.T03_n_decl_pre; A.n_decl_correct = B.T03_n_decl_correct; A.n_decl_wrong = B.T03_n_decl_wrong;
A.wrong_classes = B.T03_wrong_classes; A.n_decl_after = B.T03_n_decl_after;
A.full_ground_pct = B.T03_full_ground_pct; A.full_modal_label = B.T03_full_modal_label; A.retention_pct = B.T03_retention_pct;
A.theta_emt_deg = B.theta_mediana_deg; A.theta_ley_deg = B.theta_pred_deg;
A.err_ley_deg = abs(mod(B.theta_mediana_deg - B.theta_pred_deg + 90, 180) - 90);
A.err_T21_deg = abs(mod(B.theta_mediana_deg - B.T21_theta_pred_Pex_deg + 90, 180) - 90);
A.V2V1_falta_pct = B.V2V1_falta_pct; A.V0V1_falta_pct = B.V0V1_falta_pct;
A.valida_T23 = B.medida_valida;
LL = {'A-B', 'B-C', 'C-A'}; A.resultado4 = A.outcome;
for k = 1:height(A)
    if strcmp(A.outcome{k}, 'erronea')
        if ismember(A.tipo{k}, LL) && ismember(A.decl_class{k}, LL), A.resultado4{k} = 'par_adyacente'; else, A.resultado4{k} = 'erronea_otra'; end
    end
end
nv = nnz((A.valida_T23 == 1) ~= (E.valida == 1));
fprintf('Validez del indice frente a medida_valida de T23: %d filas distintas\n', nv);
fprintf('Fasores frente a T23: max |dif| de V2/V1 = %.2e %%, de V0/V1 = %.2e %%; max |dif| de theta_ley = %.2e grados\n', ...
    max(abs(E.V2V1_falta_pct - B.V2V1_falta_pct)), max(abs(E.V0V1_falta_pct - B.V0V1_falta_pct)), max(abs(mod(E.theta_ley_deg - B.theta_pred_deg + 90, 180) - 90)));
end

% =========================================================================
function [C, D] = controles_(E, A, R2)
% (a) faltas entre fases: mismo resultado por evento que T23 (T03_*) y que G5_005 de T52
% (b) faltas a tierra: mismo resultado que G5_005 de T52
T52 = r2_final_emt_leer(fullfile(R2, 'resultados', 'T52', 'emt_eventos.csv'));
T52 = T52(strcmp(T52.variante, 'G5_005') & ismember(T52.convertidor, {'Pex', 'DG1'}), :);
LL = ismember(E.tipo, {'A-B', 'B-C', 'C-A'});
[tf, loc] = ismember(strcat(E.id, '|', E.terminal), strcat(T52.id, '|', T52.convertidor));
G = E; G.outcome(:) = {''}; G.decl_class(:) = {''}; G.decl_ms(:) = NaN;
G.outcome(tf) = T52.resultado(loc(tf)); G.decl_class(tf) = T52.clase(loc(tf)); G.decl_ms(tf) = T52.decl_ms(loc(tf));
todos = true(height(E), 1);
defs = {'a1', 'entre fases frente a T23 (T03_*, puerta publicada)', LL, A, todos
        'a2', 'entre fases frente a T52 G5_005',                    LL, G, tf
        'b',  'a tierra frente a T52 G5_005',                      ~LL, G, tf
        'info', 'a tierra frente a T23 (puerta publicada): cambio esperado por D2', ~LL, A, todos};
C = {}; D = {};
for i = 1:size(defs, 1)
    for val = [1 0]
        R = defs{i, 4}; m = defs{i, 3} & E.valida == val; mc = m & defs{i, 5};
        d1 = mc & ~strcmp(E.outcome, R.outcome); d2 = mc & ~strcmp(E.decl_class, R.decl_class); d3 = mc & ~igual_(E.decl_ms, R.decl_ms);
        dd = d1 | d2 | d3;
        C(end+1, :) = {defs{i, 1}, defs{i, 2}, tern_(val, 'validas', 'no validas (excluidas de los recuentos)'), nnz(m), nnz(mc), nnz(m & ~defs{i, 5}), nnz(d1), nnz(d2 & ~d1), nnz(d3 & ~d1 & ~d2), nnz(dd)}; %#ok<AGROW>
        for k = find(dd)'
            D(end+1, :) = {defs{i, 1}, E.id{k}, E.conjunto{k}, E.terminal{k}, E.tipo{k}, E.valida(k), E.Rf_ohm(k), E.Rtierra_ohm(k), E.V0V1_falta_pct(k), E.V2V1_falta_pct(k), ...
                R.outcome{k}, R.decl_class{k}, R.decl_ms(k), E.outcome{k}, E.decl_class{k}, E.decl_ms(k)}; %#ok<AGROW>
        end
    end
end
C = cell2table(C, 'VariableNames', {'control', 'descripcion', 'poblacion', 'n_eventos', 'n_comparados', 'n_sin_referencia', 'dif_resultado', 'dif_solo_clase', 'dif_solo_latencia', 'dif_total'});
if isempty(D), D = cell(0, 16); end
D = cell2table(D, 'VariableNames', {'control', 'id', 'conjunto', 'terminal', 'tipo', 'valida', 'Rf_ohm', 'Rtierra_ohm', 'V0V1_falta_pct', 'V2V1_falta_pct', 'ref_outcome', 'ref_decl_class', 'ref_decl_ms', 'final_outcome', 'final_decl_class', 'final_decl_ms'});
% lo que T52 ya tenia para esos eventos con la puerta publicada (G0) y con G5_005
T52a = r2_final_emt_leer(fullfile(R2, 'resultados', 'T52', 'emt_eventos.csv')); D.T52_G0_decl_ms = nan(height(D), 1); D.T52_G5_005_decl_ms = nan(height(D), 1);
for k = 1:height(D)
    q = strcmp(T52a.id, D.id{k}) & strcmp(T52a.convertidor, D.terminal{k});
    x = T52a.decl_ms(q & strcmp(T52a.variante, 'G0')); if ~isempty(x), D.T52_G0_decl_ms(k) = x(1); end
    x = T52a.decl_ms(q & strcmp(T52a.variante, 'G5_005')); if ~isempty(x), D.T52_G5_005_decl_ms(k) = x(1); end
end
% secundarios
fprintf('Control secundario: first_label_ms de las faltas entre fases frente a T23: %d diferencias en %d eventos\n', nnz(LL & ~igual_(E.first_label_ms, A.first_label_ms)), nnz(LL));
fprintf('Control secundario: eventos con declaracion en prefalta (n_decl_pre > 0): %d de %d con el clasificador final (T23: %d); solo validos: %d (T23: %d)\n', ...
    nnz(E.n_decl_pre > 0), height(E), nnz(A.n_decl_pre > 0), nnz(E.n_decl_pre > 0 & E.valida == 1), nnz(A.n_decl_pre > 0 & E.valida == 1));
end

% =========================================================================
function U = umbrales_(T, ZBASE)
% Los tres umbrales por tipo de falta y terminal en el barrido de Rf (S2a + S2d).
TIPOS = {'A-G', 'A-B-G', 'A-B'}; TERM = {'Pex', 'DG1'}; rows = {};
for i = 1:numel(TIPOS)
    for j = 1:numel(TERM)
        B = sortrows(T(strcmp(T.tipo, TIPOS{i}) & strcmp(T.terminal, TERM{j}), :), 'Rf_ohm');
        u = umbral_serie_(B.Rf_ohm, B.outcome, B.decl_class);
        v = @(x, col) busca_(B.Rf_ohm, B.(col), x);
        rows(end+1, :) = {TIPOS{i}, TERM{j}, height(B), u.secuencia, ...
            u.ok_cont, u.ok_cont/ZBASE, u.primer_fallo, u.ok_max, ...
            u.decl_cont, u.decl_cont/ZBASE, u.decl_max, u.decl_max/ZBASE, ...
            u.nodecl_primera, u.nodecl_desde, u.nodecl_desde/ZBASE, u.monotono_ok, u.monotono_decl, u.clases_err, ...
            v(u.ok_cont, 'V0V1_falta_pct'), v(u.ok_cont, 'V2V1_falta_pct'), v(u.primer_fallo, 'V0V1_falta_pct'), v(u.primer_fallo, 'V2V1_falta_pct'), ...
            v(u.decl_max, 'V0V1_falta_pct'), v(u.decl_max, 'V2V1_falta_pct'), v(u.nodecl_desde, 'V0V1_falta_pct'), v(u.nodecl_desde, 'V2V1_falta_pct')}; %#ok<AGROW>
    end
end
U = cell2table(rows, 'VariableNames', {'tipo', 'terminal', 'n_casos', 'secuencia', ...
    'Rf_ok_contigua_ohm', 'Rf_ok_contigua_pu', 'Rf_primer_fallo_ohm', 'Rf_ok_max_ohm', ...
    'Rf_decl_contigua_ohm', 'Rf_decl_contigua_pu', 'Rf_decl_max_ohm', 'Rf_decl_max_pu', ...
    'Rf_nodecl_primera_ohm', 'Rf_nodecl_desde_ohm', 'Rf_nodecl_desde_pu', 'monotono_correcta', 'monotono_declaracion', 'clases_erroneas', ...
    'V0V1_en_ok_pct', 'V2V1_en_ok_pct', 'V0V1_en_primer_fallo_pct', 'V2V1_en_primer_fallo_pct', ...
    'V0V1_en_decl_max_pct', 'V2V1_en_decl_max_pct', 'V0V1_en_nodecl_desde_pct', 'V2V1_en_nodecl_desde_pct'});
end

function u = umbral_serie_(Rf, outcome, cls)
% Rf ordenado. ok_cont: mayor Rf tal que todas las Rf menores o iguales son correctas.
% decl_cont: idem con "declara algo". decl_max: mayor Rf con alguna declaracion.
% nodecl_desde: menor Rf tal que ni ella ni ninguna mayor se declara. Sin forzar monotonia.
ok = strcmp(outcome, 'correcta'); dec = ~strcmp(outcome, 'sin_declaracion');
u.ok_cont = contigua_(Rf, ok); u.ok_max = max([Rf(ok); NaN]); u.primer_fallo = min([Rf(~ok); NaN]);
u.decl_cont = contigua_(Rf, dec); u.decl_max = max([Rf(dec); NaN]); u.nodecl_primera = min([Rf(~dec); NaN]);
k = find(dec, 1, 'last');
if isempty(k), u.nodecl_desde = Rf(1); elseif k == numel(Rf), u.nodecl_desde = NaN; else, u.nodecl_desde = Rf(k + 1); end
u.monotono_ok = double(isequaln(u.ok_cont, u.ok_max)); u.monotono_decl = double(isequaln(u.decl_cont, u.decl_max));
s = cell(1, numel(Rf));
for k = 1:numel(Rf)
    if ok(k), c = 'C'; elseif ~dec(k), c = 'N'; else, c = ['E(' cls{k} ')']; end
    s{k} = sprintf('%g:%s', Rf(k), c);
end
u.secuencia = strjoin(s, ' ');
w = dec & ~ok; uc = unique(cls(w), 'stable'); s = cell(1, numel(uc));
for k = 1:numel(uc)
    r = Rf(w & strcmp(cls, uc{k})); s{k} = sprintf('%s (%d casos, %g a %g ohm)', uc{k}, numel(r), min(r), max(r));
end
u.clases_err = strjoin(s, '; ');
end

function r = contigua_(Rf, ok)
k = find(~ok, 1);
if isempty(k), r = max(Rf); elseif k == 1, r = NaN; else, r = Rf(k - 1); end
end

function y = busca_(Rf, col, x)
k = find(Rf == x, 1); if isempty(k), y = NaN; else, y = col(k); end
end

% =========================================================================
function W = barrido_(T, T0)
% Tabla ancha del barrido de Rf: una fila por Rf, seis series (terminal x tipo).
SER = {'Pex', 'A-G'; 'Pex', 'A-B-G'; 'Pex', 'A-B'; 'DG1', 'A-G'; 'DG1', 'A-B-G'; 'DG1', 'A-B'};
Rf = unique(T.Rf_ohm); W = table(Rf, Rf/5.775, 'VariableNames', {'Rf_ohm', 'Rf_pu'});
for s = 1:size(SER, 1)
    nm = [SER{s, 1} '_' strrep(SER{s, 2}, '-', '')];
    idx = find(strcmp(T.terminal, SER{s, 1}) & strcmp(T.tipo, SER{s, 2}));
    [tf, loc] = ismember(Rf, T.Rf_ohm(idx));
    if ~all(tf) || numel(idx) ~= numel(Rf), error('r2_final_emt_resumen:barrido', 'Barrido incompleto en %s', nm); end
    idx = idx(loc);
    W.([nm '_resultado']) = T.outcome(idx); W.([nm '_clase']) = T.decl_class(idx); W.([nm '_decl_ms']) = T.decl_ms(idx);
    W.([nm '_primera_ms']) = T.first_label_ms(idx); W.([nm '_decl_correcta_ms']) = T.correct_decl_ms(idx);
    W.([nm '_V2V1_pct']) = T.V2V1_falta_pct(idx); W.([nm '_V0V1_pct']) = T.V0V1_falta_pct(idx); W.([nm '_v0stat']) = T.full_v0stat_med(idx);
    W.([nm '_T23_resultado']) = T0.outcome(idx); W.([nm '_T23_clase']) = T0.decl_class(idx); W.([nm '_T23_decl_ms']) = T0.decl_ms(idx);
end
end

% =========================================================================
function RC = recuentos_s2_(T, T0)
GR = {'barrido (S2a+S2d)', {'S2a', 'S2d'}; 'incidencia (S2b)', {'S2b'}; 'tierra resistiva (S2c)', {'S2c'}; 'S2 completo', {'S2a', 'S2d', 'S2b', 'S2c'}};
TIPOS = {'A-G', 'A-B-G', 'A-B', 'todos'}; TERM = {'Pex', 'DG1'}; rows = {};
for g = 1:size(GR, 1)
    for j = 1:numel(TERM)
        for i = 1:numel(TIPOS)
            m = ismember(T.conjunto, GR{g, 2}) & strcmp(T.terminal, TERM{j});
            if ~strcmp(TIPOS{i}, 'todos'), m = m & strcmp(T.tipo, TIPOS{i}); end
            if ~any(m), continue, end
            a = cuenta_(T(m, :)); b = cuenta_(T0(m, :));
            rows(end+1, :) = {GR{g, 1}, TERM{j}, TIPOS{i}, nnz(m), a.c, a.e, a.n, a.cls, a.pre, a.wrong_any, a.after, b.c, b.e, b.n, b.cls, b.pre, ...
                nnz(~strcmp(T.outcome(m), T0.outcome(m)) | ~strcmp(T.decl_class(m), T0.decl_class(m)))}; %#ok<AGROW>
        end
    end
end
RC = cell2table(rows, 'VariableNames', {'grupo', 'terminal', 'tipo', 'n_eventos', 'correctas', 'erroneas', 'sin_declaracion', 'clases_erroneas', 'eventos_con_decl_prefalta', ...
    'eventos_con_alguna_decl_erronea_en_falta', 'eventos_con_decl_tras_despeje', 'T23_correctas', 'T23_erroneas', 'T23_sin_declaracion', 'T23_clases_erroneas', 'T23_eventos_con_decl_prefalta', 'eventos_que_cambian'});
end

function a = cuenta_(T)
a.c = nnz(strcmp(T.outcome, 'correcta')); a.e = nnz(strcmp(T.outcome, 'erronea')); a.n = nnz(strcmp(T.outcome, 'sin_declaracion'));
a.pre = nnz(T.n_decl_pre > 0); a.wrong_any = nnz(T.n_decl_wrong > 0); a.after = nnz(T.n_decl_after > 0);
w = strcmp(T.outcome, 'erronea'); uc = unique(T.decl_class(w)); s = cell(1, numel(uc));
for k = 1:numel(uc), s{k} = sprintf('%s (%d)', uc{k}, nnz(w & strcmp(T.decl_class, uc{k}))); end
a.cls = strjoin(s, ', ');
end

% =========================================================================
function INC = incidencia_(T, T0)
% S2b: resultado en los cuatro angulos de incidencia por tipo, terminal y Rf.
m = strcmp(T.conjunto, 'S2b'); T = T(m, :); T0 = T0(m, :);
TIPOS = {'A-G', 'A-B-G', 'A-B'}; TERM = {'Pex', 'DG1'}; rows = {};
for j = 1:numel(TERM)
    for i = 1:numel(TIPOS)
        for rf = unique(T.Rf_ohm)'
            q = strcmp(T.tipo, TIPOS{i}) & strcmp(T.terminal, TERM{j}) & T.Rf_ohm == rf;
            [B, ix] = sortrows(T(q, :), 'ang_inc_pedido_deg'); B0 = T0(q, :); B0 = B0(ix, :);
            ok = strcmp(B.outcome, 'correcta'); er = strcmp(B.outcome, 'erronea'); ok0 = strcmp(B0.outcome, 'correcta');
            rows(end+1, :) = {TERM{j}, TIPOS{i}, rf, rf/5.775, height(B), strjoin(arrayfun(@(x) sprintf('%g', x), B.ang_inc_pedido_deg', 'UniformOutput', false), ' '), ...
                max(abs(mod(B.ang_inc_medido_deg - B.ang_inc_pedido_deg + 180, 360) - 180)), ...
                nnz(ok), nnz(er), nnz(~ok & ~er), strjoin(unique(B.decl_class(er))', '/'), double(numel(unique(strcat(B.outcome, B.decl_class))) == 1), ...
                strjoin(res_(B), ' '), strtrim(sprintf('%.2f ', B.decl_ms)), strtrim(sprintf('%.2f ', B.first_label_ms)), ...
                min([B.decl_ms(ok); NaN]), max([B.decl_ms(ok); NaN]), min([B.first_label_ms(ok); NaN]), max([B.first_label_ms(ok); NaN]), ...
                mean(B.V0V1_falta_pct), mean(B.V2V1_falta_pct), ...
                strjoin(res_(B0), ' '), nnz(ok0), min([B0.decl_ms(ok0); NaN]), max([B0.decl_ms(ok0); NaN])}; %#ok<AGROW>
        end
    end
end
INC = cell2table(rows, 'VariableNames', {'terminal', 'tipo', 'Rf_ohm', 'Rf_pu', 'n_angulos', 'angulos_deg', 'error_angulo_max_deg', ...
    'n_correcta', 'n_erronea', 'n_sin_declaracion', 'clases_erroneas', 'mismo_resultado_en_todos', 'resultado_por_angulo', 'decl_ms_por_angulo', 'primera_ms_por_angulo', ...
    'decl_min_ms', 'decl_max_ms', 'primera_min_ms', 'primera_max_ms', 'V0V1_medio_pct', 'V2V1_medio_pct', 'T23_resultado_por_angulo', 'T23_n_correcta', 'T23_decl_min_ms', 'T23_decl_max_ms'});
end

function s = res_(B)
s = cell(1, height(B));
for k = 1:height(B)
    if strcmp(B.outcome{k}, 'erronea'), s{k} = ['erronea(' B.decl_class{k} ')']; else, s{k} = B.outcome{k}; end
end
end

% =========================================================================
function L = latencias_(T, T0)
% Latencias de los eventos con declaracion correcta, por grupo, terminal y tipo.
GR = {'barrido (S2a+S2d)', {'S2a', 'S2d'}, Inf; 'barrido, Rf <= 0.3 ohm', {'S2a', 'S2d'}, 0.3; 'incidencia (S2b)', {'S2b'}, Inf; 'tierra resistiva (S2c)', {'S2c'}, Inf; 'S2 completo', {'S2a', 'S2d', 'S2b', 'S2c'}, Inf};
TIPOS = {'A-G', 'A-B-G', 'A-B', 'todos'}; TERM = {'Pex', 'DG1'}; rows = {};
for g = 1:size(GR, 1)
    for j = 1:numel(TERM)
        for i = 1:numel(TIPOS)
            m = ismember(T.conjunto, GR{g, 2}) & strcmp(T.terminal, TERM{j}) & T.Rf_ohm <= GR{g, 3};
            if ~strcmp(TIPOS{i}, 'todos'), m = m & strcmp(T.tipo, TIPOS{i}); end
            a = T(m & strcmp(T.outcome, 'correcta'), :); b = T0(m & strcmp(T0.outcome, 'correcta'), :);
            if isempty(a) && isempty(b), continue, end
            rows(end+1, :) = {GR{g, 1}, TERM{j}, TIPOS{i}, nnz(m), height(a), mn_(a.first_label_ms), md_(a.first_label_ms), mx_(a.first_label_ms), mn_(a.decl_ms), md_(a.decl_ms), mx_(a.decl_ms), ...
                height(b), mn_(b.first_label_ms), mx_(b.first_label_ms), mn_(b.decl_ms), md_(b.decl_ms), mx_(b.decl_ms)}; %#ok<AGROW>
        end
    end
end
L = cell2table(rows, 'VariableNames', {'grupo', 'terminal', 'tipo', 'n_eventos', 'n_correctos', 'primera_min_ms', 'primera_mediana_ms', 'primera_max_ms', 'decl_min_ms', 'decl_mediana_ms', 'decl_max_ms', ...
    'T23_n_correctos', 'T23_primera_min_ms', 'T23_primera_max_ms', 'T23_decl_min_ms', 'T23_decl_mediana_ms', 'T23_decl_max_ms'});
end

% =========================================================================
function [NR, NU, NC] = ruido_(N, ZBASE)
% Variante de robustez del barrido de Rf: recuentos, umbrales por realizacion y cambios.
COND = {'cl05', 'snr40', 'cl05_snr40'}; TIPOS = {'A-G', 'A-B-G', 'A-B', 'todos'}; TERM = {'Pex', 'DG1', 'los dos'};
rows = {};
for c = 1:numel(COND)
    X = N(strcmp(N.condicion, COND{c}), :);
    for j = 1:numel(TERM)
        for i = 1:numel(TIPOS)
            m = true(height(X), 1);
            if j < 3, m = strcmp(X.terminal, TERM{j}); end
            if ~strcmp(TIPOS{i}, 'todos'), m = m & strcmp(X.tipo, TIPOS{i}); end
            Y = X(m, :); ev = unique(strcat(Y.id, '|', Y.terminal));
            nev_c = numel(unique(strcat(Y.id(Y.cambia == 1), '|', Y.terminal(Y.cambia == 1))));
            tr = strcat(res_(Y)', {' <- '}, res_base_(Y)'); tr = tr(Y.cambia == 1); ut = unique(tr); s = cell(1, numel(ut));
            for k = 1:numel(ut), s{k} = sprintf('%s (%d)', ut{k}, nnz(strcmp(tr, ut{k}))); end
            okb = strcmp(Y.base_outcome, 'correcta'); okk = okb & strcmp(Y.outcome, 'correcta');
            rows(end+1, :) = {COND{c}, TERM{j}, TIPOS{i}, numel(ev), height(Y), height(Y)/max(1, numel(ev)), nnz(Y.cambia), nev_c, ...
                nnz(okb), nnz(strcmp(Y.outcome, 'correcta')), nnz(strcmp(Y.base_outcome, 'erronea')), nnz(strcmp(Y.outcome, 'erronea')), ...
                nnz(strcmp(Y.base_outcome, 'sin_declaracion')), nnz(strcmp(Y.outcome, 'sin_declaracion')), strjoin(s, ', '), nnz(Y.n_decl_pre > 0), ...
                nnz(okb & ~strcmp(Y.outcome, 'correcta')), max([Y.decl_ms(okk) - Y.base_decl_ms(okk); NaN]), md_(Y.decl_ms(okk) - Y.base_decl_ms(okk)), nnz(Y.decl_ms(okk) - Y.base_decl_ms(okk) > 5), md_(Y.decl_ms(strcmp(Y.outcome, 'correcta'))), max([Y.decl_ms(strcmp(Y.outcome, 'correcta')); NaN])}; %#ok<AGROW>
        end
    end
end
NR = cell2table(rows, 'VariableNames', {'condicion', 'terminal', 'tipo', 'n_eventos', 'n_realizaciones', 'realizaciones_por_evento', 'realizaciones_que_cambian', 'eventos_que_cambian_alguna_vez', ...
    'base_correctas', 'correctas', 'base_erroneas', 'erroneas', 'base_sin_declaracion', 'sin_declaracion', 'transiciones', 'realizaciones_con_decl_prefalta', ...
    'correctas_en_base_que_dejan_de_serlo', 'retraso_max_decl_correcta_ms', 'retraso_mediana_decl_correcta_ms', 'correctas_con_retraso_mayor_de_5_ms', 'decl_correcta_mediana_ms', 'decl_correcta_max_ms'});
% umbrales por realizacion
rows = {}; T3 = {'A-G', 'A-B-G', 'A-B'}; B = N(strcmp(N.condicion, 'base'), :);
for c = 1:numel(COND)
    X = N(strcmp(N.condicion, COND{c}), :);
    rz = unique(strcat(X.canal, '|', arrayfun(@(x) sprintf('%d', x), X.rep, 'UniformOutput', false)), 'stable');
    for j = 1:2
        for i = 1:numel(T3)
            Y = sortrows(B(strcmp(B.tipo, T3{i}) & strcmp(B.terminal, TERM{j}), :), 'Rf_ohm'); u = umbral_serie_(Y.Rf_ohm, Y.outcome, Y.decl_class);
            b = [u.ok_cont, u.decl_max, u.nodecl_desde]; v = nan(numel(rz), 3);
            for r = 1:numel(rz)
                p = strsplit(rz{r}, '|');
                Y = X(strcmp(X.canal, p{1}) & X.rep == str2double(p{2}) & strcmp(X.tipo, T3{i}) & strcmp(X.terminal, TERM{j}), :);
                Y = sortrows(Y, 'Rf_ohm'); u = umbral_serie_(Y.Rf_ohm, Y.outcome, Y.decl_class);
                v(r, :) = [u.ok_cont, u.decl_max, u.nodecl_desde];
            end
            % sin ninguna Rf correcta o declarada: 0; si se declara hasta la mayor Rf del barrido: Inf
            v(isnan(v(:, 1)), 1) = 0; v(isnan(v(:, 2)), 2) = 0; v(isnan(v(:, 3)), 3) = Inf; b0 = b; b(isnan(b0(1:2))) = 0; if isnan(b0(3)), b(3) = Inf; end
            d = [~igual_(v(:, 1), b(1)), ~igual_(v(:, 2), b(2)), ~igual_(v(:, 3), b(3))];
            rows(end+1, :) = {COND{c}, TERM{j}, T3{i}, numel(rz), b(1), min(v(:, 1)), max(v(:, 1)), nnz(d(:, 1)), nnz(v(:, 1) == 0), b(1)/ZBASE, min(v(:, 1))/ZBASE, ...
                b(2), min(v(:, 2)), max(v(:, 2)), nnz(d(:, 2)), b(3), min(v(:, 3)), max(v(:, 3)), nnz(d(:, 3)), double(any(d(:)))}; %#ok<AGROW>
        end
    end
end
NU = cell2table(rows, 'VariableNames', {'condicion', 'terminal', 'tipo', 'n_realizaciones', 'base_Rf_ok_contigua_ohm', 'Rf_ok_contigua_min_ohm', 'Rf_ok_contigua_max_ohm', 'realizaciones_con_ok_distinto', 'realizaciones_sin_ninguna_correcta', ...
    'base_Rf_ok_contigua_pu', 'Rf_ok_contigua_min_pu', 'base_Rf_decl_max_ohm', 'Rf_decl_max_min_ohm', 'Rf_decl_max_max_ohm', 'realizaciones_con_decl_distinto', ...
    'base_Rf_nodecl_desde_ohm', 'Rf_nodecl_desde_min_ohm', 'Rf_nodecl_desde_max_ohm', 'realizaciones_con_nodecl_distinto', 'se_mueve_algun_umbral'});
% cambios por evento y condicion
X = N(~strcmp(N.condicion, 'base') & N.cambia == 1, :); rows = {};
key = unique(strcat(X.condicion, '|', X.id, '|', X.terminal), 'stable');
for k = 1:numel(key)
    p = strsplit(key{k}, '|'); Y = X(strcmp(X.condicion, p{1}) & strcmp(X.id, p{2}) & strcmp(X.terminal, p{3}), :);
    nt = nnz(strcmp(N.condicion, p{1}) & strcmp(N.id, p{2}) & strcmp(N.terminal, p{3}));
    r = res_(Y); ur = unique(r); s = cell(1, numel(ur));
    for q = 1:numel(ur), s{q} = sprintf('%s (%d)', ur{q}, nnz(strcmp(r, ur{q}))); end
    rb = res_base_(Y(1, :));
    rows(end+1, :) = {p{1}, p{2}, p{3}, Y.tipo{1}, Y.Rf_ohm(1), rb{1}, Y.base_decl_ms(1), height(Y), nt, strjoin(s, ', '), strjoin(unique(Y.canal)', '')}; %#ok<AGROW>
end
if isempty(rows), rows = cell(0, 11); end
NC = cell2table(rows, 'VariableNames', {'condicion', 'id', 'terminal', 'tipo', 'Rf_ohm', 'base_resultado', 'base_decl_ms', 'realizaciones_que_cambian', 'realizaciones', 'resultados_distintos', 'canales_afectados'});
if ~isempty(NC), NC = sortrows(NC, {'condicion', 'terminal', 'tipo', 'Rf_ohm'}); end
end

function s = res_base_(B)
s = cell(1, height(B));
for k = 1:height(B)
    if strcmp(B.base_outcome{k}, 'erronea'), s{k} = ['erronea(' B.base_decl_class{k} ')']; else, s{k} = B.base_outcome{k}; end
end
end

% =========================================================================
function [SF, ST] = s1_factor_(E, A)
% S1: recuento por conjunto, factor y terminal (medidas validas), y totales por terminal.
TERM = {'Pex', 'DG1'}; rows = {};
g = unique(strcat(E.conjunto, '|', E.factor), 'stable');
for j = 1:numel(TERM)
    for k = 1:numel(g)
        p = strsplit(g{k}, '|'); m = strcmp(E.conjunto, p{1}) & strcmp(E.factor, p{2}) & strcmp(E.terminal, TERM{j});
        rows(end+1, :) = fila_s1_(p{1}, p{2}, TERM{j}, E(m, :), A(m, :)); %#ok<AGROW>
    end
    for cj = [unique(E.conjunto)', {'S1 completo'}]
        if strcmp(cj{1}, 'S1 completo'), m = strcmp(E.terminal, TERM{j}); else, m = strcmp(E.conjunto, cj{1}) & strcmp(E.terminal, TERM{j}); end
        rows(end+1, :) = fila_s1_(cj{1}, 'todos', TERM{j}, E(m, :), A(m, :)); %#ok<AGROW>
    end
end
SF = cell2table(rows, 'VariableNames', {'conjunto', 'factor', 'terminal', 'n_medidas', 'n_excluidas', 'n_validas', 'correctas', 'par_adyacente', 'erronea_otra', 'sin_declaracion', 'eventos_con_decl_prefalta', ...
    'T23_correctas', 'T23_par_adyacente', 'T23_erronea_otra', 'T23_sin_declaracion', 'eventos_que_cambian'});
rows = {};
for j = 1:numel(TERM)
    B = E(strcmp(E.terminal, TERM{j}) & E.valida == 1, :);
    w = strcmp(B.outcome, 'erronea'); d = strcat(B.tipo(w), {' como '}, B.decl_class(w)); ud = unique(d); s = cell(1, numel(ud));
    for k = 1:numel(ud), s{k} = sprintf('%s (%d)', ud{k}, nnz(strcmp(d, ud{k}))); end
    for tp = [unique(B.tipo)', {'todos'}]
        q = true(height(B), 1); if ~strcmp(tp{1}, 'todos'), q = strcmp(B.tipo, tp{1}); end
        rows(end+1, :) = {TERM{j}, tp{1}, nnz(q), nnz(q & strcmp(B.resultado4, 'correcta')), nnz(q & strcmp(B.resultado4, 'par_adyacente')), nnz(q & strcmp(B.resultado4, 'erronea_otra')), ...
            nnz(q & strcmp(B.resultado4, 'sin_declaracion')), tern_(strcmp(tp{1}, 'todos'), strjoin(s, ', '), '')}; %#ok<AGROW>
    end
end
ST = cell2table(rows, 'VariableNames', {'terminal', 'tipo', 'n_validas', 'correctas', 'par_adyacente', 'erronea_otra', 'sin_declaracion', 'sentido_del_error'});
end

function f = fila_s1_(cj, fac, term, B, B0)
v = B.valida == 1; c = @(T, x) nnz(v & strcmp(T.resultado4, x));
f = {cj, fac, term, height(B), nnz(~v), nnz(v), c(B, 'correcta'), c(B, 'par_adyacente'), c(B, 'erronea_otra'), c(B, 'sin_declaracion'), nnz(v & B.n_decl_pre > 0), ...
    c(B0, 'correcta'), c(B0, 'par_adyacente'), c(B0, 'erronea_otra'), c(B0, 'sin_declaracion'), nnz(v & (~strcmp(B.outcome, B0.outcome) | ~strcmp(B.decl_class, B0.decl_class)))};
end

% =========================================================================
function SA = s1_angulo_(T, T0)
% Error entre el angulo de la elipse EMT y la ley theta = 165 - arg(V2/V1)/2, y frente a T21 (S1t).
rows = {};
for term = {'Pex', 'DG1'}
    for cj = {'S1 completo', 'S1a', 'S1b', 'S1c', 'S1t'}
        m = strcmp(T.terminal, term{1}); if ~strcmp(cj{1}, 'S1 completo'), m = m & strcmp(T.conjunto, cj{1}); end
        B = T(m, :); B0 = T0(m, :); [mx, im] = max(B.err_ley_deg);
        rows(end+1, :) = {'ley 165 - arg(V2/V1)/2', term{1}, cj{1}, 'medidas validas', height(B), median(B.err_ley_deg), mx, B.id{im}, B.V2V1_falta_pct(im), median(B0.err_ley_deg), max(B0.err_ley_deg)}; %#ok<AGROW>
    end
end
for filtro = {'validas con V2/V1 >= 2 %', 'validas'}
    m = strcmp(T.terminal, 'Pex') & strcmp(T.conjunto, 'S1t') & isfinite(T.err_T21_deg);
    if startsWith(filtro{1}, 'validas con'), m = m & T.V2V1_falta_pct >= 2; end
    B = T(m, :); B0 = T0(m, :); [mx, im] = max(B.err_T21_deg);
    rows(end+1, :) = {'prediccion de T21', 'Pex', 'S1t', filtro{1}, height(B), median(B.err_T21_deg), mx, B.id{im}, B.V2V1_falta_pct(im), median(B0.err_T21_deg), max(B0.err_T21_deg)}; %#ok<AGROW>
end
SA = cell2table(rows, 'VariableNames', {'referencia', 'terminal', 'conjunto', 'poblacion', 'n', 'error_mediana_deg', 'error_max_deg', 'id_del_maximo', 'V2V1_del_maximo_pct', 'T23_error_mediana_deg', 'T23_error_max_deg'});
end

% =========================================================================
function diag_a1_(D, R2, OUT)
% Diagnostico de las diferencias del control a1 (faltas entre fases frente a T23):
% etiquetas ventana a ventana de los primeros 25 ms con la puerta publicada y con la final.
% Es la unica parte de este archivo que vuelve a clasificar (dos pasadas por diferencia).
rows = {};
for k = 1:height(D)
    S = load(fullfile(R2, 'emt', 'casos', [D.id{k} '.mat'])); c = S.caso; t_on = c.t_falta; t_off = t_on + c.dur_falta;
    [V, t] = extract_converter_data(S, D.terminal{k}); V = double(V); t = double(t(:));
    a = {'truth', c.tipo, 't_on', t_on, 't_off', t_off, 'f', c.red_f, 't_range', [t_on - 1.0, min(t(end), t_off + 0.12)]};
    o1 = r2_decision_chain(V, t, a{:}); o2 = r2_decision_chain(V, t, r2_final_opts(), a{:});
    m = find(o1.win.t_end > t_on & o1.win.t_end <= t_on + 0.025)';
    for j = m
        rows(end+1, :) = {D.id{k}, D.terminal{k}, 1e3*(o1.win.t_end(j) - t_on), o1.win.label{j}, o1.win.pinmax(j), o1.win.run(j), o2.win.label{j}, o2.win.pinmax(j), o2.win.run(j)}; %#ok<AGROW>
    end
    g = o1.win.t_end > t_on & o1.win.t_end <= t_on + 0.025 & ~o1.win.is_kirchhoff;
    fprintf('Control a1, %s en %s: declaracion a %.2f ms con la puerta publicada y a %.2f ms con la final; en los primeros 25 ms la publicada manda %d ventanas a la rama de tierra (estadistico maximo %.4f, umbral 0.01) y la final %d\n', ...
        D.id{k}, D.terminal{k}, o1.event.decl_ms, o2.event.decl_ms, nnz(g), max([o1.win.pinmax(g); NaN]), nnz(o2.win.t_end > t_on & o2.win.t_end <= t_on + 0.025 & ~o2.win.is_kirchhoff));
end
if isempty(rows), rows = cell(0, 9); end
writetable(cell2table(rows, 'VariableNames', {'id', 'terminal', 't_fin_ventana_ms', 'etiqueta_publicada', 'estadistico_publicada', 'racha_publicada', 'etiqueta_final', 'estadistico_final', 'racha_final'}), ...
    fullfile(OUT, 'control_a1_diagnostico.csv'));
end

function t = igual_(a, b)
t = (isnan(a) & isnan(b)) | abs(a - b) <= 1e-6;
end

function T = renombra_(T, pre)
T.Properties.VariableNames = strcat(pre, T.Properties.VariableNames);
end

function s = tern_(c, a, b)
if c, s = a; else, s = b; end
end

function y = mn_(x), y = min([x(:); NaN]); end
function y = mx_(x), y = max([x(:); NaN]); end
function y = md_(x), x = x(isfinite(x)); if isempty(x), y = NaN; else, y = median(x); end, end
