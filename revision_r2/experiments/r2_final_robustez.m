%% R2_FINAL_ROBUSTEZ  Suite de robustez y cifras de las lineas 1614 y 1927-1953 con el
% clasificador FINAL de la ronda 2 (T90). Copia literal de r2_t03_mapa_cifras.m (T03) con
% tres cambios y nada mas:
%   1) toda llamada a r2_decision_chain pasa por rdc(), que antepone las opciones de
%      r2_final_opts (G5, 0.05); con un escaneo como entrada la cadena no reclasifica y
%      las opciones del clasificador no actuan (el escaneo ya se hizo con ellas);
%   2) scan_dc (DC sobre la falta, con y sin resta LS) clasifica con
%      r2_ga_classifier_param y los parametros finales en vez de ga_fault_classifier_pure;
%   3) salidas en ronda2/resultados/T90/robustez (mismos nombres que en T03).
% La columna "coincide" compara con el valor PUBLICADO (main_final.tex); el valor de T03
% (puerta publicada, convencion de la ronda 2) esta en resultados/T03/mapa_cifras.csv.
here = fileparts(mfilename('fullpath')); RAIZ = fullfile(here, '..');
addpath(fullfile(RAIZ, 'src', 'fault_classifiers')); addpath(fullfile(RAIZ, 'src', 'utils'));
OUT = fullfile(RAIZ, '_revision_TPWRD', 'ronda2', 'resultados', 'T90', 'robustez'); if ~exist(OUT, 'dir'), mkdir(OUT); end
addpath(here);
ROB = fullfile(RAIZ, '_revision_TPWRD', 'robustness');
logf = fullfile(OUT, 'mapa_log.txt'); if exist(logf, 'file'), delete(logf); end
diary(logf); cobj = onCleanup(@() diary('off'));
[FO_, FI_] = r2_final_opts(); fprintf('R2_FINAL_ROBUSTEZ  %s  MATLAB %s\nOpciones finales: %s\nHuella: %s\n', datestr(now, 31), version, FI_.text, FI_.md5);
warning('off', 'all');
[T, info] = r2_sim_cases(); fs = info.fs; f = info.f; T_ON = info.t_on; T_OFF = info.t_off;
W = 20; ST = 5;
MAP = {};   % cifra, linea, publicado, script, M_origen, clasificador_origen, convencion_origen, csv_origen, recalculado, protocolo_recalculo, coincide, nota

d = load(fullfile(info.dir, 'DataGridRL_AB.mat'));
[Vpex, tpex] = extract_converter_data(d, 'Pex');
seg = tpex >= 2.4 & tpex < 2.8; Vh = Vpex(seg, :); th = tpex(seg); th0 = th - th(1);
Npk = max(abs(Vh(:)));
segF = tpex >= 3.00 & tpex < 3.14; Vf = Vpex(segF, :); tf = tpex(segF); tf0 = tf - tf(1);
long = tpex >= 2.0 & tpex < 2.95; Vl = Vpex(long, :); tl = tpex(long);     % tramo sano largo para la cadena completa
G = 'ga_fault_classifier (GA-FuL)';
rate = @(o) o.summary.secured_pct;
chainH = @(V, tt, M) rdc(V, tt, 'gate', 'none', 'M', M);

%% ---- Tabla: control
o10 = chainH(Vh, th, 10); o8 = rdc(o10.scan, [], 'gate', 'none', 'M', 8);
oc = rdc(Vl, tl);     % cadena completa (puerta adaptativa causal)
MAP(end+1, :) = {'Tabla robustez, control: falsa operacion cruda por ventana', 'tab:robustness_suite', '20.8 %', 'robustness_suite_tpwrd_v2.m (solo consola)', 'no aplica', G, 'no aplica', '', ...
    sprintf('%.1f %%', o10.summary.raw_nonhealthy_pct), 'pure, Pex RL_AB t 2.4-2.8 s, 317 ventanas', si(abs(o10.summary.raw_nonhealthy_pct - 20.8) < 0.06), ''};
MAP(end+1, :) = {'Tabla robustez, control: asegurada', 'tab:robustness_suite', '0.0 %', 'robustness_suite_tpwrd_v2.m (solo consola)', '8', G, 'no aplica (tasa)', '', ...
    sprintf('%.1f %% (M=10); %.1f %% (M=8); cadena completa en t 2.0-2.95 s: %.1f %%, %d declaraciones', rate(o10), rate(o8), rate(oc), oc.summary.n_decl), 'pure, mismo tramo; cadena completa con puerta adaptativa causal', si(rate(o10) == 0), 'Publicado con M = 8; con M = 10 el resultado no cambia'};

%% ---- Tabla: DC decayente en sano
tau = 0.030; raw = []; s10 = []; s8 = []; sc = [];
for dc = [0.05 0.10 0.20 0.30]
    Vd = Vh; Vd(:, 1) = Vd(:, 1) + dc*Npk*exp(-th0/tau);
    a = chainH(Vd, th, 10); b = rdc(a.scan, [], 'gate', 'none', 'M', 8);
    Vd2 = Vl; m = tl >= 2.4; Vd2(m, 1) = Vd2(m, 1) + dc*Npk*exp(-(tl(m) - 2.4)/tau);
    c = rdc(Vd2, tl);
    raw(end+1) = a.summary.raw_nonhealthy_pct; s10(end+1) = rate(a); s8(end+1) = rate(b); sc(end+1) = c.summary.n_decl; %#ok<SAGROW>
end
MAP(end+1, :) = {'Tabla robustez, DC decayente 0.05-0.30 pu: cruda', 'tab:robustness_suite', '28-39 %', 'robustness_suite_tpwrd_v2.m', 'no aplica', G, 'no aplica', 'expA2_dc_defensas.csv: 28.4-38.8', ...
    sprintf('%.1f-%.1f %%', min(raw), max(raw)), 'pure, DC en fase A desde el inicio del tramo sano', si(abs(min(raw) - 28.39) < 0.06 && abs(max(raw) - 38.80) < 0.06), ''};
MAP(end+1, :) = {'Tabla robustez, DC decayente: asegurada', 'tab:robustness_suite', '0.0 %', 'robustness_suite_tpwrd_v2.m', '8', G, 'no aplica (tasa)', 'expA2_dc_defensas.csv: 0', ...
    sprintf('max %.1f %% (M=10); max %.1f %% (M=8); cadena completa: %d declaraciones', max(s10), max(s8), sum(sc)), 'pure; cadena completa con DC desde t = 2.4 s', si(max(s10) == 0), 'Publicado con M = 8'};

%% ---- Tabla: energizacion de condensadores
tau_osc = 0.010; t0 = 0.05; touch_r = []; raw = []; s10 = []; s8 = []; sc = []; det = {};
for fo = [400 600 900]
    for am = [0.10 0.20 0.30]
        osc = am*Npk*exp(-max(th0 - t0, 0)/tau_osc).*sin(2*pi*fo*(th0 - t0)).*(th0 >= t0);
        [~, i0] = min(abs(th0 - t0)); w0 = Vh(i0, :)/Npk;
        a = chainH(Vh + osc.*w0, th, 10); b = rdc(a.scan, [], 'gate', 'none', 'M', 8);
        tc = (a.win.i_start + (W - 1)/2)/fs; touch = tc >= t0 - W/fs & tc <= t0 + 3*tau_osc;
        tl0 = tl - 2.4; oscl = am*Npk*exp(-max(tl0 - t0, 0)/tau_osc).*sin(2*pi*fo*(tl0 - t0)).*(tl0 >= t0);
        [~, i0l] = min(abs(tl0 - t0)); c = rdc(Vl + oscl.*(Vl(i0l, :)/Npk), tl);
        touch_r(end+1) = 100*mean(~strcmp(a.win.label_raw(touch), 'A-B-C')); raw(end+1) = a.summary.raw_nonhealthy_pct; %#ok<SAGROW>
        s10(end+1) = rate(a); s8(end+1) = rate(b); sc(end+1) = c.summary.n_decl; %#ok<SAGROW>
        det(end+1, :) = {fo, am, touch_r(end), nnz(touch), raw(end), s8(end), s10(end), c.summary.n_decl, c.summary.decl_classes}; %#ok<SAGROW>
    end
end
writetable(cell2table(det, 'VariableNames', {'fosc_Hz','amp_pu','cruda_ventanas_transitorio_pct','n_ventanas_transitorio','cruda_registro_pct','asegurada_M8_pct','asegurada_M10_pct','cadena_completa_n_decl','cadena_completa_clases'}), fullfile(OUT, 'mapa_detalle_condensadores.csv'));
MAP(end+1, :) = {'Tabla robustez, condensadores: cruda en ventanas del transitorio', 'tab:robustness_suite', '57-89 %', 'robustness_suite_tpwrd.m (v1)', 'no aplica', G, 'centro de ventana para seleccionar las 28 ventanas', 'expB_cap_switching.csv: 57.1-89.3', ...
    sprintf('%.1f-%.1f %%', min(touch_r), max(touch_r)), 'pure, mismas 28 ventanas', si(abs(min(touch_r) - 57.14) < 0.06 && abs(max(touch_r) - 89.29) < 0.06), ''};
MAP(end+1, :) = {'Tabla robustez, condensadores: cruda en el registro completo (nota a)', 'tab:robustness_suite', '24-27 %', 'robustness_suite_tpwrd_v2.m', 'no aplica', G, 'no aplica', 'expB2_cap_persistencia.csv: 23.7-27.1', ...
    sprintf('%.1f-%.1f %%', min(raw), max(raw)), 'pure, 317 ventanas', si(abs(min(raw) - 23.66) < 0.06 && abs(max(raw) - 27.13) < 0.06), ''};
MAP(end+1, :) = {'Tabla robustez, condensadores: asegurada', 'tab:robustness_suite', '0.0 %', 'robustness_checks_residual.m (M = 10, solo consola); la suite v2 con M = 8 da 0.63 % a 900 Hz y 0.1 pu', '10 (la suite v2, 8)', G, 'no aplica (tasa)', 'expB2_cap_persistencia.csv (M = 8): max 0.63 %', ...
    sprintf('max %.2f %% (M=10); max %.2f %% (M=8); cadena completa: %d declaraciones', max(s10), max(s8), sum(sc)), 'pure; cadena completa con transitorio en t = 2.45 s', si(max(s10) == 0), 'El 0.0 % publicado solo se cumple con M = 10; el CSV guardado es de M = 8'};

%% ---- Tabla: secuencia inversa sostenida
v2w = @(k2, npk, tt) k2*npk*[cos(2*pi*f*tt), cos(2*pi*f*tt + 2*pi/3), cos(2*pi*f*tt - 2*pi/3)];
k2s = [0.005 0.01 0.02]; raw = []; s10 = []; s8 = []; sc = []; scn = [];
for k2 = k2s
    a = chainH(Vh + v2w(k2, Npk, th0), th, 10); b = rdc(a.scan, [], 'gate', 'none', 'M', 8);
    c = rdc(Vl + v2w(k2, Npk, tl), tl);
    raw(end+1) = a.summary.raw_nonhealthy_pct; s10(end+1) = rate(a); s8(end+1) = rate(b); sc(end+1) = rate(c); scn(end+1) = c.summary.n_decl; %#ok<SAGROW>
end
MAP(end+1, :) = {'Tabla robustez, V2 = 0.5 %: cruda / asegurada', 'tab:robustness_suite', '62 % / 0.0 %', 'robustness_suite_tpwrd_v2.m (cruda); robustness_checks_v2_m10.m (asegurada)', '10', G, 'no aplica (tasa)', 'expD2_v2_umbral.csv: 62.1 / 0 (M = 8)', ...
    sprintf('%.1f %% / %.1f %% (M=10); %.1f %% (M=8)', raw(1), s10(1), s8(1)), 'pure, 317 ventanas, sin puerta adaptativa', si(abs(raw(1) - 62.15) < 0.06 && s10(1) == 0), ''};
MAP(end+1, :) = {'Tabla robustez, V2 = 1 % y 2 %: cruda', 'tab:robustness_suite', '98-100 %', 'robustness_suite_tpwrd_v2.m', 'no aplica', G, 'no aplica', 'expD2_v2_umbral.csv: 97.5 / 100', ...
    sprintf('%.1f %% / %.1f %%', raw(2), raw(3)), 'pure', si(abs(raw(2) - 97.48) < 0.06 && raw(3) == 100), 'El 97.5 % se publico redondeado a 98'};
MAP(end+1, :) = {'Tabla robustez, V2 = 1 % y 2 %: asegurada', 'tab:robustness_suite', '10.7 % / 55.8 %', 'robustness_checks_v2_m10.m (solo consola)', '10', G, 'no aplica (tasa)', 'expD2_v2_umbral.csv (M = 8): 17.7 / 62.8; d1_parteB (r_c = 0.99, M = 10): 10.7 / 55.8', ...
    sprintf('%.1f %% / %.1f %% (M=10); %.1f %% / %.1f %% (M=8)', s10(2), s10(3), s8(2), s8(3)), 'pure, sin puerta adaptativa', si(abs(s10(2) - 10.73) < 0.06 && abs(s10(3) - 55.84) < 0.06), ''};

%% ---- Trazas con el protocolo de d1 (sano 2.4-2.8 s; falta 2.99-3.18 s con centro de ventana en [3.0, 3.16])
CONVS = info.converters; FF = {'DataGridRL_AB.mat', 'A-B'; 'DataGridRL_BC.mat', 'B-C'; 'DataGridRL_AC.mat', 'C-A'};
H = struct(); F = struct(); FULL = struct();
for c = 1:2
    cv = CONVS{c}; [V, t] = extract_converter_data(d, cv);
    npkc = max(abs(V(t >= 2.4 & t < 2.8, :)), [], 'all');
    for k2 = [0 0.005 0.01 0.02]
        o = rdc(V + v2w(k2, npkc, t), t, 'gate', 'none', 't_range', [2.4 2.8 - 1/fs]);
        H.(hk(cv, k2)) = o.scan;
        o = rdc(V + v2w(k2, npkc, t), t, 'gate', 'none', 't_range', [2.0 2.95]);
        FULL.(hk(cv, k2)) = o.scan;
    end
    for q = 1:3
        dq = load(fullfile(info.dir, FF{q, 1})); [V, t] = extract_converter_data(dq, cv);
        npkq = max(abs(V(t >= 2.4 & t < 2.8, :)), [], 'all');
        for k2 = [0 0.01 0.02]
            o = rdc(V + v2w(k2, npkq, t), t, 'gate', 'none', 't_range', [2.9985 3.1623]);
            F.(fk(cv, FF{q, 2}, k2)) = o.scan;
            o = rdc(V + v2w(k2, npkq, t), t, 'gate', 'none', 't_range', [2.0 3.2]);
            FULL.(fk(cv, FF{q, 2}, k2)) = o.scan;
        end
    end
end
ev = {'t_on', T_ON, 't_off', T_OFF};

% barrido de umbral fijo (parte B de d1)
RC = [0.90 0.92 0.94 0.95 0.96 0.97 0.98 0.99 0.995]; detB = {};
for rc = RC
    fo = [];
    for c = 1:2, for k2 = [0 0.005 0.01 0.02]
        o = rdc(H.(hk(CONVS{c}, k2)), [], 'gate', 'fixed', 'gate_rc', rc); fo(end+1) = rate(o); %#ok<SAGROW>
    end, end
    nok = 0;
    for c = 1:2, for q = 1:3, for k2 = [0 0.01]
        o = rdc(F.(fk(CONVS{c}, FF{q, 2}, k2)), [], 'gate', 'fixed', 'gate_rc', rc, 'truth', FF{q, 2}, ev{:});
        nok = nok + strcmp(o.event.outcome, 'correcta');
    end, end, end
    detB(end+1, :) = {rc, max(fo), nok}; %#ok<SAGROW>
end
TB = cell2table(detB, 'VariableNames', {'r_c', 'falsa_asegurada_max_pct', 'faltas_correctas_de_12'});
writetable(TB, fullfile(OUT, 'mapa_detalle_gate_fijo.csv')); disp(TB);
rc_seg = max(TB.r_c(TB.falsa_asegurada_max_pct == 0)); rc_dep = min(TB.r_c(TB.faltas_correctas_de_12 == 12));
MAP(end+1, :) = {'Umbral fijo: seguridad hasta V2 = 2 % exige r_c <= 0.96; declarar las 12 faltas exige r_c >= 0.97', '1949', '0.96 / 0.97', 'd1_roundness_gate_design.m, parte B', '10', G, 'fin de ventana', 'd1_parteB_gate_fijo.csv', ...
    sprintf('%.2f / %.2f', rc_seg, rc_dep), 'pure, protocolo d1, puerta fija de la cadena', si(rc_seg == 0.96 && rc_dep == 0.97), 'Las 12 faltas son 3 faltas RL x 2 convertidores x 2 niveles de V2 (0 y 1 %), no los 12 archivos'};

% puerta adaptativa (parte C de d1), Delta = 0.01
fo = []; lat = []; nok = 0; detC = {};
for c = 1:2, for k2 = [0 0.005 0.01 0.02]
    S = H.(hk(CONVS{c}, k2)); o = rdc(S, [], 'amb_mode', 'fixed', 'r_amb', median(S.roundness, 'omitnan'));
    oc = rdc(FULL.(hk(CONVS{c}, k2)), []);
    fo(end+1) = rate(o); %#ok<SAGROW>
    detC(end+1, :) = {CONVS{c}, 'sano', 100*k2, 'protocolo d1 (r_amb fuera de linea)', '', NaN, rate(o), o.summary.n_decl}; %#ok<SAGROW>
    detC(end+1, :) = {CONVS{c}, 'sano', 100*k2, 'cadena causal, t 2.0-2.95 s', '', NaN, rate(oc), oc.summary.n_decl}; %#ok<SAGROW>
end, end
latc = []; nokc = 0;
for c = 1:2, for q = 1:3, for k2 = [0 0.01 0.02]
    cv = CONVS{c}; ra = median(H.(hk(cv, k2)).roundness, 'omitnan');   % sano del archivo A-B, como d1
    o = rdc(F.(fk(cv, FF{q, 2}, k2)), [], 'amb_mode', 'fixed', 'r_amb', ra, 'truth', FF{q, 2}, ev{:});
    oc = rdc(FULL.(fk(cv, FF{q, 2}, k2)), [], 'truth', FF{q, 2}, ev{:});
    if k2 < 0.015
        nok = nok + strcmp(o.event.outcome, 'correcta'); lat(end+1) = o.event.decl_ms; %#ok<SAGROW>
        nokc = nokc + strcmp(oc.event.outcome, 'correcta'); latc(end+1) = oc.event.decl_ms; %#ok<SAGROW>
    end
    detC(end+1, :) = {cv, FF{q, 2}, 100*k2, 'protocolo d1 (r_amb fuera de linea)', o.event.outcome, o.event.decl_ms, o.event.retention_pct, o.event.n_decl_wrong}; %#ok<SAGROW>
    detC(end+1, :) = {cv, FF{q, 2}, 100*k2, 'cadena causal, registro t 2.0-3.2 s', oc.event.outcome, oc.event.decl_ms, oc.event.retention_pct, oc.event.n_decl_wrong + oc.event.n_decl_pre}; %#ok<SAGROW>
end, end, end
TCd = cell2table(detC, 'VariableNames', {'convertidor','caso','v2_pct','protocolo','resultado','decl_ms','asegurada_o_retencion_pct','n_decl_o_erroneas'});
writetable(TCd, fullfile(OUT, 'mapa_detalle_gate_adaptativo.csv')); disp(TCd);
sl = sort(lat);
MAP(end+1, :) = {'Puerta adaptativa: falsas aseguradas hasta V2 = 2 %', '1949 y tab:robustness_suite', '0.0 %', 'd1_roundness_gate_design.m, parte C', '10', G, 'no aplica (tasa)', 'd1_parteC_gate_adaptativo.csv: 0 en 8 trazas', ...
    sprintf('max %.1f %% (protocolo d1); cadena causal: %d declaraciones en 8 trazas', max(fo), sum(cell2mat(detC(strcmp(detC(:, 2), 'sano') & startsWith(detC(:, 4), 'cadena'), 8)))), 'pure; d1 usa r_amb fuera de linea (mediana de todo el tramo sano), no las K ventanas precedentes', si(max(fo) == 0), 'El texto dice «preceding K = 80 healthy windows»; el script no lo implementa asi'};
MAP(end+1, :) = {'Puerta adaptativa: faltas declaradas', '1949', '12/12', 'd1_roundness_gate_design.m, parte C', '10', G, 'fin de ventana', 'd1_parteC_gate_adaptativo.csv', ...
    sprintf('%d/12 (protocolo d1); %d/12 (cadena causal)', nok, nokc), 'pure', si(nok == 12 && nokc == 12), '3 faltas RL x 2 convertidores x V2 de 0 y 1 %; r_amb de las faltas tomado del tramo sano del archivo A-B'};
MAP(end+1, :) = {'Puerta adaptativa: latencias tipicas / peor caso', '1949', '15-19 ms / 56 ms', 'd1_roundness_gate_design.m, parte C', '10', G, 'fin de ventana; solo ventanas con centro >= t_on', 'd1_parteC_gate_adaptativo.csv: 14.75-18.5 y 56', ...
    sprintf('%.2f-%.2f ms / %.2f ms (protocolo d1); %.2f-%.2f ms / %.2f ms (cadena causal, registro completo)', sl(1), sl(end-1), sl(end), min(latc), max(latc(latc < 50)), max(latc)), 'pure', si(abs(sl(1) - 14.75) < 1e-6 && abs(sl(end) - 56) < 1e-6), 'Con el registro completo la racha puede empezar en ventanas que acaban 1 a 3.5 ms tras t_on y d1 las descarta; por eso la cadena da hasta 1.25 ms menos en algunos casos'};
r14 = [median(F.(fk('Pex', 'A-B', 0)).roundness), median(F.(fk('Pex', 'B-C', 0)).roundness), median(F.(fk('Pex', 'C-A', 0)).roundness)];
MAP(end+1, :) = {'Redondez de la falta en el terminal GFL', '1949', 'aprox. 0.14', 'd1_roundness_gate_design.m, parte A', 'no aplica', G, 'no aplica', 'd1_parteA_roundness.csv: 0.136', ...
    sprintf('%.3f / %.3f / %.3f (A-B, B-C, C-A)', r14), 'pure, mediana de las ventanas de falta', si(all(abs(r14 - 0.136) < 0.005)), ''};
% DG1 con 2 % ambiente
s = '';
for q = 1:3
    oc = rdc(FULL.(fk('DG1', FF{q, 2}, 0.02)), [], 'truth', FF{q, 2}, ev{:});
    s = [s, sprintf('%s: %s; ', FF{q, 2}, oc.event.outcome)]; %#ok<AGROW>
end
MAP(end+1, :) = {'GFM con 2 % ambiente: algun tipo de falta puede quedar bajo el suelo', '1949', 'cualitativo', 'd1_roundness_gate_design.m, anexo (solo consola)', '10', G, 'fin de ventana', '', s, 'pure, cadena causal, registro completo', 'cualitativo', ''};

%% ---- Suelo de detectabilidad (parte D de d1)
tsyn = (0:4399)'/fs; vb = [cos(2*pi*f*tsyn), cos(2*pi*f*tsyn - 2*pi/3), cos(2*pi*f*tsyn + 2*pi/3)]; T0S = 0.5;
KSEV = [0.01 0.02 0.03 0.04 0.05 0.10 0.15 0.20 0.30]; detD = {}; kmin = nan(2, 2);
for iv = 1:2
    k2 = (iv - 1)*0.01;
    for k = KSEV
        V = vb; dv = (k/2)*(V(:, 1) - V(:, 2)); m = tsyn >= T0S; V(m, 1) = V(m, 1) - dv(m); V(m, 2) = V(m, 2) + dv(m);
        V = V + v2w(k2, 1, tsyn);
        oo = rdc(V, tsyn, 'amb_mode', 'oracle', 'truth', 'A-B', 't_on', T0S);
        oc = rdc(oo.scan, [], 'truth', 'A-B', 't_on', T0S);
        if strcmp(oo.event.outcome, 'correcta') && isnan(kmin(iv, 1)), kmin(iv, 1) = k; end
        if strcmp(oc.event.outcome, 'correcta') && isnan(kmin(iv, 2)), kmin(iv, 2) = k; end
        detD(end+1, :) = {k, 100*k2, oo.event.outcome, oo.event.decl_ms, oc.event.outcome, oc.event.decl_ms}; %#ok<SAGROW>
    end
end
writetable(cell2table(detD, 'VariableNames', {'k_severidad','v2_pct','oraculo_resultado','oraculo_decl_ms','causal_resultado','causal_decl_ms'}), fullfile(OUT, 'mapa_detalle_detectabilidad.csv'));
MAP(end+1, :) = {'Suelo de detectabilidad: severidad minima declarada (fondo equilibrado / 1 % ambiente)', '1949', 'k = 0.02 / 0.04', 'd1_roundness_gate_design.m, parte D', '10', G, 'fin de ventana; r_amb de las K ventanas previas al inicio (conocido)', 'd1_parteD_detectabilidad.csv', ...
    sprintf('k = %.2f / %.2f (oraculo, como d1); k = %.2f / %.2f (causal)', kmin(1, 1), kmin(2, 1), kmin(1, 2), kmin(2, 2)), 'pure, faltas A-B sinteticas', si(kmin(1, 1) == 0.02 && kmin(2, 1) == 0.04), ''};

%% ---- Falta real A-B: primera correcta y declaracion (linea 1949)
of = rdc(Vf, tf0, 'gate', 'none');      % protocolo original: segmento desde 3.00 s
lab = of.win.label_raw; tcent = (of.win.i_start + (W - 1)/2)/fs;
i1 = find(strcmp(lab, 'A-B'), 1); c8 = NaN; c10 = NaN;
for i = 1:numel(lab)
    if i >= 8 && all(strcmp(lab(i-7:i), 'A-B')) && isnan(c8), c8 = tcent(i); end
    if i >= 10 && all(strcmp(lab(i-9:i), 'A-B')) && isnan(c10), c10 = tcent(i); end
end
ofull = rdc(Vpex, tpex, 'gate', 'none', 'truth', 'A-B', ev{:}, 't_range', [2.0 3.2]);
ofa = rdc(ofull.scan, [], 'truth', 'A-B', ev{:});
MAP(end+1, :) = {'Falta real A-B (GFL): primera ventana correcta', '1949', '2.6 ms', 'robustness_suite_tpwrd_v2.m (solo consola)', 'no aplica', G, 'CENTRO de ventana; rejilla que empieza en t_on', '', ...
    sprintf('%.3f ms al centro (protocolo original); %.2f ms al final de ventana en esa rejilla; %.2f ms con el registro completo', 1e3*tcent(i1), 1e3*of.win.t_end(i1), ofull.event.first_label_ms), 'pure; cadena sobre el registro completo, fin de ventana', si(abs(1e3*tcent(i1) - 2.625) < 1e-6), 'Con marca al final de ventana la primera ventana integra en la falta acaba a 4.75 ms; antes hay una ventana mixta que ya da A-B'};
MAP(end+1, :) = {'Falta real A-B (GFL): declaracion asegurada', '1949', '13.9 ms', 'robustness_checks_dc_dependability.m (fila DC = 0, solo consola)', '10', G, 'CENTRO de ventana', '', ...
    sprintf('%.3f ms al centro con M=10 (%.3f con M=8, que es lo que imprime la suite v2); %.2f ms al final de ventana (cadena, sin puerta); %.2f ms (cadena completa)', 1e3*c10, 1e3*c8, ofull.event.decl_ms, ofa.event.decl_ms), 'pure', si(abs(1e3*c10 - 13.875) < 1e-6), 'Valor a citar con la convencion de la ronda 2: 16.0 ms'};
MAP(end+1, :) = {'Tiempo de confirmacion que anade M = 10', '1927', '12.5 ms', 'aritmetica del texto', '10', 'no aplica', 'no aplica', '', '11.25 ms', 'nueve intervalos de 1.25 ms entre la primera y la decima ventana', 'no', 'Para T33: (M-1) x paso / fs = 9 x 5 / 4000 s'};

%% ---- DC sobre la propia falta (linea 1949)
detDC = {}; pub = []; rec = [];
for dc = [0 0.05 0.10 0.20 0.30]
    % protocolo original: segmento de falta, DC desde su primera muestra, historia solo dentro del segmento
    Vx = Vf; Vx(:, 1) = Vx(:, 1) + dc*Npk*exp(-tf0/0.030);
    Sr = scan_dc(Vx, tf0, fs, f, W, ST, false); Sd = scan_dc(Vx, tf0, fs, f, W, ST, true);
    tcs = (Sr.i_start + (W - 1)/2)/fs;
    a = fire_first(Sr.label_raw, tcs, 10); b = fire_first(Sd.label_raw, tcs, 10);
    % cadena: registro completo, DC desde t_on, historia real de prefalta, fin de ventana
    Vy = Vpex; m = tpex >= T_ON; Vy(m, 1) = Vy(m, 1) + dc*Npk*exp(-(tpex(m) - T_ON)/0.030);
    mm = tpex >= 2.0 & tpex <= 3.2;
    Sy = scan_dc(Vy(mm, :), tpex(mm), fs, f, W, ST, true); Sy0 = scan_dc(Vy(mm, :), tpex(mm), fs, f, W, ST, false);
    oy = rdc(Sy, [], 'gate', 'none', 'truth', 'A-B', ev{:}); oy0 = rdc(Sy0, [], 'gate', 'none', 'truth', 'A-B', ev{:});
    oya = rdc(Sy, [], 'truth', 'A-B', ev{:}); oy0a = rdc(Sy0, [], 'truth', 'A-B', ev{:});
    detDC(end+1, :) = {dc, 1e3*a, 1e3*b, oy0.event.outcome, oy0.event.decl_ms, oy.event.outcome, oy.event.decl_ms, oy0a.event.outcome, oy0a.event.decl_ms, oya.event.outcome, oya.event.decl_ms, oya.event.n_decl_wrong}; %#ok<SAGROW>
    if dc > 0, pub(end+1) = 1e3*b; rec(end+1) = oya.event.decl_ms; end %#ok<SAGROW>
end
TDC = cell2table(detDC, 'VariableNames', {'Idc_pu','orig_cruda_centro_ms','orig_LS_centro_ms','cadena_singate_cruda','cadena_singate_cruda_ms','cadena_singate_LS','cadena_singate_LS_ms','cadena_cruda','cadena_cruda_ms','cadena_LS','cadena_LS_ms','cadena_LS_n_erroneas'});
writetable(TDC, fullfile(OUT, 'mapa_detalle_dc_falta.csv')); disp(TDC);
MAP(end+1, :) = {'DC sobre la falta: retraso de la declaracion asegurada (0.05-0.30 pu, con resta LS)', '1949', '62-112 ms', 'robustness_checks_dc_dependability.m (solo consola)', '10', G, 'CENTRO de ventana; DC e historia solo desde t_on', '', ...
    sprintf('%.1f-%.1f ms al centro (protocolo original); %.1f-%.1f ms al final de ventana (cadena completa, historia de prefalta)', min(pub), max(pub), min(rec), max(rec)), 'pure; ver mapa_detalle_dc_falta.csv', si(abs(min(pub) - 62) < 1 && abs(max(pub) - 112) < 1), ''};

%% ---- Errores de instrumentacion (linea 1951)
detI = {};
for g = [0.005 0.01 0.03]
    Vg = Vh; Vg(:, 1) = Vg(:, 1)*(1 + g); a = chainH(Vg, th, 10);
    Vgf = Vf; Vgf(:, 1) = Vgf(:, 1)*(1 + g); b = rdc(Vgf, tf0, 'gate', 'none'); ret = 100*mean(strcmp(b.win.label_raw, 'A-B'));
    Ve = r2_apply_meas_error(Vpex, tpex, 'gain', [g 0 0]);
    e0 = rdc(Ve, tpex, 'gate', 'none', 'truth', 'A-B', ev{:}, 't_range', [2.0 3.2]); e1 = rdc(e0.scan, [], 'truth', 'A-B', ev{:});
    detI(end+1, :) = {'ganancia_pct', 100*g, a.summary.raw_nonhealthy_pct, rate(a), ret, e0.event.outcome, e0.event.decl_class, e0.event.decl_ms, e0.event.max_run_correct, e0.event.retention_pct, e0.event.n_decl_pre, e1.event.outcome, e1.event.decl_class, e1.event.decl_ms, e1.event.n_decl_pre, e1.event.pre_classes}; %#ok<SAGROW>
end
for dg = [0.5 1.0]
    dt = dg/360/f;
    Vp = Vh; Vp(:, 1) = interp1(th0, Vh(:, 1), th0 - dt, 'spline', 'extrap'); a = chainH(Vp, th, 10);
    Vpf = Vf; Vpf(:, 1) = interp1(tf0, Vf(:, 1), tf0 - dt, 'spline', 'extrap'); b = rdc(Vpf, tf0, 'gate', 'none'); ret = 100*mean(strcmp(b.win.label_raw, 'A-B'));
    Ve = r2_apply_meas_error(Vpex, tpex, 'phase_deg', [dg 0 0]);
    e0 = rdc(Ve, tpex, 'gate', 'none', 'truth', 'A-B', ev{:}, 't_range', [2.0 3.2]); e1 = rdc(e0.scan, [], 'truth', 'A-B', ev{:});
    detI(end+1, :) = {'fase_grados', dg, a.summary.raw_nonhealthy_pct, rate(a), ret, e0.event.outcome, e0.event.decl_class, e0.event.decl_ms, e0.event.max_run_correct, e0.event.retention_pct, e0.event.n_decl_pre, e1.event.outcome, e1.event.decl_class, e1.event.decl_ms, e1.event.n_decl_pre, e1.event.pre_classes}; %#ok<SAGROW>
end
TI = cell2table(detI, 'VariableNames', {'error','valor','sano_cruda_pct','sano_asegurada_pct','falta_retencion_orig_pct','ev_singate','ev_singate_clase','ev_singate_ms','ev_racha_max','ev_retencion_pct','ev_singate_decl_pre','ev_cadena','ev_cadena_clase','ev_cadena_ms','ev_cadena_decl_pre','ev_cadena_clases_pre'});
writetable(TI, fullfile(OUT, 'mapa_detalle_instrumentacion.csv')); disp(TI);
MAP(end+1, :) = {'Error de ganancia en una fase hasta 1 %: falsas aseguradas / retencion de la falta', '1951', '0 % / 100 %', 'robustness_sensor_errors_r12.m', '10', G, 'no aplica (tasas por ventana)', 'r12_gain_error.csv: 0 / 100', ...
    sprintf('%.1f %% / %.1f %% (1 %%); por evento (cadena completa): %s a %.2f ms', TI.sano_asegurada_pct(2), TI.falta_retencion_orig_pct(2), TI.ev_cadena{2}, TI.ev_cadena_ms(2)), 'pure; evento con el error aplicado al registro completo', si(TI.sano_asegurada_pct(2) == 0 && TI.falta_retencion_orig_pct(2) == 100), ''};
MAP(end+1, :) = {'Error de ganancia del 3 %: «rompe la puerta de excentricidad»', '1951', '97.2 % falsas aseguradas; retencion 49.5 % (en el CSV)', 'robustness_sensor_errors_r12.m', '10', G, 'no aplica (tasas por ventana)', 'r12_gain_error.csv: 100 / 97.2 / 49.5', ...
    sprintf('%.1f %% cruda / %.1f %% asegurada / retencion %.1f %%; por evento: %s (%s), %d declaraciones en prefalta (%s)', TI.sano_cruda_pct(3), TI.sano_asegurada_pct(3), TI.falta_retencion_orig_pct(3), TI.ev_cadena{3}, TI.ev_cadena_clase{3}, TI.ev_cadena_decl_pre(3), TI.ev_cadena_clases_pre{3}), 'pure', si(abs(TI.sano_asegurada_pct(3) - 97.16) < 0.06), 'Lo que falla es la puerta de secuencia cero (etiquetas de tierra), no la de excentricidad: ver clases declaradas'};
MAP(end+1, :) = {'Error de fase 0.5 y 1 grado en una fase: falsas aseguradas', '1951', '0-1.3 %', 'robustness_sensor_errors_r12.m', '10', G, 'no aplica (tasa)', 'r12_phase_error.csv: 0 / 1.26', ...
    sprintf('%.2f %% / %.2f %%', TI.sano_asegurada_pct(4), TI.sano_asegurada_pct(5)), 'pure, 317 ventanas sanas', si(TI.sano_asegurada_pct(4) == 0 && abs(TI.sano_asegurada_pct(5) - 1.26) < 0.06), ''};
MAP(end+1, :) = {'Error de fase 0.5 y 1 grado: retencion por ventana de la falta real', '1951', '64 % / 42 %', 'robustness_sensor_errors_r12.m', 'no aplica', G, 'no aplica (retencion por ventana; no mide declaracion)', 'r12_phase_error.csv: 64.22 / 42.20', ...
    sprintf('%.2f %% / %.2f %%; por evento con M = 10: 0.5 grados %s (racha max %d), 1 grado %s (racha max %d)', TI.falta_retencion_orig_pct(4), TI.falta_retencion_orig_pct(5), TI.ev_cadena{4}, TI.ev_racha_max(4), TI.ev_cadena{5}, TI.ev_racha_max(5)), 'pure; evento con el desfase aplicado al registro completo', si(abs(TI.falta_retencion_orig_pct(4) - 64.22) < 0.06 && abs(TI.falta_retencion_orig_pct(5) - 42.20) < 0.06), 'La frase «preserves security ... is therefore sufficient» no mide la declaracion: ver T02'};

%% ---- Punto de onda (linea 1951)
Tt = 1.0; tt = (0:1/fs:Tt - 1/fs)'; ea = cos(2*pi*f*tt); eb = cos(2*pi*f*tt - 2*pi/3); ec = cos(2*pi*f*tt + 2*pi/3);
detP = {};
for k = [0.3 0.6]
    for phi = 0:30:330
        t0p = 0.5 + (phi/360)/f; on = tt >= t0p; dd = (ea - eb)/2; V = [ea - k*dd.*on, eb + k*dd.*on, ec];
        i0 = find(tt >= t0p, 1) - round(fs/f);
        a = rdc(V(i0:end, :), tt(i0:end), 'gate', 'none', 'truth', 'A-B', 't_on', t0p);   % rejilla del script original
        b = rdc(V, tt, 'truth', 'A-B', 't_on', t0p);                                       % cadena completa
        detP(end+1, :) = {k, phi, a.event.first_label_ms, a.event.correct_decl_ms, a.event.outcome, b.event.first_label_ms, b.event.decl_ms, b.event.outcome}; %#ok<SAGROW>
    end
end
TPw = cell2table(detP, 'VariableNames', {'k','phi_deg','orig_primera_ms','orig_asegurada_ms','orig_resultado','cadena_primera_ms','cadena_decl_ms','cadena_resultado'});
writetable(TPw, fullfile(OUT, 'mapa_detalle_punto_de_onda.csv'));
MAP(end+1, :) = {'Punto de onda: primera correcta / asegurada / fallos', '1951', '1.0-4.9 ms / 12.3-16.2 ms / 0', 'pow_inception_sweep_r13.m', '10', G, 'fin de ventana', 'r13_pow_sweep.csv', ...
    sprintf('%.2f-%.2f ms / %.2f-%.2f ms / %d fallos (rejilla original); cadena completa: %.2f-%.2f ms / %.2f-%.2f ms / %d sin declaracion correcta', min(TPw.orig_primera_ms), max(TPw.orig_primera_ms), min(TPw.orig_asegurada_ms), max(TPw.orig_asegurada_ms), nnz(~strcmp(TPw.orig_resultado, 'correcta')), ...
    min(TPw.cadena_primera_ms), max(TPw.cadena_primera_ms), min(TPw.cadena_decl_ms), max(TPw.cadena_decl_ms), nnz(~strcmp(TPw.cadena_resultado, 'correcta'))), 'pure, 24 casos sinteticos', ...
    si(abs(min(TPw.orig_primera_ms) - 1.0) < 0.01 && abs(max(TPw.orig_primera_ms) - 4.917) < 0.01 && abs(min(TPw.orig_asegurada_ms) - 12.25) < 0.01 && abs(max(TPw.orig_asegurada_ms) - 16.167) < 0.01), ''};

%% ---- SLG en sistema aislado (linea 1953)
tC = (0:1/fs:0.4 - 1/fs)'; ea = cos(2*pi*f*tC); eb = cos(2*pi*f*tC - 2*pi/3); ec = cos(2*pi*f*tC + 2*pi/3); acc = [];
for k = [0.25 0.5 0.75 1.0]
    a = rdc([ea, eb, ec] - k*ea, tC, 'gate', 'none'); acc(end+1) = 100*mean(strcmp(a.win.label_raw, 'A-G')); %#ok<SAGROW>
end
MAP(end+1, :) = {'SLG en sistema aislado: A-G correcto por ventana, k = 0.25 a 1', '1953', '100 %', 'robustness_suite_tpwrd.m (v1), experimento C', 'no aplica', G, 'no aplica (regimen permanente sintetico)', 'expC_slg_no_aterrado.csv: 100', ...
    sprintf('%.1f / %.1f / %.1f / %.1f %%', acc), 'pure', si(all(acc == 100)), 'v0 = k pu y la corriente residual 3 w C0 V0 son expresiones cerradas del script, no salidas del clasificador'};

%% ---- Laboratorio (linea 1614)
labdir = fullfile(RAIZ, 'data', 'lab', 'Fault_emulator_Outkontaktor'); files = dir(fullfile(labdir, '*.mat'));
fsL = 10000; WL = 50; STL = 12; CH = 59:61; PERM = [3 2 1]; cyc = fsL/f;
emap = containers.Map({'ab', 'bc', 'ac'}, {'A-B', 'B-C', 'C-A'}); detL = {}; totw = 0; totf = 0; totfa = 0;
for k = 1:numel(files)
    fn = files(k).name; pref = regexp(fn, '^[a-z]+', 'match', 'once');
    if isempty(pref) || ~isKey(emap, pref), continue; end
    expected = emap(pref); dl = load(fullfile(labdir, fn));
    R = dl.results; t = R(1, :)'; t = t - t(1); V = R(CH, :)'; V = V(:, PERM);
    % localizacion e incepcion, copiadas de lab_latency_r3m1.m
    nw = floor(numel(t)/cyc); dev = zeros(nw, 1);
    for i = 1:nw, sg = (i-1)*cyc+1:i*cyc; r3 = sqrt(mean(V(sg, :).^2, 1)); dev(i) = max(r3) - min(r3); end
    base = median(dev); anom = find(dev > 5*base + 0.02, 1);
    i_region = (anom - 2)*cyc; idx_fit = i_region - 2*cyc:i_region - 1;
    tq = t(idx_fit); X = [cos(2*pi*f*tq), sin(2*pi*f*tq), ones(numel(tq), 1)]; coefs = X\V(idx_fit, :);
    idx_test = i_region:min(i_region + 3*cyc, numel(t));
    Xt = [cos(2*pi*f*t(idx_test)), sin(2*pi*f*t(idx_test)), ones(numel(idx_test), 1)];
    resid = abs(V(idx_test, :) - Xt*coefs); Vpk = max(abs(V(idx_fit, :)), [], 'all');
    hit = all(movmin(max(resid, [], 2), [0 4]) > 0.05*Vpk, 2); ii = find(hit, 1); t_incep = t(idx_test(ii));
    i0 = find(t >= t_incep, 1) - cyc; i1 = min(find(t >= t_incep, 1) + 12*cyc, numel(t) - WL);
    last = i0 + STL*floor((i1 - i0)/STL) + WL - 1;
    a = rdc(V(i0:last, :), t(i0:last), 'fs', fsL, 'W', WL, 'stride', STL, 'gate', 'none', 'truth', expected, 't_on', t_incep);
    b = rdc(V(1:last, :), t(1:last), 'fs', fsL, 'W', WL, 'stride', STL, 'truth', expected, 't_on', t_incep);   % cadena completa desde el inicio del registro
    % seguridad en prefalta, como lab_prefault_security.m: 1 s hasta 3 ciclos antes de la falta
    iEnd = max(1, (anom - 3)*cyc); iStart = max(1, iEnd - 10000);
    p = rdc(V(iStart:iEnd, :), t(iStart:iEnd), 'fs', fsL, 'W', WL, 'stride', STL, 'gate', 'none');
    pa = rdc(V(1:iEnd, :), t(1:iEnd), 'fs', fsL, 'W', WL, 'stride', STL);
    totw = totw + p.summary.n_win; totf = totf + nnz(p.win.fire); totfa = totfa + pa.summary.n_decl;
    detL(end+1, :) = {fn, expected, t_incep, a.event.first_label_ms, a.event.correct_decl_ms, a.event.outcome, a.event.decl_class, a.event.decl_ms, b.event.outcome, b.event.decl_class, b.event.decl_ms, b.event.n_decl_pre, p.summary.n_win, nnz(p.win.fire), pa.summary.n_decl}; %#ok<SAGROW>
end
if isempty(detL)
    fprintf('AVISO: no hay registros de laboratorio en %s; se omite el bloque de laboratorio.\n', labdir);
else
TL = cell2table(detL, 'VariableNames', {'registro','truth','t_incep_s','primera_ms','asegurada_correcta_ms','singate_resultado','singate_clase','singate_decl_ms','cadena_resultado','cadena_clase','cadena_decl_ms','cadena_decl_pre','prefalta_ventanas','prefalta_ventanas_aseguradas','prefalta_cadena_n_decl'});
writetable(TL, fullfile(OUT, 'mapa_detalle_laboratorio.csv')); disp(TL);
l1 = TL.primera_ms; lM = TL.asegurada_correcta_ms;
MAP(end+1, :) = {'Laboratorio: primera ventana correcta, mediana (p90; rango)', '1614', '4.1 ms (5.3; 0.5-5.3)', 'lab_latency_r3m1.m (solo consola)', 'no aplica', G, 'fin de ventana; 10 kHz, W = 50, paso 12', '', ...
    sprintf('%.2f ms (%.2f; %.2f-%.2f)', median(l1, 'omitnan'), pctl(l1, 90), min(l1), max(l1)), 'pure, misma rejilla e incepcion', si(abs(median(l1, 'omitnan') - 4.1) < 0.06 && abs(pctl(l1, 90) - 5.3) < 0.06), ''};
MAP(end+1, :) = {'Laboratorio: declaracion asegurada, mediana (p90)', '1614', '14.9 ms (16.1)', 'lab_latency_r3m1.m (solo consola)', '10', G, 'fin de ventana', '', ...
    sprintf('%.2f ms (%.2f); cadena completa: mediana %.2f ms, max %.2f ms', median(lM, 'omitnan'), pctl(lM, 90), median(TL.cadena_decl_ms, 'omitnan'), max(TL.cadena_decl_ms)), 'pure', si(abs(median(lM, 'omitnan') - 14.9) < 0.06 && abs(pctl(lM, 90) - 16.1) < 0.06), ''};
MAP(end+1, :) = {'Laboratorio: registros declarados con la clase correcta', '1614', '12/12', 'lab_latency_r3m1.m (solo consola)', '10', G, 'fin de ventana', '', ...
    sprintf('%d/%d (sin puerta, primera declaracion); %d/%d (cadena completa)', nnz(strcmp(TL.singate_resultado, 'correcta')), height(TL), nnz(strcmp(TL.cadena_resultado, 'correcta')), height(TL)), 'pure', si(nnz(strcmp(TL.singate_resultado, 'correcta')) == 12), 'El script original solo busca la racha de la clase esperada; la cadena mira la primera declaracion de cualquier clase'};
MAP(end+1, :) = {'Laboratorio: falsas declaraciones en ventanas sanas de prefalta', '1614', '0 de 8691', 'lab_prefault_security.m (solo consola)', '10', G, 'no aplica', '', ...
    sprintf('%d ventanas aseguradas de %d; cadena completa: %d declaraciones en prefalta', totf, totw, totfa), 'pure', si(totf == 0 && totw == 8691), ''};
end

%% ---- escribir el mapa
TM = cell2table(MAP, 'VariableNames', {'cifra','linea_main_final','valor_publicado','script_origen','M_origen','clasificador_origen','convencion_temporal_origen','csv_guardado','valor_recalculado','protocolo_recalculo','coincide','nota'});
writetable(TM, fullfile(OUT, 'mapa_cifras.csv'), 'QuoteStrings', true);
fprintf('\n%d cifras en mapa_cifras.csv\n', height(TM));
for i = 1:height(TM), fprintf('[%s] %s | publicado %s | recalculado %s\n', TM.coincide{i}, TM.cifra{i}, TM.valor_publicado{i}, TM.valor_recalculado{i}); end
fprintf('Huella al terminar: %s\nFIN r2_final_robustez.\n', r2_t20_md5());

%% ---------------- funciones locales
function s = si(c), if c, s = 'si'; else, s = 'no'; end, end
function k = hk(cv, k2), k = sprintf('H_%s_k%03d', cv, round(1000*k2)); end
function k = fk(cv, truth, k2), k = sprintf('F_%s_%s_k%03d', cv, strrep(truth, '-', ''), round(1000*k2)); end

function S = scan_dc(V, t, fs, f, W, stride, dcrem)
% Clasificacion por ventana con resta opcional del DC estimado por minimos
% cuadrados sobre el ciclo anterior (igual que cls_dcrem_hist de la suite v2).
    n = size(V, 1); starts = (1:stride:(n - W + 1))'; m = numel(starts);
    FO = r2_final_opts(); PF = struct('kirchhoff_thr', FO.kirchhoff_thr, 'v0_mode', FO.v0_mode);
    Hn = 80; tH = (0:Hn-1)'/fs; pinvX = pinv([cos(2*pi*f*tH), sin(2*pi*f*tH), ones(Hn, 1)]);
    S = struct('fs', fs, 'i_start', starts, 't_start', t(starts), 't_end', t(starts + W - 1), 'label_raw', {cell(m, 1)}, ...
        'roundness', nan(m, 1), 'angle_deg', nan(m, 1), 'semi_major', nan(m, 1), 'semi_minor', nan(m, 1), 'bivector_norm', nan(m, 1), 'is_line', false(m, 1));
    for j = 1:m
        s = starts(j); e = s + W - 1; Vw = V(s:e, :);
        if dcrem && (e - Hn + 1) >= 1
            hist = V(e - Hn + 1:e, :);
            for ph = 1:3, coef = pinvX*hist(:, ph); Vw(:, ph) = Vw(:, ph) - coef(3); end
        end
        try
            [ft, p] = r2_ga_classifier_param(Vw, fs, f, PF);
            S.label_raw{j} = char(ft); S.roundness(j) = p.roundness; S.angle_deg(j) = p.angle_deg;
            S.semi_major(j) = p.semi_major; S.semi_minor(j) = p.semi_minor; S.bivector_norm(j) = p.bivector_norm; S.is_line(j) = p.is_line;
        catch
            S.label_raw{j} = 'Error';
        end
    end
end

function det = fire_first(labels, tc, M)
% Primera racha de M etiquetas de falta identicas (como fire_strict del script original)
    det = NaN;
    for i = M:numel(labels)
        w = labels(i-M+1:i);
        if all(~strcmp(w, 'A-B-C')) && numel(unique(w)) == 1, det = tc(i); return, end
    end
end

function v = pctl(x, p)
% Percentil con la definicion de prctile de MATLAB (posiciones (i-0.5)/n)
    x = sort(x(~isnan(x))); n = numel(x);
    if n == 0, v = NaN; return, end
    pos = p/100*n + 0.5;
    if pos <= 1, v = x(1); elseif pos >= n, v = x(n); else, lo = floor(pos); v = x(lo) + (pos - lo)*(x(lo+1) - x(lo)); end
end

function out = rdc(V, t, varargin)
% Cadena de referencia con las opciones finales (r2_final_opts) delante.
    persistent FO
    if isempty(FO), FO = r2_final_opts(); end
    out = r2_decision_chain(V, t, FO, varargin{:});
end
