function R = r2_final_lab(varargin)
% R2_FINAL_LAB  (T90, bloque de laboratorio) Regenera los resultados finales del
% laboratorio con el clasificador final de la ronda 2 (r2_final_opts: el publicado con
% la puerta de secuencia cero G5 de umbral 0.05) y diagnostica las 13 faltas trifasicas.
%
%   r2_final_lab()                    % regenera todo en ronda2\resultados\T90\lab\
%   r2_final_lab('regen', false)      % no repite r2_lab_autokontaktor (usa los CSV ya
%                                     % presentes en la carpeta de salida)
%   r2_final_lab('outdir', carpeta)   % otra carpeta de salida
%
% Senal adoptada: conv1 referida al lado de red con el grupo Dy11 ('conv1_ref_red' de
% r2_lab_signals), en pu, 10 kHz, W = 50, paso 12, M = 10, marcas al final de la ventana.
% Toda decision sale de r2_decision_chain(V, t, r2_final_opts(), ...). No modifica nada.
%
% PARTE 1. Resultados por evento
%   a) r2_lab_autokontaktor con las opciones finales, escribiendo en la carpeta de salida
%      (lab_por_registro.csv, lab_resumen.csv, lab_severidad.csv, lab_sin_falta.csv,
%      lab_severidad_resumen.csv, lab_canales_corriente.csv, lab_log.txt).
%   b) Comparacion celda a celda con ronda2\resultados\T20\G5_005\ (reuso_G5_005.csv) y de
%      las opciones y la huella MD5 de su log con las actuales.
%   c) Segunda pasada propia por los 37 registros con falta y los 2 sin falta, solo con la
%      senal adoptada: lab_final_por_registro.csv y lab_final_resumen.csv (12 L-L del
%      paper, 12 L-L de Autokontaktor, 24 L-L, 5 + 8 trifasicas, sin falta).
% PARTE 2. Diagnostico de las trifasicas
%   trifasicas_diagnostico.csv (una fila por registro), trifasicas_ventanas.csv (serie por
%   ventana de los 13), trifasicas_ventanas_representativo.csv (un registro),
%   trifasicas_variantes.csv (barrido de Delta y otras puertas, solo diagnostico) y
%   angulo_prediccion_LL.csv (comprobacion de la formula del angulo con las 24 L-L).
% CIFRAS. cifras_lab.csv, una fila por cifra que puede ir al paper o a la carta, con su
%   CSV de origen y la linea de paper\main_rev2.tex donde esta el valor antiguo.
% Log: r2_final_lab_log.txt (opciones resueltas y huellas MD5).
%
% Ejecucion (R2024b, desde ronda2\resultados\T90\lab):
%   "C:\Program Files\MATLAB\R2024b\bin\matlab.exe" -batch "addpath('<raiz del proyecto>\experiments'); r2_final_lab"
% No hay azar: no se usan semillas.
p = struct('regen', true, 'outdir', '');
for k = 1:2:numel(varargin), p.(varargin{k}) = varargin{k+1}; end
here = fileparts(mfilename('fullpath')); addpath(here);
RAIZ = fileparts(here); R2 = fullfile(RAIZ, '_revision_TPWRD', 'ronda2');
outdir = p.outdir; if isempty(outdir), outdir = fullfile(R2, 'resultados', 'T90', 'lab'); end
if ~exist(outdir, 'dir'), mkdir(outdir); end
[o, info] = r2_final_opts();
logf = fopen(fullfile(outdir, 'r2_final_lab_log.txt'), 'w'); P = @(varargin) plog(logf, varargin{:});
P('r2_final_lab %s MATLAB %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), version);
P('Opciones finales (r2_final_opts): %s\n', jsonencode(o));
rs = info.resolved; if isa(rs.classifier, 'function_handle'), rs.classifier = func2str(rs.classifier); end
P('Opciones resueltas de la cadena (fs, W y paso se resuelven al llamar: 10 kHz, 50 y 12): %s\n', jsonencode(rs));
md5_ini = r2_t20_md5(); P('Huella de la cadena (r2_t20_md5): %s\n', md5_ini);
extra = {'r2_v0gate_stats.m', 'r2_lab_signals.m', 'r2_lab_inception.m', 'r2_lab_autokontaktor.m', 'r2_t20_defaults.m', 'r2_final_opts.m', 'r2_final_lab.m'};
md5x = cellfun(@(f) sprintf('%s=%s', f, md5file(fullfile(here, f))), extra, 'UniformOutput', false);
P('Huella de los demas archivos que intervienen: %s\n', strjoin(md5x, ' '));

fs = 1e4; f = 50; cyc = fs/f; tol = 0.25/fs; SIG = 'conv1_ref_red';
a = exp(2j*pi/3); A = [1 1 1; 1 a a^2; 1 a^2 a]/3;

%% ================================================================ PARTE 1a
if p.regen
    P('\nPARTE 1a. r2_lab_autokontaktor con las opciones finales, salida en %s\n', outdir);
    opts = r2_t20_defaults('chain', o); opts.outdir = outdir;
    r2_lab_autokontaktor(opts);
else
    P('\nPARTE 1a. No se repite r2_lab_autokontaktor (regen = false); se usan los CSV de %s\n', outdir);
end

%% ================================================================ PARTE 1b
ref = fullfile(R2, 'resultados', 'T20', 'G5_005');
files = {'lab_por_registro.csv', 'lab_resumen.csv', 'lab_severidad.csv', 'lab_severidad_resumen.csv', 'lab_sin_falta.csv', 'lab_canales_corriente.csv'};
P('\nPARTE 1b. Comparacion con %s\n', ref);
reuso_ok = true; rr = {};
for k = 1:numel(files)
    c = cmpcsv(fullfile(outdir, files{k}), fullfile(ref, files{k})); c.archivo = files{k}; rr{end+1} = c; %#ok<AGROW>
    P('  %-28s filas %d, columnas %d, celdas %d, celdas distintas %d, diferencia numerica maxima %.3g, bytes identicos %d\n', files{k}, c.filas, c.columnas, c.celdas, c.celdas_distintas, c.dif_num_max, c.bytes_identicos);
    reuso_ok = reuso_ok && c.celdas_distintas == 0;
end
RR = struct2table([rr{:}]'); RR = movevars(RR, 'archivo', 'Before', 1); RR.archivo = string(RR.archivo);
opt_ok = false; md5_ok = false; md5_ref = ''; opt_ref = '';
lg = fullfile(ref, 'lab_log.txt');
if isfile(lg)
    L = regexp(fileread(lg), '\r?\n', 'split');
    i1 = find(startsWith(L, 'Opciones de la cadena distintas de las de referencia:'), 1);
    i2 = find(startsWith(L, 'Version de la cadena (MD5):'), 1);
    if ~isempty(i1)
        opt_ref = strtrim(extractAfter(L{i1}, 'referencia:')); q = jsondecode(opt_ref);
        opt_ok = isequal(sort(fieldnames(q)), sort(fieldnames(o))) && all(cellfun(@(fn) isequal(q.(fn), o.(fn)), fieldnames(o)));
    end
    if ~isempty(i2), md5_ref = strtrim(extractAfter(L{i2}, '(MD5):')); md5_ok = strcmp(md5_ref, md5_ini); end
end
P('  Opciones del log de G5_005: %s -> iguales a r2_final_opts: %d\n', opt_ref, opt_ok);
P('  Huella MD5 del log de G5_005 igual a la actual: %d\n', md5_ok);
reutilizable = reuso_ok && opt_ok && md5_ok;
P('  => G5_005 corresponde a la configuracion final y se puede reutilizar: %d\n', reutilizable);
RR.opciones_log_iguales(:) = opt_ok; RR.md5_log_igual(:) = md5_ok; RR.reutilizable(:) = reutilizable;
writetable(RR, fullfile(outdir, 'reuso_G5_005.csv'));

%% ================================================================ PARTE 1c
T0 = readtable(fullfile(outdir, 'lab_por_registro.csv'), 'TextType', 'string');
T0 = T0(T0.senal == SIG & T0.puerta == "adaptive", :);
S0 = readtable(fullfile(outdir, 'lab_severidad.csv'), 'TextType', 'string');
truthmap = containers.Map({'ab','bc','ac','abc'}, {'A-B','B-C','C-A','A-B-C-G'});
centre = containers.Map({'A-B','B-C','C-A'}, {45, 165, 105});
P('\nPARTE 1c. Segunda pasada con la senal adoptada (%s), %d registros con falta\n', SIG, height(T0));
rows = {}; diag = {}; wins = {}; vars = {}; angs = {}; nof = {}; resolved_txt = '';
for i = 1:height(T0)
    name = char(T0.registro(i)); set = char(T0.conjunto(i)); tipo = char(T0.tipo(i)); truth = truthmap(tipo);
    r = r2_lab_signals(fullfile(RAIZ, 'data', 'lab', ['Fault_emulator_' set], [name '.mat'])); t = r.t; V = r.sig.(SIG);
    t_on = T0.t_on_s(i); k_on = round(t_on*fs) + 1; assert(abs(t(k_on) - t_on) < 1e-7, 't_on no coincide con una muestra en %s', name);
    t_off = t_on + 0.120;
    args = o; args.fs = fs; args.truth = truth; args.t_on = t_on; args.t_off = t_off; args.t_range = [max(0, t_on - 1.5), t_off + 0.1];
    out = r2_decision_chain(V, t, args); E = out.event; W = out.win;
    if isempty(resolved_txt)
        resolved_txt = sprintf('fs=%g W=%d stride=%d M=%d Delta=%g K=%d gate=%s amb_mode=%s v0_mode=%s kirchhoff_thr=%g gnd_mode=%s scale=%g classifier=%s', ...
            out.opts.fs, out.opts.W, out.opts.stride, out.opts.M, out.opts.Delta, out.opts.K, out.opts.gate, out.opts.amb_mode, out.opts.v0_mode, out.opts.kirchhoff_thr, out.opts.gnd_mode, out.opts.scale, out.opts.classifier);
        P('  Opciones resueltas en la llamada: %s\n', resolved_txt);
    end
    % coherencia con el CSV de r2_lab_autokontaktor
    assert(strcmp(E.outcome, T0.resultado(i)) && strcmp(string(E.decl_class), nz(T0.clase_declarada(i))) && eqn(E.decl_ms, T0.decl_ms(i)), 'La segunda pasada no reproduce lab_por_registro.csv en %s', name);
    % puerta de secuencia cero antigua (opciones por defecto de la cadena)
    argd = rmfield(args, {'v0_mode', 'kirchhoff_thr'}); outd = r2_decision_chain(V, t, argd); Ed = outd.event; Wd = outd.win;
    pre = W.t_end <= t_on + tol; prea = pre & W.armed;
    w = struct(); w.registro = name; w.conjunto = set; w.tipo = tipo; w.verdad = truth; w.en_paper = T0.en_paper(i);
    w.t_on_s = t_on; w.incepcion_respaldo = T0.incepcion_respaldo(i); w.P_pre_kW = T0.P_pre_kW(i);
    w.resultado = E.outcome; w.declarada = E.declared; w.clase_declarada = E.decl_class; w.decl_ms = E.decl_ms;
    w.inicio_racha_ms = E.decl_ms - 1e3*(out.opts.M - 1)*out.scan.stride/fs; w.primera_etiqueta_verdad_ms = E.first_label_ms;
    w.n_decl_falta = E.n_decl_correct + E.n_decl_wrong; w.n_decl_tras_t_off = E.n_decl_after;
    w.n_decl_prefalta = E.n_decl_pre; w.clases_prefalta = E.pre_classes;
    w.ventanas_prefalta = nnz(pre); w.ventanas_prefalta_armadas = nnz(prea);
    w.prefalta_armadas_no_sanas = nnz(prea & ~strcmp(W.label, 'A-B-C')); w.racha_max_prefalta = max([W.run(prea); 0]);
    w.desv_centro_mediana_deg = T0.desv_centro_mediana_deg(i); w.margen_frontera_min_deg = T0.margen_frontera_min_deg(i);
    w.redondez_mediana_falta = T0.redondez_mediana(i); w.umbral_redondez = T0.umbral_redondez(i);
    w.puerta_antigua_resultado = Ed.outcome; w.puerta_antigua_clase = Ed.decl_class; w.puerta_antigua_decl_ms = Ed.decl_ms;
    w.ventanas_etiqueta_distinta_puerta_antigua = nnz(~strcmp(W.label, Wd.label)); w.ventanas_cruda_distinta_puerta_antigua = nnz(~strcmp(W.label_raw, Wd.label_raw));
    w.v0_estadistico_max = max(W.pinmax); w.v0_antiguo_max = max(Wd.pinmax); w.ventanas_rama_tierra = nnz(~W.is_kirchhoff); w.ventanas_rama_tierra_antigua = nnz(~Wd.is_kirchhoff);
    rows{end+1} = w; %#ok<AGROW>
    % ---- fasores de secuencia: prefalta lejana, prefalta cercana, primeros ciclos, regimen
    ifar = (k_on - 12*cyc : k_on - 6*cyc - 1)'; inear = (k_on - 2*cyc - 20 : k_on - 21)'; ifl = (k_on + 2*cyc : k_on + 6*cyc - 1)';
    ic1 = (k_on : k_on + cyc - 1)'; ic2 = (k_on + cyc : k_on + 2*cyc - 1)';
    Sp = seqph(V(ifar,:), t(ifar), f, A); Sn = seqph(V(inear,:), t(inear), f, A); Sf = seqph(V(ifl,:), t(ifl), f, A);
    Sc1 = seqph(V(ic1,:), t(ic1), f, A); Sc2 = seqph(V(ic2,:), t(ic2), f, A);
    reg = W.t_start >= t_on + 0.040 - tol & W.t_end <= t_on + 0.120 + tol;
    if isKey(centre, truth)
        % comprobacion de la formula del angulo con las L-L: theta = 165 - arg(V2/V1)/2 (mod 180)
        g = struct(); g.registro = name; g.conjunto = set; g.verdad = truth; g.V2V1 = abs(Sf(3)/Sf(2)); g.arg_V2V1_deg = angle(Sf(3)/Sf(2))*180/pi;
        g.theta_medido_deg = axmean(W.angle_deg(reg)); g.theta_fasores_deg = mod(165 - g.arg_V2V1_deg/2, 180);
        g.diferencia_deg = wrap90(g.theta_medido_deg - g.theta_fasores_deg); g.redondez_medida = median(W.roundness(reg)); g.redondez_fasores = (1 - g.V2V1)/(1 + g.V2V1);
        angs{end+1} = g; %#ok<AGROW>
        fprintf('  %s %s\n', set, name); continue
    end
    %% ========================================================= PARTE 2 (trifasicas)
    k_cl = clear_est(V, t, k_on, cyc, f, A, abs(Sp(2))); t_cl = t(k_cl);
    prim = W.t_end > t_on + tol & W.t_end <= t_on + 0.040 + tol;                 % primeros dos ciclos (incluye ventanas a caballo)
    rest = W.t_start >= t_on + 0.120 - tol & W.t_end <= t_cl + tol;              % resto de la falta hasta el despeje
    flt = W.t_end > t_on + tol & W.t_end <= t_cl + tol;                          % toda la falta
    full = W.t_start >= t_on - tol & W.t_end <= t_cl + tol;                      % ventanas integras en la falta
    jl = find(pre, 1, 'last'); vnom = sqrt(3)/sqrt(2);
    d = struct(); d.registro = name; d.conjunto = set; d.P_pre_kW = T0.P_pre_kW(i); d.t_on_s = t_on; d.incepcion_respaldo = T0.incepcion_respaldo(i);
    d.despeje_estimado_ms = 1e3*(t_cl - t_on);
    js = find(S0.registro == name & S0.conjunto == set, 1); d.maniobra_previa_ms = S0.maniobra_previa_ms(js);
    % fundamentales
    d.V1_pre_pu = abs(Sp(2)); d.V2_pre_pu = abs(Sp(3)); d.V2V1_pre = abs(Sp(3)/Sp(2)); d.arg_V2V1_pre_deg = angle(Sp(3)/Sp(2))*180/pi;
    d.V2V1_pre_cercana = abs(Sn(3)/Sn(2));
    d.V1_ciclo1_pu = abs(Sc1(2)); d.V2V1_ciclo1 = abs(Sc1(3)/Sc1(2)); d.V1_ciclo2_pu = abs(Sc2(2)); d.V2V1_ciclo2 = abs(Sc2(3)/Sc2(2));
    d.V1_falta_pu = abs(Sf(2)); d.V1_falta_rel = abs(Sf(2))/abs(Sp(2)); d.V2_falta_pu = abs(Sf(3)); d.V2V1_falta = abs(Sf(3)/Sf(2)); d.arg_V2V1_falta_deg = angle(Sf(3)/Sf(2))*180/pi;
    d.V2_falta_menos_pre_pu = abs(Sf(3) - Sp(3)); d.V0_falta_pu = abs(Sf(1));
    % contenido no fundamental
    [hp, rp] = harm(V(ifar,:), t(ifar), f); [hf, rf] = harm(V(ifl,:), t(ifl), f);
    d.no_fund_rms_pre_pu = rp.abs; d.no_fund_rel_pre = rp.rel; d.no_fund_rms_falta_pu = rf.abs; d.no_fund_rel_falta = rf.rel;
    d.h5_pre_pu = hp(2); d.h7_pre_pu = hp(3); d.h5_falta_pu = hf(2); d.h7_falta_pu = hf(3);
    d.h5_rel_pre = hp(2)/hp(1); d.h7_rel_pre = hp(3)/hp(1); d.h5_rel_falta = hf(2)/hf(1); d.h7_rel_falta = hf(3)/hf(1);
    % redondez por tramo
    d.redondez_pre_mediana = median(W.roundness(prea)); d.redondez_pre_p05 = prctile(W.roundness(prea), 5); d.redondez_pre_min = min(W.roundness(prea));
    d.redondez_pre_fasores = (1 - d.V2V1_pre)/(1 + d.V2V1_pre);
    d.redondez_primeros_mediana = median(W.roundness(prim)); d.redondez_primeros_min = min(W.roundness(prim)); d.redondez_primeros_max = max(W.roundness(prim));
    d.redondez_regimen_mediana = median(W.roundness(reg)); d.redondez_regimen_min = min(W.roundness(reg)); d.redondez_regimen_p95 = prctile(W.roundness(reg), 95); d.redondez_regimen_max = max(W.roundness(reg));
    d.redondez_regimen_fasores = (1 - d.V2V1_falta)/(1 + d.V2V1_falta);
    % nivel ambiente y umbral
    d.r_amb_en_t_on = W.r_amb(jl); d.umbral_en_t_on = W.thr(jl);
    d.umbral_regimen_min = min(W.thr(reg)); d.umbral_regimen_mediana = median(W.thr(reg)); d.umbral_regimen_max = max(W.thr(reg));
    d.ventanas_falta_que_alimentan_ambiente = nnz(flt & W.fed);
    d.margen_umbral_menos_redondez_regimen = d.umbral_regimen_mediana - d.redondez_regimen_mediana;
    d.delta_necesaria_mediana = d.r_amb_en_t_on - d.redondez_regimen_mediana; d.delta_necesaria_min = d.r_amb_en_t_on - d.redondez_regimen_min;
    d.V2V1_suelo_puerta = (1 - d.umbral_en_t_on)/(1 + d.umbral_en_t_on);
    % magnitud y angulo
    d.semieje_mayor_rel_regimen = median(W.semi_major(reg))/vnom; d.semieje_mayor_rel_regimen_max = max(W.semi_major(reg))/vnom;
    th = W.angle_deg(reg); d.theta_regimen_deg = axmean(th); dth = wrap90(th - d.theta_regimen_deg); d.theta_regimen_min_deg = d.theta_regimen_deg + min(dth); d.theta_regimen_max_deg = d.theta_regimen_deg + max(dth);
    d.theta_fasores_deg = mod(165 - d.arg_V2V1_falta_deg/2, 180); d.sector_regimen = sector(d.theta_regimen_deg);
    d.margen_frontera_regimen_deg = min(abs(wrap90(d.theta_regimen_deg - [15 75 135])));
    d.theta_primeros_deg = axmean(W.angle_deg(prim)); d.sector_primeros = sector(d.theta_primeros_deg);
    % puerta de secuencia cero
    d.v0_estadistico_max = max(W.pinmax(pre | flt)); d.v0_umbral = out.opts.kirchhoff_thr; d.ventanas_rama_tierra = nnz(~W.is_kirchhoff);
    d.v0_antiguo_max = max(Wd.pinmax(pre | flt)); d.v0_antiguo_umbral = outd.opts.kirchhoff_thr; d.ventanas_rama_tierra_antigua = nnz(~Wd.is_kirchhoff);
    d.ventanas_etiqueta_distinta_puerta_antigua = nnz(~strcmp(W.label, Wd.label)); d.puerta_antigua_clase = Ed.decl_class; d.puerta_antigua_decl_ms = Ed.decl_ms;
    % etiquetas crudas y tras la puerta de redondez
    d.n_primeros = nnz(prim); d.cruda_primeros = label_counts(W.label_raw(prim)); d.final_primeros = label_counts(W.label(prim));
    d.n_regimen = nnz(reg); d.cruda_regimen = label_counts(W.label_raw(reg)); d.final_regimen = label_counts(W.label(reg));
    d.n_falta = nnz(flt); d.cruda_falta = label_counts(W.label_raw(flt)); d.final_falta = label_counts(W.label(flt));
    cnt = @(lab, m, s) nnz(strcmp(lab(m), s));
    d.regimen_final_CA = cnt(W.label, reg, 'C-A'); d.regimen_final_ABCG = cnt(W.label, reg, 'A-B-C-G'); d.regimen_final_AB = cnt(W.label, reg, 'A-B'); d.regimen_final_BC = cnt(W.label, reg, 'B-C'); d.regimen_final_ABC = cnt(W.label, reg, 'A-B-C');
    d.regimen_cruda_CA = cnt(W.label_raw, reg, 'C-A'); d.regimen_cruda_ABCG = cnt(W.label_raw, reg, 'A-B-C-G'); d.regimen_cruda_ABC = cnt(W.label_raw, reg, 'A-B-C');
    d.regimen_bajo_umbral = nnz(reg & W.roundness < W.thr); d.regimen_sobre_umbral = nnz(reg & W.roundness >= W.thr);
    % A-B-C-G: cuantas ventanas y por que no se declara
    abcg = flt & strcmp(W.label, 'A-B-C-G'); d.falta_final_ABCG = nnz(abcg); d.falta_final_ABCG_pct = 100*nnz(abcg)/nnz(flt);
    d.falta_cruda_ABCG = nnz(flt & strcmp(W.label_raw, 'A-B-C-G')); d.racha_max_ABCG = max([W.run(abcg); 0]); d.M = out.opts.M;
    d.ABCG_todas_sobre_umbral = all(W.roundness(abcg) >= W.thr(abcg)); d.ABCG_todas_magnitud_baja = all(W.semi_major(abcg)/vnom < 0.90);
    ja = find(abcg, 1); if isempty(ja), d.primera_ABCG_ms = NaN; else, d.primera_ABCG_ms = 1e3*(W.t_end(ja) - t_on); end
    d.integras_final_ABCG_pct = 100*nnz(full & strcmp(W.label, 'A-B-C-G'))/nnz(full);
    % racha que declara
    D = out.decl(out.decl.t_decl > t_on + tol & out.decl.t_decl <= t_cl + tol, :);
    d.n_declaraciones_falta = height(D); d.clases_declaradas_falta = strjoin(unique(D.class)', '|');
    d.clase_declarada = E.decl_class; d.decl_ms = E.decl_ms; d.inicio_racha_ms = w.inicio_racha_ms;
    d.racha_empieza_en = tern(w.inicio_racha_ms <= 40, 'primeros_ciclos', 'regimen');
    if isempty(D)
        d.racha_longitud = NaN; d.racha_fin_ms = NaN; d.racha_redondez_max = NaN; d.racha_umbral = NaN; d.tras_declarar_n = NaN; d.tras_declarar_misma_clase = NaN; d.tras_declarar_etiquetas = '';
    else
        jd = D.j(1); j1 = jd; while j1 < height(W) && W.run(j1+1) == W.run(j1) + 1, j1 = j1 + 1; end
        d.racha_longitud = W.run(j1); d.racha_fin_ms = 1e3*(W.t_end(j1) - t_on);
        j0 = jd - out.opts.M + 1; d.racha_redondez_max = max(W.roundness(j0:j1)); d.racha_umbral = W.thr(jd);
        post = false(height(W), 1); post(jd+1:end) = true; post = post & flt;
        d.tras_declarar_n = nnz(post); d.tras_declarar_misma_clase = nnz(post & strcmp(W.label, E.decl_class)); d.tras_declarar_etiquetas = label_counts(W.label(post));
    end
    d.otras_clases_declaradas = strjoin(setdiff(unique(D.class)', {E.decl_class}), '|');
    d.primera_C_A_final_ms = msfirst(W, flt & strcmp(W.label, 'C-A'), t_on);
    diag{end+1} = d; %#ok<AGROW>
    % ---- serie por ventana
    sel = W.t_end >= t_on - 0.300 & W.t_end <= t_on + 0.250; n = nnz(sel);
    tr = repmat("tras_despeje", height(W), 1); tr(pre) = "prefalta"; tr(prim) = "primeros_ciclos"; tr(flt & ~prim) = "falta_transicion"; tr(reg) = "regimen"; tr(rest) = "resto_falta";
    Wt = table(repmat(string(name), n, 1), repmat(string(set), n, 1), 1e3*(W.t_start(sel) - t_on), 1e3*(W.t_end(sel) - t_on), tr(sel), W.roundness(sel), W.r_amb(sel), W.thr(sel), ...
        W.armed(sel), W.fed(sel), W.angle_deg(sel), W.semi_major(sel)/vnom, W.semi_minor(sel)/vnom, W.pinmax(sel), Wd.pinmax(sel), string(W.label_raw(sel)), string(W.label(sel)), string(Wd.label(sel)), W.run(sel), W.fire(sel), ...
        'VariableNames', {'registro','conjunto','t_inicio_ms','t_fin_ms','tramo','redondez','r_amb','umbral','armada','alimenta_ambiente','angulo_deg','semieje_mayor_rel','semieje_menor_rel', ...
        'v0_estadistico','v0_antiguo','etiqueta_cruda','etiqueta_final','etiqueta_final_puerta_antigua','racha','asegurada'});
    wins{end+1} = Wt; %#ok<AGROW>
    % ---- variantes de la puerta de redondez sobre el mismo escaneo (solo diagnostico)
    vv = {'final (Delta 0.01)', struct(); 'Delta 0.005', struct('Delta', 0.005); 'Delta 0.015', struct('Delta', 0.015); 'Delta 0.02', struct('Delta', 0.02); 'Delta 0.03', struct('Delta', 0.03); ...
          'Delta 0.05', struct('Delta', 0.05); 'umbral fijo 0.99', struct('gate', 'fixed'); 'sin puerta de redondez', struct('gate', 'none'); 'M = 5', struct('M', 5); 'M = 20', struct('M', 20)};
    for q = 1:size(vv, 1)
        ax = args; fn = fieldnames(vv{q,2}); for u = 1:numel(fn), ax.(fn{u}) = vv{q,2}.(fn{u}); end
        ox = r2_decision_chain(out.scan, [], ax); Wx = ox.win; Ex = ox.event;
        x = struct(); x.registro = name; x.conjunto = set; x.variante = vv{q,1}; x.resultado = Ex.outcome; x.clase_declarada = Ex.decl_class; x.decl_ms = Ex.decl_ms;
        x.umbral_en_t_on = Wx.thr(jl); x.regimen_final_CA = nnz(reg & strcmp(Wx.label, 'C-A')); x.regimen_final_ABCG = nnz(reg & strcmp(Wx.label, 'A-B-C-G')); x.n_regimen = nnz(reg);
        x.racha_max_ABCG = max([Wx.run(flt & strcmp(Wx.label, 'A-B-C-G')); 0]); x.n_decl_prefalta = Ex.n_decl_pre;
        vars{end+1} = x; %#ok<AGROW>
    end
    fprintf('  %s %s (trifasica)\n', set, name);
end
% ---- registros sin falta
for set = {'Outkontaktor', 'Autokontaktor'}
    fl = dir(fullfile(RAIZ, 'data', 'lab', ['Fault_emulator_' set{1}], '*.mat'));
    for k = 1:numel(fl)
        name = fl(k).name(1:end-4); if ~isempty(regexp(name, '^[a-z]+', 'once')), continue, end
        r = r2_lab_signals(fullfile(fl(k).folder, fl(k).name)); ax = o; ax.fs = fs; out = r2_decision_chain(r.sig.(SIG), r.t, ax); W = out.win;
        q = struct(); q.registro = name; q.conjunto = set{1}; q.duracion_s = r.t(end); q.ventanas = height(W); q.ventanas_armadas = nnz(W.armed);
        q.n_declaraciones = height(out.decl); q.armadas_no_sanas = nnz(W.armed & ~strcmp(W.label, 'A-B-C')); q.armadas_no_sanas_pct = 100*q.armadas_no_sanas/q.ventanas_armadas;
        q.racha_max = max(W.run); q.v0_estadistico_max = max(W.pinmax);
        nof{end+1} = q; %#ok<AGROW>
    end
end
T = tostr(struct2table([rows{:}]')); Dg = tostr(struct2table([diag{:}]')); Wn = vertcat(wins{:}); Vr = tostr(struct2table([vars{:}]')); Ag = tostr(struct2table([angs{:}]')); N = tostr(struct2table([nof{:}]'));
writetable(T, fullfile(outdir, 'lab_final_por_registro.csv')); writetable(Dg, fullfile(outdir, 'trifasicas_diagnostico.csv'));
writetable(Wn, fullfile(outdir, 'trifasicas_ventanas.csv')); writetable(Vr, fullfile(outdir, 'trifasicas_variantes.csv'));
writetable(Ag, fullfile(outdir, 'angulo_prediccion_LL.csv')); writetable(N, fullfile(outdir, 'lab_final_sin_falta.csv'));
% registro representativo: el de la figura aparcada del paper (Falla_A-B-C.csv) si esta; si no, el de latencia mediana
rep = "abc13-feb-2025 145115"; if ~any(Dg.registro == rep), [~, im] = min(abs(Dg.decl_ms - median(Dg.decl_ms))); rep = Dg.registro(im); end
writetable(Wn(Wn.registro == rep, :), fullfile(outdir, 'trifasicas_ventanas_representativo.csv'));

%% ---------------------------------------------------------------- resumen por grupo
isLL = T.tipo ~= "abc"; out_ = T.conjunto == "Outkontaktor";
G = {'paper12_LL', isLL & out_; 'Autokontaktor_LL', isLL & ~out_; 'LL_24', isLL; 'Outkontaktor_abc', ~isLL & out_; 'Autokontaktor_abc', ~isLL & ~out_; 'abc_13', ~isLL; 'con_falta_37', true(height(T), 1)};
fid = fopen(fullfile(outdir, 'lab_final_resumen.csv'), 'w');
fprintf(fid, ['grupo,senal,n,declarados,correctos,erroneos,sin_declaracion,clases_declaradas,decl_ms_mediana,decl_ms_p90,decl_ms_min,decl_ms_max,' ...
    'inicio_racha_ms_mediana,inicio_racha_ms_p90,inicio_racha_ms_min,inicio_racha_ms_max,primera_etiqueta_verdad_ms_mediana,primera_etiqueta_verdad_ms_p90,primera_etiqueta_verdad_ms_min,primera_etiqueta_verdad_ms_max,' ...
    'registros_con_decl_prefalta,n_decl_prefalta,ventanas_prefalta,ventanas_prefalta_armadas,prefalta_armadas_no_sanas_pct,racha_max_prefalta,' ...
    'desv_centro_min_deg,desv_centro_max_deg,margen_frontera_min_deg,registros_resultado_distinto_puerta_antigua,ventanas_etiqueta_distinta_puerta_antigua,v0_estadistico_max\n']);
P('\nRESULTADO POR GRUPO, senal %s, cadena final:\n', SIG);
for ig = 1:size(G, 1)
    x = T(G{ig,2}, :); dc = x.declarada; sd = gst(x.decl_ms(dc)); si = gst(x.inicio_racha_ms(dc)); sp = gst(x.primera_etiqueta_verdad_ms);
    difp = nnz(~(x.resultado == x.puerta_antigua_resultado & nz(x.clase_declarada) == nz(x.puerta_antigua_clase) & eqn(x.decl_ms, x.puerta_antigua_decl_ms)));
    fprintf(fid, '%s,%s,%d,%d,%d,%d,%d,%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%d,%d,%d,%d,%.3f,%d,%.2f,%.2f,%.2f,%d,%d,%.3g\n', G{ig,1}, SIG, height(x), nnz(dc), ...
        nnz(x.resultado == "correcta"), nnz(x.resultado == "erronea"), nnz(x.resultado == "sin_declaracion"), label_counts(cellstr(x.verdad + ">" + nz(x.clase_declarada))), sd, si, sp, ...
        nnz(x.n_decl_prefalta > 0), sum(x.n_decl_prefalta), sum(x.ventanas_prefalta), sum(x.ventanas_prefalta_armadas), 100*sum(x.prefalta_armadas_no_sanas)/sum(x.ventanas_prefalta_armadas), max(x.racha_max_prefalta), ...
        min(x.desv_centro_mediana_deg), max(x.desv_centro_mediana_deg), min(x.margen_frontera_min_deg), difp, sum(x.ventanas_etiqueta_distinta_puerta_antigua), max(x.v0_estadistico_max));
    P('  %-18s n %2d | declarados %2d, correctos %2d, erroneos %2d, sin declarar %2d (%s) | declaracion ms med %.2f p90 %.2f max %.2f | inicio de racha med %.2f p90 %.2f max %.2f | prefalta: %d registros y %d declaraciones en %d ventanas armadas (racha maxima %d)\n', ...
        G{ig,1}, height(x), nnz(dc), nnz(x.resultado == "correcta"), nnz(x.resultado == "erronea"), nnz(x.resultado == "sin_declaracion"), label_counts(cellstr(x.verdad + ">" + nz(x.clase_declarada))), ...
        sd(1), sd(2), sd(4), si(1), si(2), si(4), nnz(x.n_decl_prefalta > 0), sum(x.n_decl_prefalta), sum(x.ventanas_prefalta_armadas), max(x.racha_max_prefalta));
end
fprintf(fid, 'sin_falta_2,%s,%d,%d,0,0,%d,,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,%d,%d,%d,%d,%.3f,%d,NaN,NaN,NaN,NaN,NaN,%.3g\n', SIG, height(N), nnz(N.n_declaraciones > 0), nnz(N.n_declaraciones == 0), ...
    nnz(N.n_declaraciones > 0), sum(N.n_declaraciones), sum(N.ventanas), sum(N.ventanas_armadas), 100*sum(N.armadas_no_sanas)/sum(N.ventanas_armadas), max(N.racha_max), max(N.v0_estadistico_max));
fclose(fid);
P('  %-18s n %2d | declaraciones %d en %d ventanas (%d armadas); racha maxima %d\n', 'sin_falta_2', height(N), sum(N.n_declaraciones), sum(N.ventanas), sum(N.ventanas_armadas), max(N.racha_max));
P('  Puerta de secuencia cero antigua (opciones por defecto) frente a la final: registros con resultado, clase o latencia distintos %d de %d; ventanas con etiqueta final distinta %d; valor maximo del estadistico final %.3g y del antiguo %.3g\n', ...
    nnz(~(T.resultado == T.puerta_antigua_resultado & nz(T.clase_declarada) == nz(T.puerta_antigua_clase) & eqn(T.decl_ms, T.puerta_antigua_decl_ms))), height(T), sum(T.ventanas_etiqueta_distinta_puerta_antigua), max(T.v0_estadistico_max), max(T.v0_antiguo_max));
P('  Formula del angulo (theta = 165 - arg(V2/V1)/2) en las %d L-L: diferencia con el angulo medido de %.2f a %.2f grados\n', height(Ag), min(Ag.diferencia_deg), max(Ag.diferencia_deg));

%% ---------------------------------------------------------------- log de las trifasicas
mm = @(v) [min(v), max(v)];
P('\nPARTE 2. DIAGNOSTICO DE LAS %d TRIFASICAS (min a max entre registros):\n', height(Dg));
P('  clase declarada: %s; declaracion %.1f a %.1f ms (mediana %.1f); inicio de la racha %.1f a %.1f ms; la racha empieza en: %s\n', label_counts(cellstr(Dg.clase_declarada)), mm(Dg.decl_ms), median(Dg.decl_ms), mm(Dg.inicio_racha_ms), label_counts(cellstr(Dg.racha_empieza_en)));
P('  V1 en falta %.3f a %.3f pu (%.3f a %.3f de la prefalta); |V2/V1| prefalta %.4f a %.4f, falta %.4f a %.4f; V2 prefalta %.4f a %.4f pu, falta %.4f a %.4f pu; arg(V2/V1) en falta %.1f a %.1f grados\n', ...
    mm(Dg.V1_falta_pu), mm(Dg.V1_falta_rel), mm(Dg.V2V1_pre), mm(Dg.V2V1_falta), mm(Dg.V2_pre_pu), mm(Dg.V2_falta_pu), mm(Dg.arg_V2V1_falta_deg));
P('  redondez: prefalta mediana %.4f a %.4f (fasores %.4f a %.4f); primeros ciclos mediana %.4f a %.4f (minimo %.4f); regimen mediana %.4f a %.4f, maximo %.4f a %.4f (fasores %.4f a %.4f)\n', ...
    mm(Dg.redondez_pre_mediana), mm(Dg.redondez_pre_fasores), mm(Dg.redondez_primeros_mediana), min(Dg.redondez_primeros_min), mm(Dg.redondez_regimen_mediana), mm(Dg.redondez_regimen_max), mm(Dg.redondez_regimen_fasores));
P('  r_amb en t_on %.4f a %.4f; umbral %.4f a %.4f; umbral menos redondez mediana de regimen %.4f a %.4f; ventanas de falta que alimentan el ambiente: %d a %d\n', mm(Dg.r_amb_en_t_on), mm(Dg.umbral_en_t_on), mm(Dg.margen_umbral_menos_redondez_regimen), mm(Dg.ventanas_falta_que_alimentan_ambiente));
P('  contenido no fundamental: prefalta %.4f a %.4f pu rms (%.2f a %.2f %% del fundamental); falta %.4f a %.4f pu rms (%.2f a %.2f %%)\n', mm(Dg.no_fund_rms_pre_pu), 100*mm(Dg.no_fund_rel_pre), mm(Dg.no_fund_rms_falta_pu), 100*mm(Dg.no_fund_rel_falta));
P('  angulo del eje en regimen %.1f a %.1f grados (fasores %.1f a %.1f); sector: %s; margen a la frontera %.1f a %.1f grados; semieje mayor/nominal %.3f a %.3f\n', mm(Dg.theta_regimen_deg), mm(Dg.theta_fasores_deg), label_counts(cellstr(Dg.sector_regimen)), mm(Dg.margen_frontera_regimen_deg), mm(Dg.semieje_mayor_rel_regimen));
P('  regimen (%d a %d ventanas): final C-A %d a %d, final A-B-C-G %d a %d; cruda A-B-C-G %d a %d\n', mm(Dg.n_regimen), mm(Dg.regimen_final_CA), mm(Dg.regimen_final_ABCG), mm(Dg.regimen_cruda_ABCG));
P('  A-B-C-G en toda la falta: %d a %d ventanas (%.1f a %.1f %%); racha maxima %d a %d frente a M = %d; todas con redondez >= umbral: %d; todas con semieje < 0.90: %d\n', ...
    mm(Dg.falta_final_ABCG), mm(Dg.falta_final_ABCG_pct), mm(Dg.racha_max_ABCG), Dg.M(1), all(Dg.ABCG_todas_sobre_umbral), all(Dg.ABCG_todas_magnitud_baja));
P('  racha que declara: longitud %d a %d ventanas, termina a %.1f a %.1f ms; tras declarar, misma clase en %d a %d de %d a %d ventanas; otras clases declaradas despues: "%s"\n', ...
    mm(Dg.racha_longitud), mm(Dg.racha_fin_ms), mm(Dg.tras_declarar_misma_clase), mm(Dg.tras_declarar_n), char(strjoin([unique(Dg.otras_clases_declaradas(Dg.otras_clases_declaradas ~= ""))', ""], ' ')));
P('  puerta de secuencia cero: estadistico final maximo %.3g (umbral %.2f), antiguo %.3g (umbral %.2f); ventanas con etiqueta distinta entre las dos puertas: %d\n', max(Dg.v0_estadistico_max), Dg.v0_umbral(1), max(Dg.v0_antiguo_max), Dg.v0_antiguo_umbral(1), sum(Dg.ventanas_etiqueta_distinta_puerta_antigua));
P('  variantes (solo diagnostico), clase declarada en los %d registros:\n', height(Dg));
uv = unique(Vr.variante, 'stable');
for q = 1:numel(uv)
    x = Vr(Vr.variante == uv(q), :);
    P('     %-24s %s | declaracion mediana %.1f ms | A-B-C-G en regimen %d a %d de %d a %d ventanas | declaraciones en prefalta %d\n', uv(q), label_counts(cellstr(nz(x.clase_declarada))), median(x.decl_ms, 'omitnan'), mm(x.regimen_final_ABCG), mm(x.n_regimen), sum(x.n_decl_prefalta));
end

%% ---------------------------------------------------------------- cifras
C = cifras(T, N, Dg, Vr, T0, S0, RR, fullfile(RAIZ, 'paper', 'main_rev2.tex'));
writetable(C, fullfile(outdir, 'cifras_lab.csv'));
P('\ncifras_lab.csv: %d cifras, %d cambian respecto a main_rev2.tex\n', height(C), nnz(C.cambia == 1));
md5_fin = r2_t20_md5();
if strcmp(md5_ini, md5_fin), P('Huella de la cadena al terminar: la misma que al empezar.\n');
else, P('AVISO: la cadena cambio en disco durante la ejecucion. MD5 al terminar: %s\n', md5_fin); end
fclose(logf);
R = struct('por_registro', T, 'sin_falta', N, 'trifasicas', Dg, 'variantes', Vr, 'cifras', C, 'reuso', RR, 'reutilizable', reutilizable);
end

%% ======================================================================
function C = cifras(T, N, Dg, Vr, T0, S0, RR, texfile)
% Una fila por cifra. csv_origen es relativo a ronda2. La linea se busca en main_rev2.tex.
if isfile(texfile)
    L = regexp(fileread(texfile), '\r?\n', 'split');
else   % sin el manuscrito no se busca la linea de cada cifra
    L = repmat({''}, 1, 1640);
    fprintf('AVISO: no esta %s; cifras_lab.csv se escribe sin la linea del manuscrito.\n', texfile);
end
SC = 'experiments/r2_final_lab.m'; B = 'resultados/T90/lab/';
rows = {};
    function add(id, desc, valor, unidad, csv, antiguo, patron, motivo, cambia)
        ln = '';
        if ~isempty(patron)
            j = find(~cellfun(@isempty, regexp(L(1400:1640), patron, 'once')), 1) + 1399;
            if isempty(j), j = find(~cellfun(@isempty, regexp(L, patron, 'once')), 1); end
            if ~isempty(j), ln = sprintf('%d', j); end
        end
        if nargin < 9, cambia = double(~strcmp(strtrim(valor), strtrim(antiguo))); end
        rows(end+1, :) = {['LAB_' id], desc, valor, unidad, [B csv], SC, cambia, antiguo, ln, motivo}; %#ok<AGROW>
    end
f1 = @(x) sprintf('%.1f', x); f2 = @(x) sprintf('%.2f', x); rg = @(v, fmt) sprintf([fmt ' a ' fmt], min(v), max(v));
isLL = T.tipo ~= "abc"; out_ = T.conjunto == "Outkontaktor";
p12 = T(isLL & out_, :); a12 = T(isLL & ~out_, :); l24 = T(isLL, :); abc = T(~isLL, :);
s12 = S0(S0.conjunto == "Outkontaktor" & S0.tipo ~= "abc", :);
RES = 'lab_final_resumen.csv'; REG = 'lab_final_por_registro.csv';
sen = 'senal conv1 referida al lado de red (antes conv1 sin referir); cadena de T03 con marcas al final de la ventana';
% ---- parrafo de latencia (linea 1616)
add('LL_PAPER_N', 'Registros L-L del paper (Outkontaktor)', sprintf('%d', height(p12)), 'registros', RES, '12', 'twelve line-to-line fault records', 'se confirma');
add('LL_PAPER_REPARTO', 'Reparto por tipo de los 12 L-L del paper (A-B, C-A, B-C)', sprintf('%d A-B, %d C-A, %d B-C', nnz(p12.verdad == "A-B"), nnz(p12.verdad == "C-A"), nnz(p12.verdad == "B-C")), 'registros', REG, '5 A-B, 3 C-A, 4 B-C', 'five A-B, three C-A, four B-C', 'se confirma');
add('SOLAPE', 'Solape entre ventanas con el paso usado (12 muestras de 50)', '76', '%', REG, '75', 'quarter-cycle windows of 50 samples, \$75\\%\$ overlap', 'el paso es floor(50/4) = 12 muestras desde la ronda 1; 75 % es un redondeo');
x = gst(p12.primera_etiqueta_verdad_ms);
add('PRIMERA_MED', sprintf('Primera ventana con la etiqueta correcta, mediana, 12 L-L del paper (exacto %.2f). NO citar: puede ser una etiqueta ambiente de prefalta', x(1)), f1(x(1)), 'ms', RES, '4.1', 'median of \$4\.1\$~ms', sen);
add('PRIMERA_P90', sprintf('Idem, percentil 90 (exacto %.2f)', x(2)), f1(x(2)), 'ms', RES, '5.3', '90th percentile \$5\.3\$~ms', 'coincide por redondeo; el significado cambia como en LAB_PRIMERA_MED');
add('PRIMERA_RANGO', sprintf('Idem, minimo a maximo (exacto %.2f a %.2f)', x(3), x(4)), sprintf('%.1f a %.1f', x(3), x(4)), 'ms', RES, '0.5 a 5.3', 'range \$0\.5\$--\$5\.3\$~ms', sen);
x = gst(p12.inicio_racha_ms);
add('INICIO_RACHA_MED', sprintf('Fin de la primera ventana de la racha que declara, mediana, 12 L-L del paper (exacto %.2f). Sustituye a la primera etiqueta correcta', x(1)), f1(x(1)), 'ms', RES, '', '', 'cifra nueva: no depende de etiquetas ambiente', 1);
add('INICIO_RACHA_P90', sprintf('Idem, percentil 90 (exacto %.2f)', x(2)), f1(x(2)), 'ms', RES, '', '', 'cifra nueva', 1);
add('INICIO_RACHA_MAX', sprintf('Idem, maximo (exacto %.2f)', x(4)), f1(x(4)), 'ms', RES, '', '', 'cifra nueva', 1);
x = gst(p12.decl_ms);
add('DECL_MED', sprintf('Declaracion asegurada (M = 10), mediana, 12 L-L del paper (exacto %.2f)', x(1)), f1(x(1)), 'ms', RES, '14.9', 'median of \$14\.9\$~ms', sen);
add('DECL_P90', sprintf('Idem, percentil 90 (exacto %.2f)', x(2)), f1(x(2)), 'ms', RES, '16.1', '90th percentile \$16\.1\$~ms', sen);
add('DECL_MIN', sprintf('Idem, minimo (exacto %.2f)', x(3)), f1(x(3)), 'ms', RES, '', '', 'cifra nueva', 1);
add('DECL_MAX', sprintf('Idem, maximo (exacto %.2f)', x(4)), f1(x(4)), 'ms', RES, '', '', 'cifra nueva', 1);
add('DECL_BAJO_UN_CICLO', 'La declaracion maxima de los 12 L-L del paper queda por debajo de un ciclo (20 ms)', tern(x(4) < 20, 'si', 'no'), '', RES, 'si', 'below one cycle', 'se confirma');
add('LL_PAPER_CORRECTAS', 'Registros L-L del paper declarados con el par correcto', sprintf('%d de %d', nnz(p12.resultado == "correcta"), height(p12)), 'registros', RES, '12 de 12', 'All twelve records were declared with the correct fault type', 'se confirma');
add('LATENCIA_NEGATIVA', 'Registros L-L del paper con declaracion anterior a la incepcion estimada', sprintf('%d', nnz(p12.decl_ms <= 0)), 'registros', REG, '0', 'no record produced a negative latency', 'se confirma');
add('PREFALTA_FALSAS', 'Declaraciones falsas en la prefalta de los 12 L-L del paper', sprintf('%d', sum(p12.n_decl_prefalta)), 'declaraciones', RES, '0', 'zero false declarations', 'se confirma');
add('PREFALTA_VENTANAS', 'Ventanas de prefalta armadas (pueden declarar) de los 12 L-L del paper: hasta 1.5 s antes de la incepcion, menos las 80 de inicializacion', sprintf('%d', sum(p12.ventanas_prefalta_armadas)), 'ventanas', RES, '8691', '8\\,691', 'convencion: rejilla de la cadena de T03 con 1.5 s de prefalta por registro y solo ventanas armadas');
add('PREFALTA_VENTANAS_TOTAL', 'Idem, contando tambien las ventanas de inicializacion', sprintf('%d', sum(p12.ventanas_prefalta)), 'ventanas', RES, '8691', '8\\,691', 'convencion (ver LAB_PREFALTA_VENTANAS)');
add('LL_PAPER_DESV_CENTRO', 'Desviacion del eje mayor respecto al centro del sector, mediana por registro, 12 L-L del paper', rg(p12.desv_centro_mediana_deg, '%.1f'), 'grados', REG, '', '', 'cifra nueva (senal referida al lado de red)', 1);
add('LL_PAPER_MARGEN_FRONTERA', 'Margen minimo por ventana a la frontera de sector, 12 L-L del paper', f1(min(p12.margen_frontera_min_deg)), 'grados', REG, '', '', 'cifra nueva', 1);
% ---- Autokontaktor y 24 L-L
x = gst(a12.decl_ms);
add('LL_AUTOK_CORRECTAS', 'Registros L-L de Autokontaktor declarados con el par correcto', sprintf('%d de %d', nnz(a12.resultado == "correcta"), height(a12)), 'registros', RES, '', '', 'cifra nueva', 1);
add('LL_AUTOK_DECL', sprintf('Declaracion asegurada, 12 L-L de Autokontaktor: mediana, p90, maximo (exacto %.2f, %.2f, %.2f)', x(1), x(2), x(4)), sprintf('%.1f, %.1f, %.1f', x(1), x(2), x(4)), 'ms', RES, '', '', 'cifra nueva', 1);
x = gst(l24.decl_ms); y = gst(l24.inicio_racha_ms);
add('LL24_CORRECTAS', 'Registros L-L (24) declarados con el par correcto', sprintf('%d de %d', nnz(l24.resultado == "correcta"), height(l24)), 'registros', RES, '', '', 'cifra nueva', 1);
add('LL24_DECL', sprintf('Declaracion asegurada, 24 L-L: mediana, p90, maximo (exacto %.2f, %.2f, %.2f)', x(1), x(2), x(4)), sprintf('%.1f, %.1f, %.1f', x(1), x(2), x(4)), 'ms', RES, '', '', 'cifra nueva', 1);
add('LL24_INICIO_RACHA', sprintf('Inicio de la racha que declara, 24 L-L: mediana, p90, maximo (exacto %.2f, %.2f, %.2f)', y(1), y(2), y(4)), sprintf('%.1f, %.1f, %.1f', y(1), y(2), y(4)), 'ms', RES, '', '', 'cifra nueva', 1);
add('LL24_DESV_CENTRO', 'Desviacion del eje respecto al centro del sector, 24 L-L', rg(l24.desv_centro_mediana_deg, '%.1f'), 'grados', REG, '', '', 'cifra nueva', 1);
% ---- trifasicas (linea 1620)
x = gst(abc.decl_ms); y = gst(abc.inicio_racha_ms);
add('ABC_N', 'Registros trifasicos (5 de Outkontaktor y 8 de Autokontaktor)', sprintf('%d', height(abc)), 'registros', RES, '1', 'In the recorded three-phase fault', 'el texto habla de un unico registro; la campana tiene 13');
add('ABC_DECLARADAS', 'Trifasicas con alguna declaracion durante la falta', sprintf('%d de %d', nnz(abc.declarada), height(abc)), 'registros', RES, '', '', 'cifra nueva', 1);
add('ABC_CLASE', 'Clase que la cadena declara en las trifasicas', char(label_counts(cellstr(nz(abc.clase_declarada)))), 'registros', RES, 'trifasica (por comparacion de semiejes)', 'identified through a comparison of semi-axes', 'la cadena de decision no asegura la clase trifasica: la redondez del hueco cae bajo el umbral adaptativo y se etiqueta por sector', 1);
add('ABC_COMO_TRIFASICA', 'Trifasicas declaradas como A-B-C-G', sprintf('%d de %d', nnz(abc.clase_declarada == "A-B-C-G"), height(abc)), 'registros', RES, '', 'identified through a comparison of semi-axes', 'idem', 1);
add('ABC_DECL', sprintf('Declaracion (como C-A) en las trifasicas: mediana, p90, maximo (exacto %.2f, %.2f, %.2f; minimo %.2f)', x(1), x(2), x(4), x(3)), sprintf('%.1f, %.1f, %.1f', x(1), x(2), x(4)), 'ms', RES, '', '', 'cifra nueva', 1);
add('ABC_INICIO_RACHA', sprintf('Inicio de la racha que declara en las trifasicas: mediana, p90, maximo (exacto %.2f, %.2f, %.2f; minimo %.2f)', y(1), y(2), y(4), y(3)), sprintf('%.1f, %.1f, %.1f', y(1), y(2), y(4)), 'ms', RES, '', '', 'cifra nueva', 1);
DG = 'trifasicas_diagnostico.csv';
add('ABC_V1', 'Tension de secuencia directa durante la falta trifasica', rg(Dg.V1_falta_pu, '%.3f'), 'pu', DG, '', '', 'cifra nueva', 1);
add('ABC_V1_REL', 'Idem, relativa a la prefalta', rg(Dg.V1_falta_rel, '%.3f'), '', DG, '', '', 'cifra nueva', 1);
add('ABC_V2V1_FALTA', '|V2/V1| residual durante la falta trifasica', rg(100*Dg.V2V1_falta, '%.2f'), '%', DG, '', '', 'cifra nueva', 1);
add('ABC_V2V1_PRE', '|V2/V1| en prefalta de los registros trifasicos', rg(100*Dg.V2V1_pre, '%.2f'), '%', DG, '', '', 'cifra nueva', 1);
add('ABC_V2_FALTA', 'Tension de secuencia inversa durante la falta trifasica', rg(100*Dg.V2_falta_pu, '%.2f'), '% de la nominal', DG, '', '', 'cifra nueva', 1);
add('ABC_V2_PRE', 'Tension de secuencia inversa en prefalta', rg(100*Dg.V2_pre_pu, '%.2f'), '% de la nominal', DG, '', '', 'cifra nueva', 1);
add('ABC_REDONDEZ_REGIMEN', 'Redondez en regimen de falta trifasica, mediana por registro', rg(Dg.redondez_regimen_mediana, '%.3f'), '', DG, '', '', 'cifra nueva', 1);
add('ABC_REDONDEZ_FASORES', 'Redondez que predice la componente fundamental, (1-k)/(1+k) con k = |V2/V1|', rg(Dg.redondez_regimen_fasores, '%.3f'), '', DG, '', '', 'cifra nueva', 1);
add('ABC_UMBRAL', 'Umbral adaptativo de redondez aplicado (r_amb - 0.01), fijado en prefalta', rg(Dg.umbral_en_t_on, '%.3f'), '', DG, '', '', 'cifra nueva', 1);
add('ABC_R_AMB', 'Nivel ambiente de redondez r_amb en la incepcion', rg(Dg.r_amb_en_t_on, '%.3f'), '', DG, '', '', 'cifra nueva', 1);
add('ABC_ANGULO', 'Angulo del eje mayor en regimen de falta trifasica (sector C-A: 75 a 135)', rg(Dg.theta_regimen_deg, '%.1f'), 'grados', DG, '0', 'the inclination angle may be asserted to be zero', 'el desequilibrio residual fija un eje medible dentro del sector C-A', 1);
add('ABC_ABCG_PCT', 'Ventanas de la falta trifasica con etiqueta final A-B-C-G', rg(Dg.falta_final_ABCG_pct, '%.1f'), '%', DG, '', '', 'cifra nueva', 1);
add('ABC_ABCG_RACHA_MAX', 'Racha maxima de A-B-C-G en una falta trifasica (hacen falta 10)', sprintf('%d', max(Dg.racha_max_ABCG)), 'ventanas', DG, '', '', 'cifra nueva', 1);
add('ABC_NO_FUND_PRE', 'Contenido no fundamental en prefalta, relativo al fundamental', rg(100*Dg.no_fund_rel_pre, '%.1f'), '%', DG, '', '', 'cifra nueva', 1);
add('ABC_NO_FUND_FALTA', 'Contenido no fundamental en la falta trifasica, relativo al fundamental', rg(100*Dg.no_fund_rel_falta, '%.1f'), '%', DG, '', '', 'cifra nueva', 1);
add('ABC_PUERTA_V0', 'Registros trifasicos cuyo resultado cambia entre la puerta de secuencia cero final (G5, 0.05) y la antigua', sprintf('%d de %d', nnz(~(abc.resultado == abc.puerta_antigua_resultado & nz(abc.clase_declarada) == nz(abc.puerta_antigua_clase) & eqn(abc.decl_ms, abc.puerta_antigua_decl_ms))), height(abc)), 'registros', REG, '', '', 'los registros suman cero por construccion: la puerta de secuencia cero no interviene', 1);
x = Vr(Vr.variante == "Delta 0.02", :);
add('ABC_DELTA_002', 'Diagnostico: clase declarada en las trifasicas si el margen fuera Delta = 0.02 (no es la configuracion final)', char(label_counts(cellstr(nz(x.clase_declarada)))), 'registros', 'trifasicas_variantes.csv', '', '', 'solo diagnostico', 1);
% ---- seguridad
add('PREFALTA_37_FALSAS', 'Declaraciones en prefalta en los 37 registros con falta', sprintf('%d', sum(T.n_decl_prefalta)), 'declaraciones', RES, '', '', 'cifra nueva', 1);
add('PREFALTA_37_VENTANAS', 'Ventanas de prefalta armadas en los 37 registros con falta', sprintf('%d', sum(T.ventanas_prefalta_armadas)), 'ventanas', RES, '', '', 'cifra nueva', 1);
add('SIN_FALTA_DECL', 'Declaraciones en los 2 registros sin falta (10 s cada uno)', sprintf('%d', sum(N.n_declaraciones)), 'declaraciones', 'lab_final_sin_falta.csv', '', '', 'cifra nueva', 1);
add('SIN_FALTA_VENTANAS', 'Ventanas armadas en los 2 registros sin falta', sprintf('%d', sum(N.ventanas_armadas)), 'ventanas', 'lab_final_sin_falta.csv', '', '', 'cifra nueva', 1);
add('PUERTA_V0_37', 'Registros con falta cuyo resultado, clase o latencia cambia entre la puerta de secuencia cero final y la antigua', sprintf('%d de %d', nnz(~(T.resultado == T.puerta_antigua_resultado & nz(T.clase_declarada) == nz(T.puerta_antigua_clase) & eqn(T.decl_ms, T.puerta_antigua_decl_ms))), height(T)), 'registros', REG, '', '', 'la puerta G5 no cambia ningun resultado del laboratorio', 1);
add('REUSO_G5_005', 'resultados/T20/G5_005 coincide celda a celda con la regeneracion final (6 CSV), con las mismas opciones y la misma huella MD5', tern(all(RR.reutilizable), 'si', 'no'), '', 'reuso_G5_005.csv', '', '', 'comprobacion', 0);
% ---- montaje (linea 1446)
SV = 'lab_severidad.csv';
add('I_PRE', 'Corriente de carga previa, rms por fase, 12 L-L del paper', rg(s12.I_pre_rms_A, '%.1f'), 'A', SV, '3.6', 'pre-fault load currents of \$3\.6\$~A', 'el texto solo recogia el nivel de 4 de los 12 registros');
add('P_PRE', 'Potencia previa de los 12 L-L del paper', rg(s12.P_pre_kW, '%.1f'), 'kW', SV, '2.5 (kVA)', '2\.5\$~kVA', 'idem');
add('I_PICO', 'Pico de la corriente de falta en bornes del convertidor, 12 L-L del paper', rg(s12.I_pico_falta_A, '%.1f'), 'A', SV, '35 a 47', '35--47~A', 'idem');
add('VDC', 'Tension continua en prefalta, 12 L-L del paper', sprintf('%.0f', median(s12.Vdc_pre_V)), 'V', SV, '594', 'DC-link voltage of 594~V', 'se confirma');
add('V2V1_LL', '|V2/V1| en regimen de falta, 24 L-L, senal adoptada', rg(T0.V2V1_senal(T0.tipo ~= "abc"), '%.2f'), '', 'lab_por_registro.csv', '', '', 'cifra nueva', 1);
C = cell2table(rows, 'VariableNames', {'id','descripcion','valor','unidad','csv_origen','script','cambia','valor_antiguo','linea_main_rev2','motivo'});
for c = ["id","descripcion","valor","unidad","csv_origen","script","valor_antiguo","linea_main_rev2","motivo"], C.(c) = string(C.(c)); end
end

%% ======================================================================
function c = cmpcsv(fa, fb)
% Comparacion celda a celda de dos CSV (NaN igual a NaN; numeros con igualdad exacta).
c = struct('filas', NaN, 'columnas', NaN, 'celdas', NaN, 'celdas_distintas', NaN, 'dif_num_max', NaN, 'bytes_identicos', false);
if ~isfile(fa) || ~isfile(fb), return, end
c.bytes_identicos = isequal(fileread(fa), fileread(fb));
A = readcell(fa, 'Delimiter', ','); B = readcell(fb, 'Delimiter', ',');
c.filas = size(A, 1) - 1; c.columnas = size(A, 2); c.celdas = c.filas*c.columnas;
if ~isequal(size(A), size(B)), c.celdas_distintas = Inf; return, end
nd = 0; md = 0;
for i = 1:numel(A)
    x = A{i}; y = B{i}; nx = isnumeric(x) && isscalar(x); ny = isnumeric(y) && isscalar(y);
    if nx && ny
        if isnan(x) && isnan(y), continue, end
        dd = abs(x - y); if ~(dd == 0), nd = nd + 1; end
        if isfinite(dd), md = max(md, dd); end
    elseif ~isequal(s1(x), s1(y))
        nd = nd + 1;
    end
end
c.celdas_distintas = nd; c.dif_num_max = md;
end
function s = s1(x), if isa(x, 'missing'), s = ""; else, s = string(x); end, end
function h = md5file(fn)
if ~isfile(fn), h = 'NO_EXISTE'; return, end
fid = fopen(fn, 'r'); b = fread(fid, Inf, '*uint8'); fclose(fid);
m = java.security.MessageDigest.getInstance('MD5'); m.update(b); h = sprintf('%02x', typecast(m.digest(), 'uint8'));
end
function [S, ph] = seqph(v, t, f, A)
% Fasores de 50 Hz por minimos cuadrados (con termino continuo) y componentes de secuencia [V0; V1; V2]
X = [cos(2*pi*f*t), sin(2*pi*f*t), ones(numel(t), 1)];
c = X \ v; ph = (c(1,:) - 1j*c(2,:)).'; S = A*ph;
end
function [h, r] = harm(v, t, f)
% Amplitud (raiz de la media cuadratica entre fases) del fundamental y de los armonicos 5 y 7, y
% contenido no fundamental: lo que queda al restar el ajuste de 50 Hz mas continua.
X = [cos(2*pi*f*t), sin(2*pi*f*t), ones(numel(t), 1)]; c = X \ v; res = v - X*c;
fund = v - res - c(3,:);
r.abs = sqrt(mean(res.^2, 'all')); r.rel = sqrt(sum(res.^2, 'all')/sum(fund.^2, 'all'));
hh = [1 5 7]; Xh = ones(numel(t), 1);
for k = hh, Xh = [Xh, cos(2*pi*k*f*t), sin(2*pi*k*f*t)]; end %#ok<AGROW>
ch = Xh \ v; h = zeros(1, numel(hh));
for k = 1:numel(hh), h(k) = sqrt(mean(ch(2*k,:).^2 + ch(2*k+1,:).^2)); end
end
function k_cl = clear_est(V, t, k_on, cyc, f, A, V1pre)
% Despeje estimado: centro del primer ciclo, a partir de 40 ms, cuya secuencia directa supera
% 0.6 de la de prefalta (resolucion de 1 ms).
k_cl = min(k_on + 8*cyc, size(V, 1));
for k = k_on + 2*cyc : 10 : min(k_on + 14*cyc, size(V, 1) - cyc)
    S = seqph(V(k:k+cyc-1,:), t(k:k+cyc-1), f, A);
    if abs(S(2)) > 0.6*V1pre, k_cl = k + cyc/2; return, end
end
end
function s = sector(th)
% Igual que sector_ de r2_decision_chain y classify_by_angle del clasificador
th = mod(th, 180);
if th >= 15 && th <= 75, s = 'A-B'; elseif th >= 75 && th <= 135, s = 'C-A'; else, s = 'B-C'; end
end
function v = msfirst(W, m, t_on), j = find(m, 1); if isempty(j), v = NaN; else, v = 1e3*(W.t_end(j) - t_on); end, end
function m = axmean(th), th = th(isfinite(th)); m = mod(angle(mean(exp(2j*deg2rad(th))))/2*180/pi, 180); end
function y = wrap90(x), y = mod(x + 90, 180) - 90; end
function s = gst(x), x = x(isfinite(x)); if isempty(x), s = [NaN NaN NaN NaN]; else, s = [median(x), prctile(x, 90), min(x), max(x)]; end, end
function s = nz(s), s = string(s); s(ismissing(s)) = ""; end
function tf = eqn(a, b), tf = (isnan(a) & isnan(b)) | abs(a - b) < 1e-6; end
function v = tern(c, a, b), if c, v = a; else, v = b; end, end
function s = label_counts(lab)
if isempty(lab), s = ''; return, end
lab = cellstr(lab); lab(cellfun(@isempty, lab)) = {'-'};
[u, ~, j] = unique(lab); c = accumarray(j(:), 1); [c, o] = sort(c, 'descend'); u = u(o);
s = strjoin(arrayfun(@(i) sprintf('%s:%d', u{i}, c(i)), 1:numel(u), 'UniformOutput', false), ' ');
end
function T = tostr(T)
for c = string(T.Properties.VariableNames)
    if iscell(T.(c)), T.(c) = string(T.(c)); elseif ischar(T.(c)), T.(c) = string(cellstr(T.(c))); end
end
end
function plog(fid, varargin), fprintf(varargin{:}); fprintf(fid, varargin{:}); end
