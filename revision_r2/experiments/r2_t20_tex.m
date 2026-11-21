function r2_t20_tex(opts)
% R2_T20_TEX  Genera, a partir de los CSV de resultados\T20, la tabla compacta para
% el paper (textos\T20_tabla.tex) y la lista de cifras que usan los parrafos del
% manuscrito (resultados\T20\T20_cifras.csv: clave, valor, archivo de origen).
% No calcula nada nuevo: solo agrega sintef_por_registro.csv, lab_por_registro.csv,
% lab_severidad.csv y lab_sin_falta.csv. Se llama desde r2_t20_run_all.
if nargin < 1 || isempty(opts), opts = r2_t20_defaults(); end
d = opts.outdir;
S = readtable(fullfile(d, 'sintef_por_registro.csv'), 'TextType', 'string', 'Delimiter', ',');
L = readtable(fullfile(d, 'lab_por_registro.csv'), 'TextType', 'string', 'Delimiter', ',');
V = readtable(fullfile(d, 'lab_severidad.csv'), 'TextType', 'string', 'Delimiter', ',');
N = readtable(fullfile(d, 'lab_sin_falta.csv'), 'TextType', 'string', 'Delimiter', ',');
S = nomiss(S); L = nomiss(L); V = nomiss(V); N = nomiss(N);
sig1 = string(opts.lab.signals{1});
L1 = L(L.senal == sig1 & L.puerta == "adaptive", :);
f1 = 'sintef_por_registro.csv'; f2 = 'lab_por_registro.csv'; f3 = 'lab_severidad.csv';

%% ---------------------------------------------------------------- filas
R = {};   % {grupo, etiqueta, mascara/tabla, firma}
m = @(c, g) S.cat == c & S.ground == g;
rng2 = @(x) sprintf('%.2f--%.2f', min(x), max(x));
rngd = @(x) rngdeg(x);
R(end+1,:) = {'s', 'SLG, isolated', S(m("SLG", "Pet0"), :), ['$V_0$ ' rng2(S.V0_rel(m("SLG", "Pet0")))]};
R(end+1,:) = {'s', 'SLG, Petersen', S(m("SLG", "Pet1"), :), ['$V_0$ ' rng2(S.V0_rel(m("SLG", "Pet1")))]};
R(end+1,:) = {'s', 'SLG 1~k$\Omega$, isolated', S(m("HISLG", "Pet0"), :), ['$V_0$ ' rng2(S.V0_rel(m("HISLG", "Pet0")))]};
R(end+1,:) = {'s', 'SLG 1~k$\Omega$, Petersen', S(m("HISLG", "Pet1"), :), ['$V_0$ ' rng2(S.V0_rel(m("HISLG", "Pet1")))]};
ll1 = S.cat == "LL" & S.Srated == "S245"; ll2 = S.cat == "LL" & S.Srated == "S490";
R(end+1,:) = {'s', 'L-L, lower rating', S(ll1, :), ['$\Delta\theta$ ' rngd(S.dev_centre_st_deg(ll1))]};
R(end+1,:) = {'s', 'L-L, higher rating', S(ll2, :), ['$\Delta\theta$ ' rngd(S.dev_centre_st_deg(ll2))]};
lp = L1.en_paper == 1; la = L1.conjunto == "Autokontaktor" & L1.tipo ~= "abc"; l3 = L1.tipo == "abc";
R(end+1,:) = {'l', 'L-L, original 12', L1(lp, :), ['$\Delta\theta$ ' rngd(L1.desv_centro_mediana_deg(lp))]};
R(end+1,:) = {'l', 'L-L, further 12', L1(la, :), ['$\Delta\theta$ ' rngd(L1.desv_centro_mediana_deg(la))]};
R(end+1,:) = {'l', 'Three-phase', L1(l3, :), sprintf('$|V_2/V_1|$ %.3f--%.3f', min(L1.V2V1_senal(l3)), max(L1.V2V1_senal(l3)))};

tex = fullfile(opts.textdir, tname(opts, 'tabla'));
fid = fopen(tex, 'w');
fprintf(fid, '%% %s  (generado por experiments/r2_t20_tex.m a partir de resultados/T20/*.csv; no editar a mano)\n', tname(opts, 'tabla'));
fprintf(fid, '%% DESTINO en main_rev2.tex: seccion de resultados, despues de la tabla de precision (subsec:comparison) y antes de la\n');
fprintf(fid, '%% discusion de robustez. Label: tab:external_validation. Una columna. Va con los parrafos de T20_manuscrito.tex.\n');
fprintf(fid, '%% CIFRA: sintef_por_registro.csv (filas EMT) y lab_por_registro.csv (filas de laboratorio, senal %s, puerta adaptive).\n', sig1);
fprintf(fid, '%% Si falta espacio (D8), la tabla puede ir solo en la carta; los parrafos citan las cifras sin depender de ella.\n');
fprintf(fid, '\\begin{table}[t]\n\\centering\n\\caption{\\rtwo{Event-level results of the unchanged decision chain on two data sets not used for its design. $n$: events (records $\\times$ measured buses); OK: secured declaration of the true class; W: secured declaration of another class; ND: no declaration; $t_d$: declaration latency, median (maximum), referred to the end of the window. $V_0$ in per unit of the pre-fault positive-sequence voltage; $\\Delta\\theta$: offset of the ellipse major axis from the sector center, in degrees. The neutral treatment does not affect the line-to-line events, so each L-L row of the EMT set holds 72 distinct cases counted twice.}}\n');
fprintf(fid, '\\label{tab:external_validation}\n\\footnotesize\n\\setlength{\\tabcolsep}{2pt}\n\\begin{tabular}{@{}lrrrrll@{}}\n\\toprule\n');
fprintf(fid, 'Fault set & $n$ & OK & W & ND & $t_d$ (ms) & Signature \\\\\n\\midrule\n');
grp = '';
for i = 1:size(R, 1)
    x = R{i,3};
    if ~strcmp(grp, R{i,1})
        if i > 1, fprintf(fid, '\\midrule\n'); end
        if R{i,1} == 's', fprintf(fid, '\\multicolumn{7}{@{}l}{\\emph{EMT data set (two buses per record)}} \\\\\n');
        else, fprintf(fid, '\\multicolumn{7}{@{}l}{\\emph{Laboratory records (grid-referred voltages)}} \\\\\n'); end
        grp = R{i,1};
    end
    if R{i,1} == 's', oc = x.outcome; dm = x.decl_ms; else, oc = x.resultado; dm = x.decl_ms; end
    ok = oc == "correcta";
    if any(ok), td = sprintf('%.1f (%.1f)', median(dm(ok)), max(dm(ok))); else, td = '--'; end
    fprintf(fid, '%s & %d & %d & %d & %d & %s & %s \\\\\n', R{i,2}, height(x), nnz(ok), nnz(oc == "erronea"), nnz(oc == "sin_declaracion"), td, R{i,4});
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
fclose(fid);

%% ---------------------------------------------------------------- cifras
K = {};
    function put(k, v, src), K(end+1,:) = {k, v, src}; end %#ok<AGROW>
num = @(v, fmt) sprintf(fmt, v);
put('sintef_registros', num(height(S)/2, '%d'), f1); put('sintef_eventos', num(height(S), '%d'), f1);
for c = ["SLG", "HISLG", "LL"]
    for g = ["Pet0", "Pet1"]
        x = S(S.cat == c & S.ground == g, :); ok = x.outcome == "correcta"; p = lower(c + "_" + g);
        put(p + "_n", num(height(x), '%d'), f1); put(p + "_ok", num(nnz(ok), '%d'), f1);
        put(p + "_err", num(nnz(x.outcome == "erronea"), '%d'), f1); put(p + "_nd", num(nnz(x.outcome == "sin_declaracion"), '%d'), f1);
        put(p + "_V0_min", num(min(x.V0_rel), '%.2f'), f1); put(p + "_V0_max", num(max(x.V0_rel), '%.2f'), f1);
        put(p + "_decl_med_ms", num(median(x.decl_ms(ok)), '%.1f'), f1); put(p + "_decl_max_ms", num(max([x.decl_ms(ok); NaN]), '%.1f'), f1);
        put(p + "_decl_cualquier_clase_med_ms", num(median(x.decl_ms, 'omitnan'), '%.1f'), f1); put(p + "_decl_cualquier_clase_max_ms", num(max(x.decl_ms), '%.1f'), f1);
        if c ~= "LL"
            put(p + "_argV0_min", num(min(x.argV0_rel_faulted_deg), '%.0f'), f1); put(p + "_argV0_max", num(max(x.argV0_rel_faulted_deg), '%.0f'), f1);
            put(p + "_ventanas_tierra_pct_min", num(min(x.full_ground_pct), '%.1f'), f1);
            put(p + "_fase_falta_rel_min", num(min(x.vratio_faulted_min), '%.2f'), f1); put(p + "_fase_falta_rel_max", num(max(x.vratio_faulted_min), '%.2f'), f1);
            put(p + "_fase_mas_baja_es_la_de_falta", num(nnz(x.lowest_phase == x.faulted), '%d'), f1);
            llg = (x.truth == "A-G" & x.decl_class == "A-B-G") | (x.truth == "B-G" & x.decl_class == "B-C-G") | (x.truth == "C-G" & x.decl_class == "C-A-G");
            put(p + "_declaradas_LLG_fase_y_siguiente", num(nnz(llg), '%d'), f1);
            put(p + "_declaradas_cualquier_clase_tierra", num(nnz(ismember(x.decl_class, ["A-G","B-G","C-G","A-B-G","B-C-G","C-A-G"])), '%d'), f1);
        end
    end
end
x = S(S.cat == "SLG", :); ok = x.outcome == "correcta";
put('slg_n', num(height(x), '%d'), f1); put('slg_ok', num(nnz(ok), '%d'), f1); put('slg_decl_med_ms', num(median(x.decl_ms(ok)), '%.1f'), f1);
put('slg_decl_p90_ms', num(prctile(x.decl_ms(ok), 90), '%.1f'), f1); put('slg_decl_max_ms', num(max(x.decl_ms(ok)), '%.1f'), f1);
put('slg_retrasadas_n', num(nnz(x.variant == "retrasada"), '%d'), f1); put('slg_retrasadas_ok', num(nnz(x.variant == "retrasada" & ok), '%d'), f1);
for s = ["S245", "S490"]
    x = S(S.cat == "HISLG" & S.ground == "Pet1" & S.Srated == s, :); p = "hislg_pet1_" + lower(s);
    put(p + "_n", num(height(x), '%d'), f1); put(p + "_ok", num(nnz(x.outcome == "correcta"), '%d'), f1);
    put(p + "_argV0_min", num(min(x.argV0_rel_faulted_deg), '%.0f'), f1); put(p + "_argV0_max", num(max(x.argV0_rel_faulted_deg), '%.0f'), f1);
    x = S(S.cat == "LL" & S.Srated == s, :); ok = x.outcome == "correcta"; p = "ll_" + lower(s);
    put(p + "_n", num(height(x), '%d'), f1); put(p + "_ok", num(nnz(ok), '%d'), f1); put(p + "_err", num(nnz(x.outcome == "erronea"), '%d'), f1);
    put(p + "_nd", num(nnz(x.outcome == "sin_declaracion"), '%d'), f1);
    put(p + "_dev_min", num(min(x.dev_centre_st_deg), '%.1f'), f1); put(p + "_dev_med", num(median(x.dev_centre_st_deg), '%.1f'), f1); put(p + "_dev_max", num(max(x.dev_centre_st_deg), '%.1f'), f1);
    put(p + "_decl_med_ms", num(median(x.decl_ms(ok)), '%.1f'), f1); put(p + "_decl_max_ms", num(max([x.decl_ms(ok); NaN]), '%.1f'), f1);
    put(p + "_margen_min_deg", num(min(x.margin_min_full_deg), '%.1f'), f1);
    for b = ["DG2", "Pex"], y = x(x.bus == b, :); put(p + "_" + lower(b) + "_ok", sprintf('%d/%d', nnz(y.outcome == "correcta"), height(y)), f1); end
end
x = S(S.cat == "LL", :);
for g = ["Pet0", "Pet1"]
    for sr = ["S245", "S490"]
        y = x(x.ground == g & x.Srated == sr, :); put("ll_" + lower(sr) + "_" + lower(g) + "_ok", sprintf('%d/%d', nnz(y.outcome == "correcta"), height(y)), f1);
    end
end
y0 = sortrows(x(x.ground == "Pet0", :), {'bus','faulted','Srated','pu','limit'}); y1 = sortrows(x(x.ground == "Pet1", :), {'bus','faulted','Srated','pu','limit'});
put('ll_pet0_pet1_resultados_distintos', num(nnz(y0.outcome ~= y1.outcome | y0.decl_class ~= y1.decl_class), '%d'), f1);
put('ll_pet0_pet1_dif_angulo_max_deg', num(max(abs(y0.dev_centre_st_deg - y1.dev_centre_st_deg)), '%.2f'), f1);
put('ll_redondez_medida', sprintf('%.2f a %.2f', min(x.roundness_st_med), max(x.roundness_st_med)), f1);
put('ll_redondez_fundamental', sprintf('%.2f a %.2f', min(x.roundness_pred), max(x.roundness_pred)), f1);
put('ll_desv_fasores_deg', sprintf('%.1f a %.1f', min(x.dev_centre_pred_deg), max(x.dev_centre_pred_deg)), f1);
put('ll_desv_clasificador_deg', sprintf('%.1f a %.1f', min(x.dev_centre_st_deg), max(x.dev_centre_st_deg)), f1);
put('ll_contenido_no_fundamental_falta_pu', sprintf('%.3f a %.3f', min(x.resid50_rel), max(x.resid50_rel)), f1);
put('ll_contenido_no_fundamental_prefalta_pu', sprintf('%.4f a %.4f', min(x.resid50_pre_rel), max(x.resid50_pre_rel)), f1);
put('ll_dispersion_intraregistro_iqr_max_deg', num(max(x.theta_st_iqr), '%.1f'), f1);
for sr = ["S245", "S490"]
    y = x(x.Srated == sr, :); put("ll_" + lower(sr) + "_desv_fasores_deg", sprintf('%.1f a %.1f', min(y.dev_centre_pred_deg), max(y.dev_centre_pred_deg)), f1);
end
put('ll_V2V1_min', num(min(x.V2V1), '%.2f'), f1); put('ll_V2V1_max', num(max(x.V2V1), '%.2f'), f1);
put('ll_ventanas_tierra_pct_max', num(max(x.full_ground_pct), '%.1f'), f1); put('ll_V0_rel_max', num(max(x.V0_rel), '%.1e'), f1);
put('ll_pred_menos_cls_mediana_abs', num(median(abs(x.pred_minus_cls_deg), 'omitnan'), '%.1f'), f1);
P = readtable(fullfile(d, 'sintef_angulo_LL_pareado.csv'), 'TextType', 'string', 'Delimiter', ',');
for i = 1:height(P)
    put("ll_pareado_" + P.factor(i) + "_" + strrep(P.contraste(i), " ", "_"), sprintf('%+.1f [%+.1f, %+.1f] n=%d', P.dif_mediana_deg(i), P.dif_min_deg(i), P.dif_max_deg(i), P.n_pares(i)), 'sintef_angulo_LL_pareado.csv');
end
put('pre_ventanas', num(sum(S.pre_n_win), '%d'), f1); put('pre_ventanas_armadas', num(sum(S.pre_n_armed), '%d'), f1); put('pre_ventanas_inicializacion', num(sum(S.pre_n_init), '%d'), f1); put('pre_registros_con_declaracion', num(nnz(S.n_decl_pre > 0), '%d'), f1);
put('pre_registros_con_declaracion_sin_puerta', num(nnz(S.n_decl_pre_nogate > 0), '%d'), f1);
w = @(p, n) sum(p.*n)/sum(n);
for b = ["DG2", "Pex"]
    y = S(S.bus == b, :);
    put("pre_" + lower(b) + "_no_sana_final_pct", num(w(y.pre_nonhealthy_pct, y.pre_n_win), '%.1f'), f1);
    put("pre_" + lower(b) + "_no_sana_final_armadas_pct", num(w(y.pre_armed_nonhealthy_pct, y.pre_n_armed), '%.1f'), f1);
    put("pre_" + lower(b) + "_no_sana_final_desde_0p2s_pct", num(w(y.pre_late_nonhealthy_pct, y.pre_late_n_win), '%.1f'), f1); put("pre_" + lower(b) + "_no_sana_cruda_pct", num(w(y.pre_raw_nonhealthy_pct, y.pre_n_win), '%.1f'), f1);
end
put('pre_tierra_cruda_pct', num(w(S.pre_raw_ground_pct, S.pre_n_win), '%.3f'), f1); put('pre_tierra_cruda_desde_0p2s_pct', num(w(S.pre_late_raw_ground_pct, S.pre_late_n_win), '%.3f'), f1);
put('pre_pinmax_max', num(max(S.pinmax_pre_max), '%.1e'), f1); put('pre_V0_rel_max', num(max(S.V0pre_rel), '%.1e'), f1);
put('post_registros_con_declaracion', num(nnz(S.n_decl_after > 0), '%d'), f1);
for sx = ["volts", "pu_nom"]
    if ismember("outcome_" + sx, string(S.Properties.VariableNames))
        put("escala_" + sx + "_resultados_distintos", num(nnz(S.("outcome_" + sx) ~= S.outcome | S.("decl_class_" + sx) ~= S.decl_class), '%d'), f1);
        dm = S.("decl_ms_" + sx); put("escala_" + sx + "_latencias_distintas", num(nnz(~(dm == S.decl_ms | (isnan(dm) & isnan(S.decl_ms)))), '%d'), f1);
        put("escala_" + sx + "_ventanas_distintas", sprintf('%d de %d', sum(S.("ndiff_final_" + sx)), sum(S.n_win)), f1);
        put("escala_" + sx + "_ventanas_distintas_en_falta", num(sum(S.("ndiff_final_fault_" + sx)), '%d'), f1);
    end
end
% laboratorio
grpL = {'lab_paper12', lp; 'lab_auto_LL', la; 'lab_LL24', L1.tipo ~= "abc"; 'lab_abc', l3};
for i = 1:4
    x = L1(grpL{i,2}, :); ok = x.resultado == "correcta"; p = string(grpL{i,1});
    put(p + "_n", num(height(x), '%d'), f2); put(p + "_ok", num(nnz(ok), '%d'), f2); put(p + "_err", num(nnz(x.resultado == "erronea"), '%d'), f2);
    put(p + "_clases", strjoin(unique(x.clase_declarada(x.clase_declarada ~= "" & ~ismissing(x.clase_declarada))), '|'), f2);
    put(p + "_decl_med_ms", num(median(x.decl_ms(ok)), '%.1f'), f2); put(p + "_decl_p90_ms", num(prctile(x.decl_ms(ok), 90), '%.1f'), f2); put(p + "_decl_max_ms", num(max([x.decl_ms(ok); NaN]), '%.1f'), f2);
    put(p + "_primera_med_ms", num(median(x.primera_correcta_ms, 'omitnan'), '%.2f'), f2); put(p + "_primera_p90_ms", num(prctile(x.primera_correcta_ms, 90), '%.2f'), f2);
    put(p + "_primera_min_ms", num(min(x.primera_correcta_ms), '%.2f'), f2); put(p + "_primera_max_ms", num(max(x.primera_correcta_ms), '%.2f'), f2);
    put(p + "_inicio_racha_med_ms", num(median(x.inicio_racha_ms(ok)), '%.2f'), f2); put(p + "_inicio_racha_max_ms", num(max([x.inicio_racha_ms(ok); NaN]), '%.2f'), f2);
    put(p + "_dev_min", num(min(x.desv_centro_mediana_deg), '%.1f'), f2); put(p + "_dev_max", num(max(x.desv_centro_mediana_deg), '%.1f'), f2);
    put(p + "_margen_frontera_min", num(min(x.margen_frontera_min_deg), '%.1f'), f2);
    put(p + "_redondez_med_min", num(min(x.redondez_mediana), '%.3f'), f2); put(p + "_redondez_med_max", num(max(x.redondez_mediana), '%.3f'), f2);
    put(p + "_umbral_min", num(min(x.umbral_redondez), '%.3f'), f2); put(p + "_umbral_max", num(max(x.umbral_redondez), '%.3f'), f2);
    put(p + "_V2V1_min", num(min(x.V2V1_senal), '%.3f'), f2); put(p + "_V2V1_max", num(max(x.V2V1_senal), '%.3f'), f2);
    put(p + "_V2V1_suelo_min", num(min(x.V2V1_suelo), '%.3f'), f2); put(p + "_V2V1_suelo_max", num(max(x.V2V1_suelo), '%.3f'), f2);
    put(p + "_decl_prefalta", num(nnz(x.n_decl_prefalta > 0), '%d'), f2);
    put(p + "_ventanas_prefalta", num(sum(x.pre_lejana_n + x.pre_100ms_n), '%d'), f2);
    put(p + "_pre_lejana_no_sana_pct", sprintf('%.1f a %.1f', min(x.pre_lejana_no_sana_pct), max(x.pre_lejana_no_sana_pct)), f2);
    put(p + "_pre_100ms_no_sana_pct", sprintf('%.1f a %.1f', min(x.pre_100ms_no_sana_pct), max(x.pre_100ms_no_sana_pct)), f2);
end
for sg = string(opts.lab.signals(2:end))
    x = L(L.senal == sg & L.puerta == "adaptive" & L.tipo ~= "abc", :); ok = x.resultado == "correcta";
    put("lab_LL24_" + sg + "_ok", sprintf('%d/%d', nnz(ok), height(x)), f2); put("lab_LL24_" + sg + "_decl_med_ms", num(median(x.decl_ms(ok)), '%.1f'), f2);
    put("lab_LL24_" + sg + "_decl_max_ms", num(max(x.decl_ms(ok)), '%.1f'), f2);
    put("lab_LL24_" + sg + "_dev", sprintf('%+.1f a %+.1f', min(x.desv_centro_mediana_deg), max(x.desv_centro_mediana_deg)), f2);
    put("lab_LL24_" + sg + "_margen_frontera_min", num(min(x.margen_frontera_min_deg), '%.1f'), f2);
    y = L(L.senal == sg & L.puerta == "adaptive" & L.tipo == "abc", :); put("lab_abc_" + sg + "_resultados", strjoin(unique(y.resultado + ":" + y.clase_declarada), ' '), f2);
end
put('lab_decl_prefalta_total', num(nnz(L.n_decl_prefalta > 0), '%d'), f2);
put('lab_ventanas_prefalta_37_registros', num(sum(L1.pre_lejana_n + L1.pre_100ms_n), '%d'), f2);
for cj = ["Outkontaktor", "Autokontaktor"]
    y = L1(L1.conjunto == cj, :); put("lab_" + lower(cj) + "_pre_lejana_no_sana_pct", sprintf('%.1f a %.1f', min(y.pre_lejana_no_sana_pct), max(y.pre_lejana_no_sana_pct)), f2);
    put("lab_" + lower(cj) + "_pre_100ms_no_sana_pct", sprintf('%.1f a %.1f', min(y.pre_100ms_no_sana_pct), max(y.pre_100ms_no_sana_pct)), f2);
end
y = L1(L1.tipo == "abc", :); put('lab_abc_decl_ms', sprintf('%.1f a %.1f', min(y.decl_ms), max(y.decl_ms)), f2);
put('lab_abc_redondez_mediana', sprintf('%.3f a %.3f', min(y.redondez_mediana), max(y.redondez_mediana)), f2);
ll = V.tipo ~= "abc"; a3 = V.tipo == "abc";
put('lab_LL_V2V1_box2', rng2(V.box2_V2V1(ll)), f3); put('lab_LL_lazo_rel_box2', rng2(V.box2_lazo_rel(ll)), f3);
put('lab_LL_arg_residual_box2', sprintf('%.1f a %.1f', min(V.box2_arg_residual_deg(ll)), max(V.box2_arg_residual_deg(ll))), f3);
put('lab_LL_Zf_sobre_Zs', sprintf('%.2f--%.2f, arg %+.0f a %+.0f', min(V.box2_Zf_sobre_Zs_abs(ll)), max(V.box2_Zf_sobre_Zs_abs(ll)), min(V.box2_Zf_sobre_Zs_arg_deg(ll)), max(V.box2_Zf_sobre_Zs_arg_deg(ll))), f3);
put('lab_abc_V1_rel_box2', sprintf('%.3f--%.3f', min(V.box2_V1_rel(a3)), max(V.box2_V1_rel(a3))), f3); put('lab_abc_V2V1_box2', sprintf('%.3f--%.3f', min(V.box2_V2V1(a3)), max(V.box2_V2V1(a3))), f3);
for px = ["box2", "conv1_ref_red"]
    put("lab_abc_" + px + "_V2_prefalta_pu", sprintf('%.4f a %.4f', min(V.(px + "_V2_pre_pu")(a3)), max(V.(px + "_V2_pre_pu")(a3))), f3);
    put("lab_abc_" + px + "_V2_falta_pu", sprintf('%.4f a %.4f', min(V.(px + "_V2_falta_pu")(a3)), max(V.(px + "_V2_falta_pu")(a3))), f3);
    put("lab_abc_" + px + "_V2_falta_menos_prefalta_pu", sprintf('%.4f a %.4f', min(V.(px + "_V2_falta_menos_pre_pu")(a3)), max(V.(px + "_V2_falta_menos_pre_pu")(a3))), f3);
    put("lab_abc_" + px + "_V2V1_prefalta", sprintf('%.4f a %.4f', min(V.(px + "_V2V1_pre")(a3)), max(V.(px + "_V2V1_pre")(a3))), f3);
end
put('lab_abc_Zf_sobre_Zs', sprintf('%.3f--%.3f', min(V.box2_Zf_sobre_Zs_abs(a3)), max(V.box2_Zf_sobre_Zs_abs(a3))), f3);
for tp = ["ab", "bc", "ac"]
    y = V(V.tipo == tp, :); put("lab_" + tp + "_V2V1_dispersion_Out_Auto", sprintf('Out %.4f--%.4f, Auto %.4f--%.4f', min(y.box2_V2V1(y.conjunto == "Outkontaktor")), max(y.box2_V2V1(y.conjunto == "Outkontaktor")), ...
        min(y.box2_V2V1(y.conjunto == "Autokontaktor")), max(y.box2_V2V1(y.conjunto == "Autokontaktor"))), f3);
end
put('lab_P_pre_kW', sprintf('%.1f--%.1f', min(V.P_pre_kW), max(V.P_pre_kW)), f3);
au = V.conjunto == "Autokontaktor";
put('lab_auto_maniobra_registros', sprintf('%d de %d', nnz(~isnan(V.maniobra_previa_ms(au))), nnz(au)), f3);
put('lab_auto_maniobra_ms', sprintf('%.0f a %.0f', min(V.maniobra_previa_ms(au)), max(V.maniobra_previa_ms(au))), f3);
put('lab_auto_escalon_pct', sprintf('%.1f a %.1f', min(V.box2_paso_previo_pct(au)), max(V.box2_paso_previo_pct(au))), f3);
put('lab_out_maniobra_registros', sprintf('%d de %d', nnz(~isnan(V.maniobra_previa_ms(~au))), nnz(~au)), f3);
put('lab_incepcion_respaldo', num(nnz(V.incepcion_respaldo), '%d'), f3);
n1 = N(N.senal == sig1, :);
put('lab_sin_falta_ventanas_armadas', num(sum(n1.n_win - n1.n_init), '%d'), 'lab_sin_falta.csv');
put('lab_sin_falta_registros', num(height(n1), '%d'), 'lab_sin_falta.csv'); put('lab_sin_falta_ventanas', num(sum(n1.n_win), '%d'), 'lab_sin_falta.csv');
put('lab_sin_falta_declaraciones', num(sum(N.n_decl), '%d'), 'lab_sin_falta.csv'); put('lab_sin_falta_no_sana_pct', sprintf('%.1f a %.1f', min(n1.nonhealthy_pct), max(n1.nonhealthy_pct)), 'lab_sin_falta.csv');
KT = cell2table(K, 'VariableNames', {'clave', 'valor', 'origen'});
writetable(KT, fullfile(d, 'T20_cifras.csv'));
fprintf('Tabla: %s\nCifras: %s (%d)\n', tex, fullfile(d, 'T20_cifras.csv'), height(KT));
end

function s = rngdeg(x)
% intervalo de grados; un solo valor si el recorrido es menor de 1 grado
if max(x) - min(x) < 1, s = sprintf('$%+.1f$', median(x)); else, s = sprintf('$%+.0f$ to $%+.0f$', min(x), max(x)); end
end

function T = nomiss(T)
% cadenas vacias en lugar de <missing> (las comparaciones con missing siempre son falsas)
for c = string(T.Properties.VariableNames)
    if isstring(T.(c)), x = T.(c); x(ismissing(x)) = ""; T.(c) = x; end
end
end

function n = tname(opts, what)
if isempty(opts.tag), n = ['T20_' what '.tex']; else, n = ['T20_' opts.tag '_' what '.tex']; end
end
