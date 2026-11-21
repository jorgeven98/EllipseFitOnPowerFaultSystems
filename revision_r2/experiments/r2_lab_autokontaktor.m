function [T, S] = r2_lab_autokontaktor(opts)
% R2_LAB_AUTOKONTAKTOR  (T20; registros de laboratorio, peticiones R1.2, R4.2a, E2)
% Pasa la cadena de decision de referencia (r2_decision_chain, T03) por los 39
% registros del laboratorio: Fault_emulator_Outkontaktor (12 L-L del paper y 5
% trifasicos) y Fault_emulator_Autokontaktor (12 L-L, 8 trifasicos y 2 sin falta).
%
%   r2_lab_autokontaktor()         % cadena de referencia
%   r2_lab_autokontaktor(opts)     % opts = r2_t20_defaults(...)
%
% Senales (r2_lab_signals): conv1 referida al lado de red (PRINCIPAL), meas_box2 y
% conv1 fase-estrella (lo publicado). 10 kHz, W = 50, paso 12, M = 10.
% Incepcion: r2_lab_inception sobre conv1 fase-estrella (el criterio del paper y de
% T12), contrastada con un segundo estimador en meas_box2; si difieren mas de 1 ms
% se usa el segundo y el registro queda marcado. Un solo t_on por registro para las
% tres senales; t_off = t_on + 120 ms (ventanas dentro de la falta, como T12).
% Por registro: tipo segun el nombre, severidad medida, cociente aparente de
% impedancias, que miden los canales de corriente, resultado por evento, latencias,
% angulo de la elipse y margen respecto al suelo de la puerta de redondez.
% Salidas en resultados\T20\ (o resultados\T20\<tag>\). Solo lectura de datos.
if nargin < 1 || isempty(opts), opts = r2_t20_defaults(); end
here = fileparts(mfilename('fullpath')); addpath(here);
RAIZ = opts.raiz; outdir = opts.outdir;
logf = fopen(fullfile(outdir, 'lab_log.txt'), 'w'); P = @(varargin) plog(logf, varargin{:});
P('r2_lab_autokontaktor %s MATLAB %s tag "%s"\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), version, opts.tag);
P('Opciones de la cadena distintas de las de referencia: %s\n', jsonencode(opts.chain));
md5_ini = r2_t20_md5(); P('Version de la cadena (MD5): %s\n', md5_ini);
fs = 1e4; f = 50; cyc = fs/f; tol = 0.25/fs;
a = exp(2j*pi/3); A = [1 1 1; 1 a a^2; 1 a^2 a]/3;
truthmap = containers.Map({'ab','bc','ac','abc'}, {'A-B','B-C','C-A','A-B-C-G'});
centre = containers.Map({'A-B','B-C','C-A'}, {45, 165, 105});
argnom = containers.Map({'A-B','B-C','C-A'}, {-120, 0, 120});       % arg(V2/V1) de una falta franca
pairix = containers.Map({'A-B','B-C','C-A'}, {[1 2], [2 3], [3 1]});
sets = {'Outkontaktor', 'Autokontaktor'};
signals = opts.lab.signals; gates = {'adaptive', 'none'};
paper12 = @(set, tipo) strcmp(set, 'Outkontaktor') && ~strcmp(tipo, 'abc');

rows = {}; sev = {}; nof = {};
for is = 1:2
    fl = dir(fullfile(RAIZ, 'data', 'lab', ['Fault_emulator_' sets{is}], '*.mat'));
    for k = 1:numel(fl)
        r = r2_lab_signals(fullfile(fl(k).folder, fl(k).name)); t = r.t;
        name = fl(k).name(1:end-4); tipo = regexp(name, '^[a-z]+', 'match', 'once'); if isempty(tipo), tipo = 'sin_falta'; end
        %% ---- registro sin falta: tasas sobre todo el registro
        if strcmp(tipo, 'sin_falta')
            for s = 1:numel(signals)
                args = opts.chain; args.fs = fs; out = r2_decision_chain(r.sig.(signals{s}), t, args);
                q = out.summary; q.registro = name; q.conjunto = sets{is}; q.senal = signals{s}; q.dur_s = t(end);
                S0 = seqph(r.sig.(signals{s})(20001:60000,:), t(20001:60000), f, A); q.V1_pu = abs(S0(2)); q.V2V1 = abs(S0(3)/S0(2));
                q.P_kW = mean(r.Psi); q.thr_med = median(out.win.thr(out.win.armed)); q.roundness_min = min(out.win.roundness);
                nof{end+1} = q; %#ok<AGROW>
            end
            continue
        end
        truth = truthmap(tipo);
        %% ---- incepcion
        [k_ref, ~] = r2_lab_inception(r.sig.conv1_fase, t, fs, f);          % criterio del paper y de T12
        [k_chk, k_cd, ncy] = incep_check(r.sig.box2, t, fs, f, A);            % contraste en meas_box2
        fallback = isnan(k_ref) || abs(k_ref - k_chk) > 10;
        if fallback, k_on = k_chk; else, k_on = k_ref; end
        t_on = t(k_on); t_off = t_on + opts.lab.t_off_ms/1e3;
        ifar = k_on - 12*cyc : k_on - 6*cyc - 1; inear = k_on - 2*cyc - 20 : k_on - 21; ifl = k_on + 2*cyc : k_on + 6*cyc - 1;
        %% ---- severidad e impedancia aparente (por registro)
        q = struct(); q.registro = name; q.conjunto = sets{is}; q.tipo = tipo; q.verdad = truth; q.en_paper = paper12(sets{is}, tipo);
        q.k_ref = k_ref; q.k_chk = k_chk; q.k_cd = k_cd; q.k_on = k_on; q.t_on_s = t_on; q.incepcion_respaldo = fallback;
        q.dif_ref_chk_ms = 1e3*(k_ref - k_chk)/fs; q.dur_falta_ciclos = ncy;
        [q.maniobra_previa_ms, q.maniobra_previa_pico_rel] = prestep(r.sig.box2, k_on, fs, f);
        q.P_pre_kW = mean(r.Psi(ifar)); q.I_pre_rms_A = mean(rms(r.Ic(ifar,:))); q.Vdc_pre_V = mean(r.Vdc(ifar));
        q.I_pico_falta_A = max(abs(r.Ic(k_on : k_on + 8*cyc, :)), [], 'all');
        for s = {'box2', 'conv1_ref_red'}
            V = r.sig.(s{1}); px = s{1};
            [Sf0, ph0] = seqph(V(ifar,:), t(ifar), f, A); Sn0 = seqph(V(inear,:), t(inear), f, A); [Sf, ph] = seqph(V(ifl,:), t(ifl), f, A);
            q.([px '_V1_pre_pu']) = abs(Sf0(2)); q.([px '_V2V1_pre']) = abs(Sf0(3)/Sf0(2));
            q.([px '_V2_pre_pu']) = abs(Sf0(3)); q.([px '_V2_falta_pu']) = abs(Sf(3)); q.([px '_V2_falta_menos_pre_pu']) = abs(Sf(3) - Sf0(3));
            q.([px '_paso_previo_pct']) = 100*(abs(Sn0(2))/abs(Sf0(2)) - 1);
            q.([px '_V1_rel']) = abs(Sf(2))/abs(Sf0(2)); q.([px '_V2V1']) = abs(Sf(3)/Sf(2)); q.([px '_argV2V1_deg']) = angle(Sf(3)/Sf(2))*180/pi;
            q.([px '_fase_min_rel']) = min(abs(ph)./abs(ph0));
            if isKey(argnom, truth)
                pr = pairix(truth); q.([px '_lazo_rel']) = abs(ph(pr(1)) - ph(pr(2)))/abs(ph0(pr(1)) - ph0(pr(2)));
                V2r = Sf(3)*exp(-1j*deg2rad(argnom(truth)));              % V2 referida al par en falta
                q.([px '_arg_residual_deg']) = angle(V2r/Sf(2))*180/pi;  % 0 en una falta franca
                z = Sf(2)/V2r - 1;                                        % Zf/Z2 si la red inversa es pasiva
            else
                q.([px '_lazo_rel']) = mean(abs(ph)./abs(ph0)); q.([px '_arg_residual_deg']) = NaN;
                rho = Sf(2)/Sn0(2); z = rho/(1 - rho);                    % Zf/Zs por fase (divisor de tension)
            end
            q.([px '_Zf_sobre_Zs_abs']) = abs(z); q.([px '_Zf_sobre_Zs_arg_deg']) = angle(z)*180/pi;
        end
        % canales de corriente: que miden ICH2/ICH3 y la corriente de meas_box2
        q.ICH2_pre_A = mean(r.ICH(ifar,2)); q.ICH3_pre_A = mean(r.ICH(ifar,3)); q.ICH2_falta_media_A = mean(r.ICH(ifl,2));
        pI = ph50(r.ICH(ifl,:), t(ifl), f); q.ICH2_falta_50Hz_A = abs(pI(2)); q.ICH3_falta_50Hz_A = abs(pI(3)); q.ICH_otros_max_A = max(abs(r.ICH(ifl,[1 4 5 6])), [], 'all');
        q.Pdc_pre_kW = mean(r.Vdc(ifar).*r.ICH(ifar,2))/1e3; q.Pdc_falta_kW = mean(r.Vdc(ifl).*r.ICH(ifl,2))/1e3; q.Pac_falta_kW = mean(r.Psi(ifl));
        Ipred = (r.Mv'*r.Ic')'; g = sum(r.Ib(ifar,:).*Ipred(ifar,:), 'all')/sum(Ipred(ifar,:).^2, 'all');
        pb = ph50(r.Ib(ifl,:), t(ifl), f); pd = ph50(r.Ib(ifl,:) - g*Ipred(ifl,:), t(ifl), f); pc = ph50(Ipred(ifl,:), t(ifl), f);
        q.Ibox2_falta_max_A = max(abs(pb)); q.Iconv_ref_falta_max_A = max(abs(pc)); q.Ibox2_menos_conv_ref_max_A = max(abs(pd)); q.ganancia_conv_ref = g;
        sev{end+1} = q; %#ok<AGROW>
        %% ---- clasificacion con la cadena, tres senales
        for s = 1:numel(signals)
            V = r.sig.(signals{s}); scan = [];
            [Sf0, ~] = seqph(V(ifar,:), t(ifar), f, A); Sf = seqph(V(ifl,:), t(ifl), f, A);
            for g2 = 1:2
                args = opts.chain; args.gate = gates{g2}; args.truth = truth; args.t_on = t_on; args.t_off = t_off;
                if g2 == 2 && isfield(opts.chain, 'gate') && strcmp(opts.chain.gate, 'none'), continue, end
                if isempty(scan)
                    args.fs = fs; args.t_range = [max(0, t_on - opts.lab.pre_s), t_off + 0.1];
                    out = r2_decision_chain(V, t, args); scan = out.scan;
                else
                    out = r2_decision_chain(scan, [], args);
                end
                E = out.event; W = out.win;
                full = W.t_start >= t_on - tol & W.t_end <= t_off + tol;
                pre = W.t_end <= t_on + tol; pfar = W.t_end <= t_on - 0.1 + tol & W.armed; pnear = pre & ~(W.t_end <= t_on - 0.1 + tol);
                w = struct(); w.registro = name; w.conjunto = sets{is}; w.tipo = tipo; w.verdad = truth; w.en_paper = q.en_paper;
                w.senal = signals{s}; w.puerta = gates{g2}; w.t_on_s = t_on; w.incepcion_respaldo = fallback; w.P_pre_kW = q.P_pre_kW;
                w.resultado = E.outcome; w.clase_declarada = E.decl_class; w.decl_ms = E.decl_ms; w.primera_correcta_ms = E.first_label_ms;
                w.primera_correcta_cruda_ms = E.first_label_raw_ms; w.inicio_racha_ms = E.decl_ms - 1e3*(out.opts.M - 1)*out.scan.stride/fs; w.racha_correcta_max = E.max_run_correct; w.retencion_pct = E.retention_pct;
                w.n_decl_erroneas = E.n_decl_wrong; w.clases_erroneas = E.wrong_classes; w.n_decl_prefalta = E.n_decl_pre; w.clases_prefalta = E.pre_classes;
                dpre = out.decl(out.decl.t_decl <= t_on + tol, :);
                if isempty(dpre), w.decl_prefalta_ms = NaN; else, w.decl_prefalta_ms = 1e3*(dpre.t_decl(1) - t_on); end
                w.n_ventanas_falta = nnz(full); lab = W.label_raw(full); w.etiquetas_falta = label_counts(lab);
                w.cambios_etiqueta = nnz(~strcmp(lab(1:end-1), lab(2:end)));
                th = W.angle_deg(full); thm = mod(angle(mean(exp(2j*deg2rad(th))))/2*180/pi, 180); dth = wrap90(th - thm);
                w.theta_mediana_deg = thm; w.theta_min_deg = thm + min(dth); w.theta_max_deg = thm + max(dth);
                if isKey(centre, truth)
                    mg = min(abs(wrap90(th - [15 75 135])), [], 2); dc = wrap90(th - centre(truth));
                    w.desv_centro_mediana_deg = median(dc); w.dist_centro_mediana_deg = median(abs(dc)); w.margen_frontera_min_deg = min(mg); w.margen_frontera_mediana_deg = median(mg);
                else
                    w.desv_centro_mediana_deg = NaN; w.dist_centro_mediana_deg = NaN; w.margen_frontera_min_deg = NaN; w.margen_frontera_mediana_deg = NaN;
                end
                w.redondez_mediana = median(W.roundness(full)); w.redondez_max = max(W.roundness(full));
                w.semieje_mayor_rel_mediana = median(W.semi_major(full))/(sqrt(3)/sqrt(2));
                jl = find(pre, 1, 'last'); w.umbral_redondez = W.thr(jl); w.r_ambiente = W.r_amb(jl);
                w.margen_redondez = w.umbral_redondez - w.redondez_max;
                w.V2V1_senal = abs(Sf(3)/Sf(2)); w.V1_rel_senal = abs(Sf(2))/abs(Sf0(2));
                w.V2V1_suelo = (1 - w.umbral_redondez)/(1 + w.umbral_redondez);     % redondez = (1-k)/(1+k)
                w.pinmax_max = max(W.pinmax(pre | full));
                w.pre_lejana_n = nnz(pfar); w.pre_lejana_no_sana_cruda_pct = 100*mean(~strcmp(W.label_raw(pfar), 'A-B-C')); w.pre_lejana_no_sana_pct = 100*mean(~strcmp(W.label(pfar), 'A-B-C'));
                w.pre_100ms_n = nnz(pnear); w.pre_100ms_no_sana_cruda_pct = 100*mean(~strcmp(W.label_raw(pnear), 'A-B-C')); w.pre_100ms_no_sana_pct = 100*mean(~strcmp(W.label(pnear), 'A-B-C'));
                w.pre_100ms_etiquetas = label_counts(W.label(pnear)); w.pre_redondez_min = min(W.roundness(pre & W.armed));
                rows{end+1} = w; %#ok<AGROW>
            end
        end
        fprintf('  %s %s\n', sets{is}, name);
    end
end
T = tostr(struct2table([rows{:}]')); S = tostr(struct2table([sev{:}]')); N = tostr(struct2table([nof{:}]'));
writetable(T, fullfile(outdir, 'lab_por_registro.csv')); writetable(S, fullfile(outdir, 'lab_severidad.csv')); writetable(N, fullfile(outdir, 'lab_sin_falta.csv'));

%% ------------------------------------------------------------ incepcion
P('\nINCEPCION. Registros con falta: %d. r2_lab_inception (conv1) frente al contraste en meas_box2:\n', height(S));
P('  diferencia (ms): mediana %.2f, min %.2f, max %.2f; registros con respaldo (|dif| > 1 ms): %d\n', median(S.dif_ref_chk_ms, 'omitnan'), min(S.dif_ref_chk_ms), max(S.dif_ref_chk_ms), nnz(S.incepcion_respaldo));
for i = find(S.incepcion_respaldo)'
    P('     %-13s %-24s k_ref = %d, k_chk = %d (%.1f ms antes de la falta real): se usa k_chk\n', S.conjunto(i), S.registro(i), S.k_ref(i), S.k_chk(i), -S.dif_ref_chk_ms(i));
end
P('  duracion de la falta (ciclos con firma de falta): de %d a %d\n', min(S.dur_falta_ciclos), max(S.dur_falta_ciclos));
for is = 1:2
    x = S(S.conjunto == sets{is}, :);
    P('  %-13s maniobra previa detectada en %d de %d registros; instante de %.1f a %.1f ms respecto a la falta; transitorio de hasta %.3f del pico; escalon de |V1| de %+.2f a %+.2f %%\n', ...
        sets{is}, nnz(~isnan(x.maniobra_previa_ms)), height(x), min(x.maniobra_previa_ms), max(x.maniobra_previa_ms), max(x.maniobra_previa_pico_rel), min(x.box2_paso_previo_pct), max(x.box2_paso_previo_pct));
end

%% ------------------------------------------------------------- severidad
P('\nSEVERIDAD MEDIDA en meas_box2, regimen de falta (ciclos 3 a 6), por conjunto y tipo [min, max]:\n');
fid = fopen(fullfile(outdir, 'lab_severidad_resumen.csv'), 'w');
fprintf(fid, ['conjunto,tipo,n,P_pre_kW_min,P_pre_kW_max,paso_previo_pct_min,paso_previo_pct_max,V2V1_min,V2V1_max,lazo_rel_min,lazo_rel_max,V1_rel_min,V1_rel_max,' ...
    'arg_residual_deg_min,arg_residual_deg_max,Zf_sobre_Zs_abs_min,Zf_sobre_Zs_abs_max,Zf_sobre_Zs_arg_deg_min,Zf_sobre_Zs_arg_deg_max,' ...
    'conv1ref_V2V1_min,conv1ref_V2V1_max,conv1ref_arg_residual_deg_min,conv1ref_arg_residual_deg_max\n']);
tg = {'ab', 'bc', 'ac', 'abc'};
for is = 1:2
    for q = 1:4
        x = S(S.conjunto == sets{is} & S.tipo == tg{q}, :); if isempty(x), continue, end
        mm = @(v) [min(v), max(v)];
        fprintf(fid, '%s,%s,%d,%.2f,%.2f,%.2f,%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.2f,%.2f,%.4f,%.4f,%.2f,%.2f,%.4f,%.4f,%.2f,%.2f\n', sets{is}, tg{q}, height(x), mm(x.P_pre_kW), ...
            mm(x.box2_paso_previo_pct), mm(x.box2_V2V1), mm(x.box2_lazo_rel), mm(x.box2_V1_rel), mm(x.box2_arg_residual_deg), mm(x.box2_Zf_sobre_Zs_abs), mm(x.box2_Zf_sobre_Zs_arg_deg), ...
            mm(x.conv1_ref_red_V2V1), mm(x.conv1_ref_red_arg_residual_deg));
        P('  %-13s %-3s n = %d | P previa %.1f a %.1f kW | paso previo %+.2f a %+.2f %% | V2/V1 %.4f a %.4f | lazo %.4f a %.4f | V1 %.4f a %.4f | arg residual %+.2f a %+.2f | Zf/Zs %.3f a %.3f, arg %+.1f a %+.1f\n', ...
            sets{is}, tg{q}, height(x), mm(x.P_pre_kW), mm(x.box2_paso_previo_pct), mm(x.box2_V2V1), mm(x.box2_lazo_rel), mm(x.box2_V1_rel), mm(x.box2_arg_residual_deg), ...
            mm(x.box2_Zf_sobre_Zs_abs), mm(x.box2_Zf_sobre_Zs_arg_deg));
    end
end
fclose(fid);
% canales de corriente
c1 = corrcoef(S.Pdc_pre_kW, S.P_pre_kW); pf = polyfit(S.P_pre_kW, S.Pdc_pre_kW, 1);
P('\nCANALES DE CORRIENTE:\n');
P('  Egston ICH2 e ICH3 son corriente continua: ICH2 previa de %.2f a %.2f A, ICH3 = -ICH2 (suma maxima %.2f A); componente de 50 Hz en falta %.2f A como maximo.\n', ...
    min(S.ICH2_pre_A), max(S.ICH2_pre_A), max(abs(S.ICH2_pre_A + S.ICH3_pre_A)), max([S.ICH2_falta_50Hz_A; S.ICH3_falta_50Hz_A]));
P('  Vdc*ICH2 frente a la potencia alterna del convertidor (P_si) en prefalta, %d registros: r = %.5f, Pdc = %.3f*P + %.3f kW. En falta: Pdc de %.2f a %.2f kW, P de %.2f a %.2f kW.\n', ...
    height(S), c1(1,2), pf(1), pf(2), min(S.Pdc_falta_kW), max(S.Pdc_falta_kW), min(S.Pac_falta_kW), max(S.Pac_falta_kW));
P('  => ICH2/ICH3 son la alimentacion de continua del convertidor, no la rama de falta. Los otros cuatro canales ICH no pasan de %.2f A.\n', max(S.ICH_otros_max_A));
P('  Corriente de meas_box2 frente a la de conv1 referida al lado de red (Mv''*i), en falta: |I box2| hasta %.1f A, diferencia fasorial maxima %.1f A (mediana %.1f A).\n', ...
    max(S.Ibox2_falta_max_A), max(S.Ibox2_menos_conv_ref_max_A), median(S.Ibox2_menos_conv_ref_max_A));
P('  => meas_box2 mide la corriente hacia el transformador y el convertidor, no la de la rama de falta. Ningun canal registra la corriente de falta:\n');
P('     la impedancia de falta en ohmios no se puede estimar. Solo queda el cociente aparente Zf/Zs deducido de las tensiones.\n');
writetable(S(:, {'registro','conjunto','tipo','P_pre_kW','Pdc_pre_kW','ICH2_pre_A','ICH3_pre_A','ICH2_falta_media_A','ICH2_falta_50Hz_A','ICH3_falta_50Hz_A','ICH_otros_max_A', ...
    'Pdc_falta_kW','Pac_falta_kW','Ibox2_falta_max_A','Iconv_ref_falta_max_A','Ibox2_menos_conv_ref_max_A','ganancia_conv_ref'}), fullfile(outdir, 'lab_canales_corriente.csv'));

%% --------------------------------------------------------------- resumen
fid = fopen(fullfile(outdir, 'lab_resumen.csv'), 'w');
fprintf(fid, ['grupo,tipo,senal,puerta,n,correctas,erroneas,sin_declaracion,clases_declaradas,primera_correcta_ms_mediana,primera_correcta_ms_p90,primera_correcta_ms_min,primera_correcta_ms_max,' ...
    'decl_ms_mediana,decl_ms_p90,decl_ms_min,decl_ms_max,inicio_racha_ms_mediana,inicio_racha_ms_max,desv_centro_mediana_deg,desv_centro_min_deg,desv_centro_max_deg,margen_frontera_min_deg,' ...
    'redondez_mediana,redondez_max,umbral_redondez_min,margen_redondez_min,V2V1_min,V2V1_suelo_max,registros_con_decl_prefalta,pre_lejana_no_sana_pct,pre_100ms_no_sana_pct,ventanas_pre_lejana\n']);
grp = {'paper12', T.en_paper; 'Outkontaktor', T.conjunto == "Outkontaktor"; 'Autokontaktor', T.conjunto == "Autokontaktor"; 'todos', true(height(T), 1)};
tys = {'LL', T.tipo ~= "abc"; 'abc', T.tipo == "abc"};
P('\nRESULTADO POR EVENTO:\n  %-13s %-3s %-13s %-8s | %2s %4s %4s %4s | 1a correcta ms med/p90/max | declaracion ms med/p90/max | desv. centro med [min, max] | margen frontera min | redondez max / umbral min | decl. prefalta\n', ...
    'grupo', 'tip', 'senal', 'puerta', 'n', 'ok', 'err', 'sin');
for ig = 1:size(grp, 1)
    for it = 1:2
        for s = 1:numel(signals)
            for g2 = 1:2
                x = T(grp{ig,2} & tys{it,2} & T.senal == signals{s} & T.puerta == gates{g2}, :); if isempty(x), continue, end
                ok = x.resultado == "correcta"; wn = @(p, w) sum(p.*w)/sum(w);
                fprintf(fid, '%s,%s,%s,%s,%d,%d,%d,%d,%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.5f,%d,%.3f,%.3f,%d\n', grp{ig,1}, tys{it,1}, signals{s}, gates{g2}, ...
                    height(x), nnz(ok), nnz(x.resultado == "erronea"), nnz(x.resultado == "sin_declaracion"), label_counts(cellstr(x.verdad + ">" + x.clase_declarada)), ...
                    median(x.primera_correcta_ms, 'omitnan'), prc(x.primera_correcta_ms, 90), min(x.primera_correcta_ms), max(x.primera_correcta_ms), ...
                    median(x.decl_ms(ok)), prc(x.decl_ms(ok), 90), min([x.decl_ms(ok); NaN]), max([x.decl_ms(ok); NaN]), median(x.inicio_racha_ms(ok)), max([x.inicio_racha_ms(ok); NaN]), median(x.desv_centro_mediana_deg), min(x.desv_centro_mediana_deg), max(x.desv_centro_mediana_deg), ...
                    min(x.margen_frontera_min_deg), median(x.redondez_mediana), max(x.redondez_max), min(x.umbral_redondez), min(x.margen_redondez), min(x.V2V1_senal), max(x.V2V1_suelo), ...
                    nnz(x.n_decl_prefalta > 0), wn(x.pre_lejana_no_sana_pct, x.pre_lejana_n), wn(x.pre_100ms_no_sana_pct, x.pre_100ms_n), sum(x.pre_lejana_n));
                P('  %-13s %-3s %-13s %-8s | %2d %4d %4d %4d | %5.2f %5.2f %5.2f | %5.2f %5.2f %5.2f | %+5.1f [%+5.1f, %+5.1f] | %5.1f | %.3f / %.3f | %d\n', grp{ig,1}, tys{it,1}, signals{s}, gates{g2}, ...
                    height(x), nnz(ok), nnz(x.resultado == "erronea"), nnz(x.resultado == "sin_declaracion"), median(x.primera_correcta_ms, 'omitnan'), prc(x.primera_correcta_ms, 90), max(x.primera_correcta_ms), ...
                    median(x.decl_ms(ok)), prc(x.decl_ms(ok), 90), max([x.decl_ms(ok); NaN]), median(x.desv_centro_mediana_deg), min(x.desv_centro_mediana_deg), max(x.desv_centro_mediana_deg), ...
                    min(x.margen_frontera_min_deg), max(x.redondez_max), min(x.umbral_redondez), nnz(x.n_decl_prefalta > 0));
            end
        end
    end
end
fclose(fid);
P('\nREGISTROS SIN FALTA (todo el registro):\n');
for i = 1:height(N)
    P('  %-22s %-13s %.1f s, %d ventanas: etiqueta cruda no sana %.2f %%, final %.2f %%, declaraciones %d %s | V2/V1 %.4f, redondez minima %.4f, umbral mediano %.4f\n', ...
        N.registro(i), N.senal(i), N.dur_s(i), N.n_win(i), N.raw_nonhealthy_pct(i), N.nonhealthy_pct(i), N.n_decl(i), N.decl_classes(i), N.V2V1(i), N.roundness_min(i), N.thr_med(i));
end

%% --------------------------------------------- paridad con T12 (12 L-L)
t12 = fullfile(RAIZ, '_revision_TPWRD', 'ronda2', 'resultados', 'T12', 'lab_dy_classification.csv');
if isfile(t12) && isempty(fieldnames(opts.chain))
    Q = readtable(t12, 'TextType', 'string'); Q.senal(Q.senal == "box2_fase") = "box2";
    fid = fopen(fullfile(outdir, 'lab_paridad_T12.csv'), 'w');
    fprintf(fid, 'registro,senal,puerta,resultado_igual,clase_igual,dif_decl_ms,dif_primera_ms,dif_theta_deg,dif_margen_min_deg\n'); nbad = 0; mx = [0 0 0 0];
    for i = 1:height(Q)
        j = find(T.registro == Q.registro(i) & T.senal == Q.senal(i) & T.puerta == Q.puerta(i) & T.conjunto == "Outkontaktor");
        if isempty(j), continue, end
        dd = [T.decl_ms(j) - Q.decl_ms(i), T.primera_correcta_ms(j) - Q.primera_correcta_ms(i), T.theta_mediana_deg(j) - Q.theta_mediana_deg(i), T.margen_frontera_min_deg(j) - Q.margen_frontera_min_deg(i)];
        e1 = T.resultado(j) == Q.resultado(i); e2 = T.clase_declarada(j) == Q.clase_declarada(i);
        fprintf(fid, '%s,%s,%s,%d,%d,%.3f,%.3f,%.3f,%.3f\n', Q.registro(i), Q.senal(i), Q.puerta(i), e1, e2, dd);
        nbad = nbad + ~(e1 && e2 && abs(dd(1)) < 0.006 && abs(dd(2)) < 0.006); mx = max(mx, abs(dd));
    end
    fclose(fid);
    P('\nPARIDAD CON T12 (lab_dy_classification.csv, %d filas): %d filas distintas; diferencia maxima de declaracion %.3f ms, de primera etiqueta %.3f ms, de angulo %.2f grados, de margen %.2f grados.\n', height(Q), nbad, mx);
end
md5_fin = r2_t20_md5();
if strcmp(md5_ini, md5_fin), P('\nVersion de la cadena al terminar: la misma que al empezar.\n');
else, P('\nAVISO: la cadena cambio en disco durante la ejecucion. MD5 al terminar: %s\n', md5_fin); end
fclose(logf);
end

%% ======================================================================
function [k_chk, k_cd, ncy] = incep_check(V, t, fs, f, A)
% Segundo estimador de la incepcion, independiente de r2_lab_inception y robusto a la
% maniobra previa de los registros Autokontaktor (escalon de tension de un 2 % entre
% 12 y 47 ms antes de la falta, con un transitorio breve que supera el 5 % del pico).
% 1) ciclo con firma de falta: |V2/V1| > 0.3 o |V1| < 0.6 del nivel mediano;
% 2) desde el ciclo anterior, primera muestra en que la tension se aparta mas del
%    25 % del pico de su valor un ciclo antes, sostenido 5 muestras (k_cd);
% 3) afinado con el criterio del paper (5 % respecto a la extrapolacion sinusoidal),
%    ajustando la senoide en los dos ciclos que acaban 2 ms antes de k_cd y probando
%    solo los 4 ms anteriores a k_cd.
cyc = round(fs/f); nw = floor(numel(t)/cyc); v1 = zeros(nw, 1); k2 = zeros(nw, 1);
for i = 1:nw
    ix = (i-1)*cyc + 1 : i*cyc; S = seqph(V(ix,:), t(ix), f, A); v1(i) = abs(S(2)); k2(i) = abs(S(3)/S(2));
end
flt = k2 > 0.3 | v1 < 0.6*median(v1); c = find(flt, 1); ncy = nnz(flt);
k_chk = NaN; k_cd = NaN; if isempty(c) || c < 6, return, end
Vpk = median(v1);
d = [zeros(cyc, 3); abs(V(cyc+1:end,:) - V(1:end-cyc,:))]; m = max(d, [], 2)/Vpk;
i0 = (c-2)*cyc + 1; hit = movmin(m(i0:end), [0 4]) > 0.25; k_cd = i0 - 1 + find(hit, 1);
ifit = k_cd - 2*cyc - 20 : k_cd - 21; X = [cos(2*pi*f*t(ifit)), sin(2*pi*f*t(ifit)), ones(numel(ifit), 1)]; cf = X \ V(ifit,:);
it = k_cd - 40 : k_cd + 10; Xt = [cos(2*pi*f*t(it)), sin(2*pi*f*t(it)), ones(numel(it), 1)];
res = max(abs(V(it,:) - Xt*cf), [], 2)/max(abs(V(ifit,:)), [], 'all');
h2 = movmin(res, [0 4]) > 0.05; j = find(h2, 1);
if isempty(j), k_chk = k_cd; else, k_chk = it(j); end
end

function [t_ms, pk] = prestep(V, k_on, fs, f)
% Maniobra previa a la falta: primera muestra, entre 120 y 3 ms antes de k_on, en que
% la diferencia con el ciclo anterior supera 3 veces su maximo en los 120 ms previos
% (sostenido 3 muestras). Devuelve el instante relativo (ms) y el pico de esa diferencia
% (relativo al pico de tension) hasta 3 ms antes de la falta.
cyc = round(fs/f); Vpk = max(abs(V(k_on - 12*cyc : k_on - 7*cyc, :)), [], 'all');
d = [zeros(cyc, 3); abs(V(cyc+1:end,:) - V(1:end-cyc,:))]; m = max(d, [], 2)/Vpk;
base = max(m(k_on - 2400 : k_on - 1201)); seg = k_on - 1200 : k_on - 30;
j = find(movmin(m(seg), [0 2]) > max(3*base, 0.02), 1); pk = max(m(seg));
if isempty(j), t_ms = NaN; else, t_ms = 1e3*(seg(j) - k_on)/fs; end
end

function [S, ph] = seqph(v, t, f, A)
X = [cos(2*pi*f*t), sin(2*pi*f*t), ones(numel(t), 1)];
c = X \ v; ph = (c(1,:) - 1j*c(2,:)).'; S = A*ph;
end
function p = ph50(x, t, f)
X = [cos(2*pi*f*t), sin(2*pi*f*t), ones(numel(t), 1)]; c = X \ x; p = c(1,:) - 1j*c(2,:);
end
function s = label_counts(lab)
if isempty(lab), s = ''; return, end
lab = cellstr(lab); lab(cellfun(@isempty, lab)) = {'-'};
[u, ~, j] = unique(lab); c = accumarray(j(:), 1); [c, o] = sort(c, 'descend'); u = u(o);
s = strjoin(arrayfun(@(i) sprintf('%s:%d', u{i}, c(i)), 1:numel(u), 'UniformOutput', false), ' ');
end
function T = tostr(T)
for c = string(T.Properties.VariableNames)
    if iscell(T.(c)), T.(c) = string(T.(c)); end
end
end
function p = prc(x, q)
x = x(isfinite(x)); if isempty(x), p = NaN; else, p = prctile(x, q); end
end
function y = wrap90(x), y = mod(x + 90, 180) - 90; end
function plog(fid, varargin), fprintf(varargin{:}); fprintf(fid, varargin{:}); end
