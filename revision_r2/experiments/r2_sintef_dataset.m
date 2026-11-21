function T = r2_sintef_dataset(opts)
% R2_SINTEF_DATASET  (T20; peticiones R1.2, R1.1, E1, R4.2a, E2)
% Pasa la cadena de decision de referencia (r2_decision_chain, T03) por todos los
% registros unicos del dataset EMT de SINTEF (Database_Santiago), en los dos buses
% de cada archivo (DG2 y Pex).
%
%   r2_sintef_dataset()            % cadena de referencia
%   r2_sintef_dataset(opts)        % opts = r2_t20_defaults(...); ver ese archivo
%
% Por registro y bus guarda: resultado por evento (declaracion, clase, latencias),
% fasores de secuencia en prefalta y en falta, V0 relativo, angulo de la elipse
% (clasificador y prediccion 165 - arg(V2/V1)/2), redondez, valor de la puerta de
% secuencia cero, tasas de etiquetas en prefalta y tras el despeje, y la comparacion
% entre tres escalados (pu con la tension de prefalta, voltios, pu de base fija).
% Resumenes por factor, matriz de clases, angulo L-L por factor con diferencias
% pareadas, lectura de tierra por neutro y tasas de prefalta.
%
% Convenios (todos los instantes al final de la ventana):
%   t_on  = 1.000 s + retardo del registro (t(1); 0 salvo en dealyedLogging)
%   t_off = t_on + 0.100 s (despeje por temporizador; no hay marcador en los .mat)
%   prefalta: ventanas con t_end <= t_on; el registro entero se procesa desde t(1)
%   regimen de falta: de t_on + 40 ms a t_off (tres ultimos ciclos)
%   escala principal: cada registro se divide por el pico de secuencia directa de
%   los cuatro ciclos que acaban un ciclo antes de t_on (pu de prefalta)
% Solo lectura sobre el dataset. No se guardan ondas, solo magnitudes derivadas.
% Con opts.reuse_scan = true guarda (y reutiliza) la clasificacion por ventana de la
% escala principal en sintef_scans_cache.mat, para probar otra puerta de redondez u
% otra persistencia sin reclasificar. Es una cache: no va al repositorio.
% Salidas en resultados\T20\ (o resultados\T20\<tag>\).
if nargin < 1 || isempty(opts), opts = r2_t20_defaults(); end
here = fileparts(mfilename('fullpath')); addpath(here);
outdir = opts.outdir; BASE = opts.sintef.base;
logf = fopen(fullfile(outdir, 'sintef_log.txt'), 'w');
P = @(varargin) plog(logf, varargin{:});
P('r2_sintef_dataset %s MATLAB %s tag "%s"\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), version, opts.tag);
P('Opciones de la cadena distintas de las de referencia: %s\n', jsonencode(opts.chain));
md5_ini = r2_t20_md5(); P('Version de la cadena (MD5): %s\n', md5_ini);
fs = 4000; f = 50; spc = fs/f; tol = 0.25/fs;
a = exp(2j*pi/3); A = [1 1 1; 1 a a^2; 1 a^2 a]/3;
GROUND = {'A-G','B-G','C-G','A-B-G','B-C-G','C-A-G'};
centre = containers.Map({'A-B','B-C','C-A'}, {45, 165, 105});
scales = opts.sintef.scales; nsc = numel(scales);

%% ------------------------------------------------------------ 1. inventario
all = dir(fullfile(BASE, '**', '*.mat'));
[~, o] = sort(strcat({all.folder}, filesep, {all.name})); all = all(o);
n = numel(all); rel = cell(n,1); md = cell(n,1); isorig = false(n,1); score = zeros(n,1); L = cell(n,1);
for k = 1:n
    rel{k} = erase(all(k).folder, [BASE filesep]);
    isorig(k) = contains(lower(rel{k}), 'original');
    fid = fopen(fullfile(all(k).folder, all(k).name), 'r'); b = fread(fid, Inf, '*uint8'); fclose(fid);
    h = java.security.MessageDigest.getInstance('MD5'); h.update(b);
    md{k} = sprintf('%02x', typecast(h.digest(), 'uint8'));
    L{k} = parse_label(all(k).name, rel{k});
    score(k) = L{k}.folder_match;
end
keep = ~isorig; reason = repmat({''}, n, 1); reason(isorig) = {'carpeta original (version previa sustituida)'};
[um, ~, im] = unique(md);
for g = 1:numel(um)
    idx = find(im == g & ~isorig);
    if numel(idx) > 1
        [~, b] = max(score(idx)); home = idx(b);
        for j = setdiff(idx(:)', home)
            keep(j) = false;
            reason{j} = sprintf('duplicado exacto (MD5) de %s', fullfile(rel{home}, all(home).name));
        end
    end
end
fid = fopen(fullfile(outdir, 'sintef_exclusiones.csv'), 'w');
fprintf(fid, 'carpeta,archivo,md5,motivo\n');
for k = find(~keep)', fprintf(fid, '%s,%s,%s,%s\n', rel{k}, all(k).name, md{k}, reason{k}); end
fclose(fid);
recs = find(keep);
P('Archivos .mat: %d. Excluidos: %d en carpetas original, %d duplicados exactos. Registros unicos: %d (%d ejecuciones con los dos buses).\n', ...
    n, nnz(isorig), nnz(~keep & ~isorig), numel(recs), 2*numel(recs));
bad = recs(score(recs) < 4);
if ~isempty(bad), P('AVISO: %d registros cuyo nombre no coincide con su carpeta (se usa el nombre).\n', numel(bad)); end
if ~isempty(opts.sintef.pick), recs = recs(ismember({all(recs).name}, opts.sintef.pick)); end
if ~isempty(opts.sintef.limit), recs = recs(round(linspace(1, numel(recs), min(numel(recs), opts.sintef.limit)))); end

%% ------------------------------------------- 2. escaneos guardados (opcional)
scanfile = fullfile(outdir, 'sintef_scans_cache.mat'); SC = struct();   % cache opcional (unos 95 MB), regenerable
ids = strcat(rel(recs), filesep, {all(recs).name}');
if opts.reuse_scan && isfile(scanfile)
    q = load(scanfile);
    if isequaln(q.key, opts.scan_key) && isequal(q.ids, ids) && isequal(q.scales, scales)
        SC = q.SC; P('Se reutiliza la clasificacion por ventana de %s\n', scanfile);
    else
        P('La clasificacion guardada no corresponde a estas opciones: se reclasifica.\n');
    end
end
have = isfield(SC, 'scan'); if ~have, SC.scan = cell(numel(recs), 2); end
buses = {'DG2', 'Pex'};

%% ------------------------------------------------------------ 3. registros
rows = cell(2*numel(recs), 1); nr = 0; tic;
for ir = 1:numel(recs)
    k = recs(ir); lab = L{k};
    d = load(fullfile(all(k).folder, all(k).name)); G = d.GridConvs;
    for ib = 1:2
        bus = buses{ib};
        if ib == 1, ts = G.DG2.Vac; else, ts = G.Pex.Vabc; end
        v = double(ts.Data); t = double(ts.Time(:));
        t_on = opts.sintef.t_fault + t(1); t_off = t_on + opts.sintef.dur;
        k_on = find(t <= t_on + tol, 1, 'last');
        ipre = k_on - 5*spc + 1 : k_on - spc; ifl = k_on + 1 : k_on + 5*spc; ist = k_on + 2*spc + 1 : k_on + 5*spc;
        [S0, ph0] = seqph(v(ipre,:), t(ipre), f, A); S = seqph(v(ist,:), t(ist), f, A); Sf = seqph(v(ifl,:), t(ifl), f, A);
        V1pre = abs(S0(2));
        r = struct();
        r.registro = ids{ir}; r.archivo = all(k).name; r.bus = bus; r.cat = lab.cat; r.faulted = lab.faulted; r.truth = lab.truth;
        r.Srated = lab.Srated; r.pu = lab.pu; r.ground = lab.ground; r.limit = lab.limit; r.variant = lab.variant;
        r.N = numel(t); r.t0 = t(1); r.t_on = t_on; r.t_off = t_off;
        % incepcion estimada (solo informativa): desviacion > 5 % respecto al ciclo anterior
        dcy = [zeros(spc,3); abs(v(spc+1:end,:) - v(1:end-spc,:))];
        ki = find(max(dcy, [], 2)/V1pre > 0.05 & t > t_on - 0.05, 1);
        if isempty(ki), r.t_inc_est = NaN; else, r.t_inc_est = t(ki); end
        % tension de la fase en falta en el instante de incepcion (punto de onda)
        fa = lab.faulted - 'A' + 1;
        if numel(fa) == 1, vref = v(k_on, fa); else, vref = v(k_on, fa(1)) - v(k_on, fa(2)); vref = vref/sqrt(3); end
        r.pow_v_pu = vref/V1pre;
        % fasores
        r.V1pre_V = V1pre; r.V2pre_rel = abs(S0(3))/V1pre; r.V0pre_rel = abs(S0(1))/V1pre;
        r.V0_rel = abs(S(1))/V1pre; r.V0_rel_full = abs(Sf(1))/V1pre; r.V1_rel = abs(S(2))/V1pre; r.V2_rel = abs(S(3))/V1pre;
        r.V2V1 = abs(S(3)/S(2)); r.argV2V1_deg = angle(S(3)/S(2))*180/pi;
        r.theta_pred_deg = mod(165 - 0.5*r.argV2V1_deg, 180);
        r.roundness_pred = abs(1 - r.V2V1)/(1 + r.V2V1);            % redondez de la elipse de la componente fundamental
        Xf = [cos(2*pi*f*t(ist)), sin(2*pi*f*t(ist)), ones(numel(ist), 1)]; ef = v(ist,:) - Xf*(Xf\v(ist,:));
        r.resid50_rel = sqrt(mean(ef(:).^2))/(V1pre/sqrt(2));          % contenido no fundamental en regimen de falta (rms, pu)
        Xp = [cos(2*pi*f*t(ipre)), sin(2*pi*f*t(ipre)), ones(numel(ipre), 1)]; ep = v(ipre,:) - Xp*(Xp\v(ipre,:));
        r.resid50_pre_rel = sqrt(mean(ep(:).^2))/(V1pre/sqrt(2));
        if numel(fa) == 1
            r.argV0_rel_faulted_deg = angle(-S(1)/ph0(fa))*180/pi;      % 0 = V0 opuesta a la tension de prefalta de la fase en falta
        else
            r.argV0_rel_faulted_deg = NaN;
        end
        rf = rms(v(ist,:))./rms(v(ipre,:));
        oth = setdiff(1:3, fa);
        r.vratio_faulted_min = min(rf(fa)); r.vratio_faulted_max = max(rf(fa)); r.vratio_others_max = max(rf(oth)); r.vratio_others_min = min(rf(oth));
        [~, jm] = min(rf); r.lowest_phase = char('A' + jm - 1);
        % ---- cadena, tres escalados
        scv = [V1pre, 1, opts.sintef.v_nom]; names = {'pu_pre', 'volts', 'pu_nom'};
        W1 = [];
        for is = 1:nsc
            sc = scv(strcmp(names, scales{is}));
            args = opts.chain; args.truth = lab.truth; args.t_on = t_on; args.t_off = t_off;
            if have && is == 1 && ~isempty(SC.scan{ir, ib})
                scan = SC.scan{ir, ib}; scan.label_raw = cellstr(scan.label_raw);
                out = r2_decision_chain(scan, [], args);
            else
                args.fs = fs; args.scale = sc;
                out = r2_decision_chain(v, t, args);
                scan = out.scan; scan.label_raw = categorical(scan.label_raw);
                if is == 1 && opts.reuse_scan, SC.scan{ir, ib} = scan; end
            end
            E = out.event; W = out.win; sx = scales{is};
            if is == 1
                W1 = W; r.scale_V = sc;
                r = add_event(r, E, '');
                full = W.t_start >= t_on - tol & W.t_end <= t_off + tol;
                st = W.t_start >= t_on + 0.04 - tol & W.t_end <= t_off + tol;
                pre = W.t_end <= t_on + tol; prel = pre & W.t_start >= t(1) + 0.2 - tol;
                aft = W.t_start >= t_off - tol;
                isg = ismember(W.label, GROUND); isgraw = ismember(W.label_raw, GROUND);
                r.n_full = nnz(full); r.full_ground_pct = 100*mean(isg(full)); r.full_correct_pct = 100*mean(strcmp(W.label(full), lab.truth));
                r.full_healthy_pct = 100*mean(strcmp(W.label(full), 'A-B-C'));
                r.full_modal_label = E.full_modal_label;
                r.full_labels = label_counts(W.label(full));
                r.full_kirchhoff_fail_pct = 100*mean(~W.is_kirchhoff(full)); r.full_degenerate_pct = 100*mean(W.deg_decides(full));
                r.pinmax_full_med = median(W.pinmax(full)); r.pinmax_full_min = min(W.pinmax(full)); r.pinmax_full_max = max(W.pinmax(full));
                r.roundness_st_med = median(W.roundness(st)); r.roundness_st_max = max(W.roundness(st));
                r.semi_major_st_med = median(W.semi_major(st));
                % angulo de la elipse (ventanas de regimen que pasan la puerta de Kirchhoff)
                sa = st & W.is_kirchhoff & isfinite(W.angle_deg);
                [r.theta_st_deg, r.theta_st_iqr, r.theta_st_min, r.theta_st_max, r.n_theta_st] = axial_stats(W.angle_deg(sa));
                fa2 = full & W.is_kirchhoff & isfinite(W.angle_deg);
                [r.theta_full_deg, r.theta_full_iqr, r.theta_full_min, r.theta_full_max] = axial_stats(W.angle_deg(fa2));
                if isKey(centre, lab.truth)
                    c0 = centre(lab.truth);
                    r.dev_centre_st_deg = wrap90(r.theta_st_deg - c0); r.dev_centre_pred_deg = wrap90(r.theta_pred_deg - c0);
                    r.margin_min_st_deg = min([min(abs(wrap90(W.angle_deg(sa) - [15 75 135])), [], 2); NaN]);
                    r.margin_min_full_deg = min([min(abs(wrap90(W.angle_deg(fa2) - [15 75 135])), [], 2); NaN]);
                    r.pred_minus_cls_deg = wrap90(r.theta_pred_deg - r.theta_st_deg);
                else
                    r.dev_centre_st_deg = NaN; r.dev_centre_pred_deg = NaN; r.margin_min_st_deg = NaN; r.margin_min_full_deg = NaN; r.pred_minus_cls_deg = NaN;
                end
                % prefalta (1 s sano) y tras el despeje
                r.pre_n_win = nnz(pre); r.pre_n_init = nnz(pre & ~W.armed); r.pre_n_armed = nnz(pre & W.armed);
                r.pre_armed_nonhealthy_pct = 100*mean(~strcmp(W.label(pre & W.armed), 'A-B-C'));
                r.pre_raw_nonhealthy_pct = 100*mean(~strcmp(W.label_raw(pre), 'A-B-C')); r.pre_nonhealthy_pct = 100*mean(~strcmp(W.label(pre), 'A-B-C'));
                r.pre_raw_ground_pct = 100*mean(isgraw(pre)); r.pre_kirchhoff_fail_pct = 100*mean(~W.is_kirchhoff(pre));
                r.pre_late_n_win = nnz(prel); r.pre_late_raw_nonhealthy_pct = 100*mean(~strcmp(W.label_raw(prel), 'A-B-C'));
                r.pre_late_nonhealthy_pct = 100*mean(~strcmp(W.label(prel), 'A-B-C')); r.pre_late_raw_ground_pct = 100*mean(isgraw(prel));
                r.pre_labels = label_counts(W.label(pre & W.armed));
                r.pinmax_pre_med = median(W.pinmax(pre)); r.pinmax_pre_max = max(W.pinmax(pre)); r.pinmax_pre_late_max = max(W.pinmax(prel));
                r.roundness_pre_late_med = median(W.roundness(prel)); r.roundness_pre_late_min = min(W.roundness(prel));
                r.thr_at_fault = W.thr(find(pre, 1, 'last')); r.r_amb_at_fault = W.r_amb(find(pre, 1, 'last'));
                r.n_decl_pre = E.n_decl_pre; r.pre_classes = E.pre_classes;
                r.aft_n_win = nnz(aft); r.aft_nonhealthy_pct = 100*mean(~strcmp(W.label(aft), 'A-B-C')); r.aft_ground_pct = 100*mean(isg(aft));
                r.aft_pinmax_max = max(W.pinmax(aft)); r.n_decl_after = E.n_decl_after; r.after_classes = E.after_classes;
                % sin puerta de redondez (etiquetas crudas + persistencia)
                a2 = args; a2.gate = 'none'; a2.fs = fs; sc2 = scan; sc2.label_raw = cellstr(sc2.label_raw);
                o2 = r2_decision_chain(sc2, [], a2); E2 = o2.event;
                r.outcome_nogate = E2.outcome; r.decl_class_nogate = E2.decl_class; r.decl_ms_nogate = E2.decl_ms; r.n_decl_pre_nogate = E2.n_decl_pre;
                r.pre_classes_nogate = E2.pre_classes;
            else
                r.(['scale_' sx]) = sc;
                r.(['outcome_' sx]) = E.outcome; r.(['decl_class_' sx]) = E.decl_class; r.(['decl_ms_' sx]) = E.decl_ms;
                r.(['first_label_ms_' sx]) = E.first_label_ms; r.(['n_decl_pre_' sx]) = E.n_decl_pre; r.(['n_decl_after_' sx]) = E.n_decl_after;
                r.(['ndiff_raw_' sx]) = nnz(~strcmp(W.label_raw, W1.label_raw)); r.(['ndiff_final_' sx]) = nnz(~strcmp(W.label, W1.label));
                infault = W.t_end > t_on + tol & W.t_end <= t_off + tol;
                r.(['ndiff_final_fault_' sx]) = nnz(~strcmp(W.label(infault), W1.label(infault)));
                r.(['ndiff_line_' sx]) = nnz(W.is_line ~= W1.is_line);
            end
        end
        r.n_win = height(W1);
        nr = nr + 1; rows{nr} = r;
    end
    if mod(ir, 50) == 0, fprintf('  %d/%d registros, %.0f s\n', ir, numel(recs), toc); end
end
T = struct2table([rows{1:nr}]');
for c = string(T.Properties.VariableNames)
    if iscell(T.(c)), T.(c) = string(T.(c)); end
end
writetable(T, fullfile(outdir, 'sintef_por_registro.csv'));
key = opts.scan_key; %#ok<NASGU>
if opts.reuse_scan && ~have
    save(scanfile, 'SC', 'key', 'ids', 'scales', '-v7.3');
end
P('Ejecuciones procesadas: %d (%.0f s). Ventanas por registro: %d a %d.\n', height(T), toc, min(T.n_win), max(T.n_win));

%% ------------------------------------------------------- 4. comprobaciones
P('\nIncepcion estimada menos t_on (ms), por categoria y bus [min, mediana, max]:\n');
for c = unique(T.cat)'
    for b = buses
        m = T.cat == c & T.bus == b{1}; x = 1e3*(T.t_inc_est(m) - T.t_on(m));
        P('  %-6s %-4s n = %3d  [%7.2f %7.2f %7.2f]  sin estimacion: %d\n', c, b{1}, nnz(m), min(x), median(x, 'omitnan'), max(x), nnz(isnan(x)));
    end
end
P('Prefalta: |V0|/|V1| maximo %.2e, |V2|/|V1| mediana %.4f y maximo %.4f; pico de secuencia directa de %.1f a %.1f V.\n', ...
    max(T.V0pre_rel), median(T.V2pre_rel), max(T.V2pre_rel), min(T.V1pre_V), max(T.V1pre_V));

%% ---------------------------------------------------- 5. resumen por factor
T.pus = compose('%.4g', T.pu); T.pus = string(T.pus);
facs = {'todos', 'Srated', 'pus', 'ground', 'limit', 'faulted', 'variant'};
fid = fopen(fullfile(outdir, 'sintef_resumen_factor.csv'), 'w');
fprintf(fid, ['cat,bus,factor,nivel,n,correctas,erroneas,sin_declaracion,clases_declaradas,decl_ms_mediana,decl_ms_p90,decl_ms_max,' ...
    'primera_etiqueta_ms_mediana,primera_etiqueta_ms_max,sin_primera_etiqueta,registros_con_decl_prefalta,retencion_pct_mediana,retencion_pct_min\n']);
for c = unique(T.cat)'
    for b = buses
        m0 = T.cat == c & T.bus == b{1};
        for q = 1:numel(facs)
            if strcmp(facs{q}, 'todos'), lv = "todos"; else, lv = unique(T.(facs{q})(m0))'; end
            for l = lv
                if strcmp(facs{q}, 'todos'), m = m0; else, m = m0 & T.(facs{q}) == l; end
                x = T(m, :); ok = x.outcome == "correcta";
                fprintf(fid, '%s,%s,%s,%s,%d,%d,%d,%d,%s,%.2f,%.2f,%.2f,%.2f,%.2f,%d,%d,%.1f,%.1f\n', c, b{1}, strrep(facs{q}, 'pus', 'pu'), l, height(x), nnz(ok), ...
                    nnz(x.outcome == "erronea"), nnz(x.outcome == "sin_declaracion"), label_counts(cellstr(x.decl_class(x.declared))), ...
                    median(x.decl_ms(ok)), pct(x.decl_ms(ok), 90), max([x.decl_ms(ok); NaN]), median(x.first_label_ms, 'omitnan'), max(x.first_label_ms), ...
                    nnz(isnan(x.first_label_ms)), nnz(x.n_decl_pre > 0), median(x.full_correct_pct), min(x.full_correct_pct));
            end
        end
    end
end
fclose(fid);
% matriz clase declarada x verdad, por categoria, neutro y bus
fid = fopen(fullfile(outdir, 'sintef_matriz_clases.csv'), 'w');
fprintf(fid, 'cat,ground,bus,verdad,clase_declarada,n\n');
[g, c1, c2, c3, c4, c5] = findgroups(T.cat, T.ground, T.bus, T.truth, T.decl_class); cnt = splitapply(@numel, T.pu, g);
for i = 1:numel(cnt)
    dc = c5(i); if strlength(dc) == 0, dc = "sin_declaracion"; end
    fprintf(fid, '%s,%s,%s,%s,%s,%d\n', c1(i), c2(i), c3(i), c4(i), dc, cnt(i));
end
fclose(fid);
P('\nRESULTADO POR EVENTO (cadena completa, escala %s):\n', scales{1});
P('  %-6s %-4s %-5s | %4s %9s %8s %8s | %-24s | decl ms med/p90/max | 1a etiqueta ms med/max\n', 'cat', 'bus', 'neutro', 'n', 'correctas', 'erroneas', 'sin decl', 'clases declaradas');
for c = unique(T.cat)'
    for b = buses
        for gq = unique(T.ground)'
            x = T(T.cat == c & T.bus == b{1} & T.ground == gq, :); ok = x.outcome == "correcta";
            P('  %-6s %-4s %-5s | %4d %9d %8d %8d | %-24s | %6.2f %6.2f %6.2f | %6.2f %6.2f\n', c, b{1}, gq, height(x), nnz(ok), nnz(x.outcome == "erronea"), ...
                nnz(x.outcome == "sin_declaracion"), label_counts(cellstr(x.truth + ">" + x.decl_class)), median(x.decl_ms(ok)), pct(x.decl_ms(ok), 90), max([x.decl_ms(ok); NaN]), ...
                median(x.first_label_ms, 'omitnan'), max(x.first_label_ms));
        end
    end
end

%% ---------------------------------------------------------- 6. angulo L-L
LL = T(T.cat == "LL", :);
if ~isempty(LL)
    fid = fopen(fullfile(outdir, 'sintef_angulo_LL.csv'), 'w');
    fprintf(fid, ['bus,factor,nivel,n,desv_centro_mediana_deg,desv_p25_deg,desv_p75_deg,desv_min_deg,desv_max_deg,margen_frontera_min_regimen_deg,' ...
        'margen_frontera_min_falta_deg,V2V1_mediana,argV2V1_mediana_deg,redondez_mediana,pred_menos_clasificador_maxabs_deg,dispersion_intra_registro_iqr_max_deg,' ...
        'desv_fasores_min_deg,desv_fasores_max_deg,redondez_medida_min,redondez_medida_max,redondez_fasores_min,redondez_fasores_max,resid50_rel_min,resid50_rel_max\n']);
    facsL = {'todos', 'faulted', 'Srated', 'pus', 'ground', 'limit'};
    for b = [buses, {'ambos'}]
        if strcmp(b{1}, 'ambos'), m0 = true(height(LL), 1); else, m0 = LL.bus == b{1}; end
        for q = 1:numel(facsL)
            if strcmp(facsL{q}, 'todos'), lv = "todos"; else, lv = unique(LL.(facsL{q})(m0))'; end
            for l = lv
                if strcmp(facsL{q}, 'todos'), m = m0; else, m = m0 & LL.(facsL{q}) == l; end
                x = LL(m, :); dv = x.dev_centre_st_deg;
                fprintf(fid, '%s,%s,%s,%d,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.3f,%.2f,%.4f,%.2f,%.2f,%.2f,%.2f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n', b{1}, strrep(facsL{q}, 'pus', 'pu'), l, height(x), ...
                    median(dv, 'omitnan'), pct(dv, 25), pct(dv, 75), min(dv), max(dv), min(x.margin_min_st_deg), min(x.margin_min_full_deg), median(x.V2V1), ...
                    median(x.argV2V1_deg), median(x.roundness_st_med), max(abs(x.pred_minus_cls_deg)), max(x.theta_st_iqr), ...
                    min(x.dev_centre_pred_deg), max(x.dev_centre_pred_deg), min(x.roundness_st_med), max(x.roundness_st_med), min(x.roundness_pred), max(x.roundness_pred), min(x.resid50_rel), max(x.resid50_rel));
            end
        end
    end
    fclose(fid);
    % diferencias pareadas: mismo registro salvo el factor
    fid = fopen(fullfile(outdir, 'sintef_angulo_LL_pareado.csv'), 'w');
    fprintf(fid, 'factor,contraste,n_pares,dif_mediana_deg,dif_min_deg,dif_max_deg,dif_abs_mediana_deg,dif_abs_max_deg\n');
    allf = {'bus', 'faulted', 'Srated', 'pus', 'ground', 'limit'};
    con = {'bus', "DG2", "Pex"; 'Srated', "S245", "S490"; 'ground', "Pet0", "Pet1"; 'limit', "MB", "VI"; 'pus', "0.5", "1"; ...
           'faulted', "AB", "BC"; 'faulted', "AB", "AC"};
    P('\nANGULO L-L. Desviacion respecto al centro de sector (regimen de falta), diferencias pareadas por factor:\n');
    for q = 1:size(con, 1)
        dd = paired(LL, allf, con{q,1}, con{q,2}, con{q,3}, 'dev_centre_st_deg');
        fprintf(fid, '%s,%s menos %s,%d,%.2f,%.2f,%.2f,%.2f,%.2f\n', strrep(con{q,1}, 'pus', 'pu'), con{q,3}, con{q,2}, numel(dd), median(dd), min(dd), max(dd), median(abs(dd)), max(abs(dd)));
        P('  %-8s %-4s menos %-4s: n = %3d pares, mediana %+6.2f, [%+6.2f, %+6.2f] grados\n', strrep(con{q,1}, 'pus', 'pu'), con{q,3}, con{q,2}, numel(dd), median(dd), min(dd), max(dd));
    end
    % recorrido con la carga dentro de cada celda (resto de factores fijos)
    oth = setdiff(allf, {'pus'}, 'stable'); kk = join(LL{:, oth}, '|', 2); [gk, ~] = findgroups(kk);
    rg = splitapply(@(x) max(x) - min(x), LL.dev_centre_st_deg, gk); nn = splitapply(@numel, LL.dev_centre_st_deg, gk);
    fprintf(fid, 'pu,recorrido entre 0.5 y 1.0 por celda,%d,%.2f,%.2f,%.2f,%.2f,%.2f\n', numel(rg), median(rg), min(rg), max(rg), median(rg), max(rg));
    fclose(fid);
    P('  carga: recorrido de la desviacion entre los %d niveles dentro de cada celda: mediana %.2f, maximo %.2f grados (%d celdas)\n', mode(nn), median(rg), max(rg), numel(rg));
    P('  Global: desviacion mediana %+.2f grados, [%+.2f, %+.2f]; margen minimo a la frontera en regimen %.2f grados y en toda la falta %.2f grados\n', ...
        median(LL.dev_centre_st_deg, 'omitnan'), min(LL.dev_centre_st_deg), max(LL.dev_centre_st_deg), min(LL.margin_min_st_deg), min(LL.margin_min_full_deg));
    P('  Redondez medida por ventana (mediana en regimen) de %.3f a %.3f; redondez de la componente fundamental |1-k|/(1+k) de %.3f a %.3f; contenido no fundamental en falta de %.3f a %.3f pu (prefalta: %.4f a %.4f)\n', ...
        min(LL.roundness_st_med), max(LL.roundness_st_med), min(LL.roundness_pred), max(LL.roundness_pred), min(LL.resid50_rel), max(LL.resid50_rel), min(LL.resid50_pre_rel), max(LL.resid50_pre_rel));
    P('  Desviacion segun los fasores a 50 Hz: de %+.1f a %+.1f grados; segun el clasificador: de %+.1f a %+.1f grados\n', min(LL.dev_centre_pred_deg), max(LL.dev_centre_pred_deg), min(LL.dev_centre_st_deg), max(LL.dev_centre_st_deg));
    P('  |V2/V1| de %.3f a %.3f; arg(V2/V1) menos el nominal del par: de %+.1f a %+.1f grados; prediccion menos clasificador: max %.2f grados\n', ...
        min(LL.V2V1), max(LL.V2V1), -2*max(LL.dev_centre_pred_deg), -2*min(LL.dev_centre_pred_deg), max(abs(LL.pred_minus_cls_deg)));
end

%% ------------------------------------------------- 7. lectura de la tierra
GF = T(T.cat ~= "LL", :);
if ~isempty(GF)
    fid = fopen(fullfile(outdir, 'sintef_tierra.csv'), 'w');
    fprintf(fid, ['cat,ground,bus,n,V0_rel_mediana,V0_rel_min,V0_rel_max,argV0_resp_fase_mediana_deg,argV0_min_deg,argV0_max_deg,' ...
        'ventanas_tierra_pct_mediana,ventanas_tierra_pct_min,ventanas_fase_correcta_pct_mediana,ventanas_fase_correcta_pct_min,' ...
        'correctas,erroneas,sin_declaracion,clases_declaradas,decl_ms_mediana,decl_ms_max,primera_etiqueta_ms_mediana,primera_etiqueta_ms_max,' ...
        'tension_fase_falta_rel_mediana,tension_fases_sanas_rel_max,fase_mas_baja_es_la_de_falta,pinmax_falta_mediana,pinmax_falta_min\n']);
    P('\nLECTURA DE LA TIERRA (SLG y HISLG):\n');
    for c = unique(GF.cat)'
        for gq = unique(GF.ground)'
            for b = buses
                x = GF(GF.cat == c & GF.ground == gq & GF.bus == b{1}, :); ok = x.outcome == "correcta";
                am = angle(mean(exp(1j*deg2rad(x.argV0_rel_faulted_deg))))*180/pi; ad = wrap180(x.argV0_rel_faulted_deg - am);
                fprintf(fid, '%s,%s,%s,%d,%.3f,%.3f,%.3f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%d,%d,%d,%s,%.2f,%.2f,%.2f,%.2f,%.3f,%.3f,%d,%.3f,%.3f\n', c, gq, b{1}, height(x), ...
                    median(x.V0_rel), min(x.V0_rel), max(x.V0_rel), am, am + min(ad), am + max(ad), median(x.full_ground_pct), min(x.full_ground_pct), ...
                    median(x.full_correct_pct), min(x.full_correct_pct), nnz(ok), nnz(x.outcome == "erronea"), nnz(x.outcome == "sin_declaracion"), ...
                    label_counts(cellstr(x.truth + ">" + x.decl_class)), median(x.decl_ms(ok)), max([x.decl_ms(ok); NaN]), median(x.first_label_ms, 'omitnan'), max(x.first_label_ms), ...
                    median(x.vratio_faulted_min), max(x.vratio_others_max), nnz(x.lowest_phase == x.faulted), median(x.pinmax_full_med), min(x.pinmax_full_min));
                P('  %-6s %-5s %-4s n = %3d | V0/V1pre %.3f [%.3f, %.3f], arg %+6.1f [%+6.1f, %+6.1f] | ventanas con tierra %5.1f %% (min %5.1f), fase correcta %5.1f %% (min %5.1f) | %d correctas, %d erroneas, %d sin decl\n', ...
                    c, gq, b{1}, height(x), median(x.V0_rel), min(x.V0_rel), max(x.V0_rel), am, am + min(ad), am + max(ad), median(x.full_ground_pct), min(x.full_ground_pct), ...
                    median(x.full_correct_pct), min(x.full_correct_pct), nnz(ok), nnz(x.outcome == "erronea"), nnz(x.outcome == "sin_declaracion"));
            end
        end
    end
    fclose(fid);
end

%% ----------------------------------------------- 8. prefalta y tras despeje
fid = fopen(fullfile(outdir, 'sintef_prefalta.csv'), 'w');
fprintf(fid, ['bus,registros,ventanas_prefalta,ventanas_inicializacion,ventanas_armadas,etiqueta_final_no_sana_armadas_pct,etiqueta_cruda_no_sana_pct,etiqueta_final_no_sana_pct,etiqueta_cruda_tierra_pct,' ...
    'desde_0p2s_ventanas,desde_0p2s_cruda_no_sana_pct,desde_0p2s_final_no_sana_pct,desde_0p2s_cruda_tierra_pct,' ...
    'registros_con_declaracion_prefalta,declaraciones_prefalta,clases_prefalta,registros_con_declaracion_prefalta_sin_puerta,' ...
    'pinmax_prefalta_max,pinmax_desde_0p2s_max,V0_rel_prefalta_max,redondez_prefalta_mediana,redondez_prefalta_min,' ...
    'registros_con_declaracion_tras_despeje,clases_tras_despeje,ventanas_tras_despeje_no_sanas_pct\n']);
P('\nPREFALTA (1 s sano por registro) Y TRAS EL DESPEJE:\n');
for b = [buses, {'ambos'}]
    if strcmp(b{1}, 'ambos'), x = T; else, x = T(T.bus == b{1}, :); end
    wsum = @(p, w) sum(p.*w)/sum(w);
    fprintf(fid, '%s,%d,%d,%d,%d,%.3f,%.3f,%.3f,%.3f,%d,%.3f,%.3f,%.3f,%d,%d,%s,%d,%.3e,%.3e,%.3e,%.5f,%.5f,%d,%s,%.2f\n', b{1}, height(x), sum(x.pre_n_win), sum(x.pre_n_init), ...
        sum(x.pre_n_armed), wsum(x.pre_armed_nonhealthy_pct, x.pre_n_armed), ...
        wsum(x.pre_raw_nonhealthy_pct, x.pre_n_win), wsum(x.pre_nonhealthy_pct, x.pre_n_win), wsum(x.pre_raw_ground_pct, x.pre_n_win), ...
        sum(x.pre_late_n_win), wsum(x.pre_late_raw_nonhealthy_pct, x.pre_late_n_win), wsum(x.pre_late_nonhealthy_pct, x.pre_late_n_win), wsum(x.pre_late_raw_ground_pct, x.pre_late_n_win), ...
        nnz(x.n_decl_pre > 0), sum(x.n_decl_pre), label_counts(cellstr(x.pre_classes(x.n_decl_pre > 0))), nnz(x.n_decl_pre_nogate > 0), ...
        max(x.pinmax_pre_max), max(x.pinmax_pre_late_max), max(x.V0pre_rel), median(x.roundness_pre_late_med), min(x.roundness_pre_late_min), ...
        nnz(x.n_decl_after > 0), label_counts(cellstr(x.after_classes(x.n_decl_after > 0))), wsum(x.aft_nonhealthy_pct, x.aft_n_win));
    P('  %-5s %4d registros, %7d ventanas: etiqueta cruda no sana %.3f %%, final no sana %.3f %%, de tierra %.3f %% | desde 0.2 s: %.3f %%, %.3f %%, %.3f %% | declaraciones en prefalta: %d registros (%d sin puerta) | pinmax max %.2e | tras despeje: %d registros con declaracion\n', ...
        b{1}, height(x), sum(x.pre_n_win), wsum(x.pre_raw_nonhealthy_pct, x.pre_n_win), wsum(x.pre_nonhealthy_pct, x.pre_n_win), wsum(x.pre_raw_ground_pct, x.pre_n_win), ...
        wsum(x.pre_late_raw_nonhealthy_pct, x.pre_late_n_win), wsum(x.pre_late_nonhealthy_pct, x.pre_late_n_win), wsum(x.pre_late_raw_ground_pct, x.pre_late_n_win), ...
        nnz(x.n_decl_pre > 0), nnz(x.n_decl_pre_nogate > 0), max(x.pinmax_pre_max), nnz(x.n_decl_after > 0));
end
fclose(fid);

%% --------------------------------------------------------------- 9. escalado
fid = fopen(fullfile(outdir, 'sintef_escalado.csv'), 'w');
fprintf(fid, ['escala,cat,bus,n,correctas,erroneas,sin_declaracion,registros_con_declaracion_prefalta,registros_con_resultado_distinto,' ...
    'registros_con_latencia_distinta,ventanas_etiqueta_cruda_distinta,ventanas_etiqueta_final_distinta,ventanas_falta_etiqueta_final_distinta,' ...
    'ventanas_con_rama_linea_distinta,ventanas_totales\n']);
P('\nESCALADO (respecto a %s):\n', scales{1});
for is = 1:nsc
    sx = scales{is};
    for c = unique(T.cat)'
        for b = buses
            x = T(T.cat == c & T.bus == b{1}, :);
            if is == 1
                oc = x.outcome; dc = x.decl_class; dm = x.decl_ms; np = x.n_decl_pre; z = zeros(height(x), 1); nd = [z z z z];
            else
                oc = x.(['outcome_' sx]); dc = x.(['decl_class_' sx]); dm = x.(['decl_ms_' sx]); np = x.(['n_decl_pre_' sx]);
                nd = [x.(['ndiff_raw_' sx]), x.(['ndiff_final_' sx]), x.(['ndiff_final_fault_' sx]), x.(['ndiff_line_' sx])];
            end
            ddm = ~(dm == x.decl_ms | (isnan(dm) & isnan(x.decl_ms)));
            fprintf(fid, '%s,%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n', sx, c, b{1}, height(x), nnz(oc == "correcta"), nnz(oc == "erronea"), nnz(oc == "sin_declaracion"), ...
                nnz(np > 0), nnz(oc ~= x.outcome | dc ~= x.decl_class), nnz(ddm), sum(nd(:,1)), sum(nd(:,2)), sum(nd(:,3)), sum(nd(:,4)), sum(x.n_win));
            P('  %-7s %-6s %-4s: %3d correctas, %3d erroneas, %3d sin decl | resultado distinto en %d registros, latencia distinta en %d | ventanas distintas: cruda %d, final %d (en falta %d) de %d\n', ...
                sx, c, b{1}, nnz(oc == "correcta"), nnz(oc == "erronea"), nnz(oc == "sin_declaracion"), nnz(oc ~= x.outcome | dc ~= x.decl_class), nnz(ddm), sum(nd(:,1)), sum(nd(:,2)), sum(nd(:,3)), sum(x.n_win));
        end
    end
end
fclose(fid);
md5_fin = r2_t20_md5();
if strcmp(md5_ini, md5_fin), P('\nVersion de la cadena al terminar: la misma que al empezar.\n');
else, P('\nAVISO: la cadena cambio en disco durante la ejecucion. MD5 al terminar: %s\n', md5_fin); end
fclose(logf);
end

%% ======================================================================
function L = parse_label(name, relpath)
[~, nm] = fileparts(name); tok = strsplit(nm, '_');
L.cat = tok{1}; L.faulted = tok{2}; L.Srated = tok{3}; L.pu = str2double(erase(tok{4}, 'pu')); L.ground = tok{5}; L.limit = tok{6};
L.variant = 'normal'; if contains(lower(relpath), 'ealyed'), L.variant = 'retrasada'; end
switch L.faulted
    case {'A','B','C'}, L.truth = [L.faulted '-G'];
    case 'AB', L.truth = 'A-B';
    case 'BC', L.truth = 'B-C';
    case {'AC','CA'}, L.truth = 'C-A';
    otherwise, error('fases no reconocidas: %s', name);
end
p = strsplit(relpath, filesep); m = 0;
if numel(p) >= 5
    m = strcmp(p{1}, L.cat) + strcmp(p{2}, L.faulted) + strcmp(p{3}, L.Srated) + strcmp(p{5}, L.ground);
    s = erase(p{4}, 'pu'); if startsWith(s, '0') && ~contains(s, '.'), s = ['0.' s(2:end)]; end
    m = m + (abs(str2double(s) - L.pu) < 1e-9);
end
L.folder_match = m;
end

function [S, ph] = seqph(v, t, f, A)
X = [cos(2*pi*f*t), sin(2*pi*f*t), ones(numel(t), 1)];
c = X \ v; ph = (c(1,:) - 1j*c(2,:)).'; S = A*ph;
end

function r = add_event(r, E, suf)
fn = {'outcome','declared','decl_class','decl_ms','correct_decl_ms','first_label_ms','first_label_raw_ms','n_decl_correct','n_decl_wrong', ...
      'wrong_classes','max_run_correct','retention_raw_pct'};
for k = 1:numel(fn), r.([fn{k} suf]) = E.(fn{k}); end
end

function s = label_counts(lab)
if isempty(lab), s = ''; return, end
lab = cellstr(lab); lab(cellfun(@isempty, lab)) = {'-'};
[u, ~, j] = unique(lab); c = accumarray(j(:), 1); [c, o] = sort(c, 'descend'); u = u(o);
s = strjoin(arrayfun(@(i) sprintf('%s:%d', u{i}, c(i)), 1:numel(u), 'UniformOutput', false), ' ');
end

function [m, iq, mn, mx, n] = axial_stats(th)
% Estadistica de angulos axiales (periodo 180): media circular, IQR, minimo y maximo
th = th(isfinite(th)); n = numel(th);
if n == 0, m = NaN; iq = NaN; mn = NaN; mx = NaN; return, end
m = mod(angle(mean(exp(2j*deg2rad(th))))/2*180/pi, 180);
d = wrap90(th - m); m = mod(m + median(d), 180); d = wrap90(th - m);
iq = pct(d, 75) - pct(d, 25); mn = m + min(d); mx = m + max(d);
end

function d = paired(X, allf, fac, l1, l2, var)
% Diferencia de var entre los niveles l2 y l1 de fac, con el resto de factores iguales
oth = setdiff(allf, {fac}, 'stable');
a = X(X.(fac) == l1, :); b = X(X.(fac) == l2, :);
ka = join(a{:, oth}, '|', 2); kb = join(b{:, oth}, '|', 2);
[tf, loc] = ismember(ka, kb); d = b.(var)(loc(tf)) - a.(var)(tf); d = d(isfinite(d));
end

function y = wrap90(x), y = mod(x + 90, 180) - 90; end
function y = wrap180(x), y = mod(x + 180, 360) - 180; end
function p = pct(x, q)
x = sort(x(isfinite(x))); if isempty(x), p = NaN; return, end
if isscalar(x), p = x; return, end
p = interp1(linspace(0.5/numel(x), 1 - 0.5/numel(x), numel(x)), x, q/100, 'linear', 'extrap');
p = min(max(p, x(1)), x(end));
end
function plog(fid, varargin), fprintf(varargin{:}); fprintf(fid, varargin{:}); end
