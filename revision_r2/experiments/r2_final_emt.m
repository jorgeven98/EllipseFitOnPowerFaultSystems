function r2_final_emt(varargin)
% R2_FINAL_EMT  Postproceso DEFINITIVO de la campana EMT de T23 con el clasificador
% final de la ronda 2 (T90, bloque EMT). Peticiones R1.2 (S1) y R4.2, E2 (S2).
%
%   r2_final_emt                      % lo regenera todo
%   r2_final_emt('partes', {'resumen','tex'})   % solo resumenes y textos (sin reclasificar)
%   r2_final_emt('ids', {'S2a_AG_001'}, 'out_dir', carpeta)   % prueba de humo
%
% Clasificador: el publicado con la puerta de secuencia cero G5 de umbral 0.05
% (decision D2). Las opciones salen SIEMPRE de r2_final_opts y la decision de
% r2_decision_chain; aqui no se reimplementa ni la puerta ni la persistencia.
% No simula nada: lee las ondas de ronda2\emt\casos\ y el indice
% ronda2\emt\indice_casos.csv (solo lectura). Terminales: Pex (GFL, bus de la
% falta) y DG1 (GFM). Se excluyen los casos XC y, de los recuentos, las medidas
% no validas del indice (valida_Pex, valida_DG1).
% Misma ventana temporal y mismos t_on y t_off que r2_emt_postprocess (T23):
% t_on = caso.t_falta, t_off = t_on + caso.dur_falta, f = caso.red_f,
% t_range = [t_on - 1, min(t(end), t_off + 0.12)].
%
% Partes:
%   'base'     una decision por caso y terminal -> eventos_todos.csv, S1_final.csv,
%              S2_final.csv, excluidas.csv
%   'ruido'    barrido de Rf (S2a + S2d) con ruido de 40 dB y el peor caso
%              diferencial de clase 0.5 de T52 (retardo de 40/60 grados y +1 % de
%              ganancia en un mismo canal; canales A, B y C por turno), 5 semillas
%              fijas -> S2_ruido_eventos.csv
%   'tex'      r2_final_emt_tex.py (uv): textos\T23_tabla.tex y T23_figura.tex
%   'resumen'  r2_final_emt_resumen: controles, umbrales, recuentos, latencias,
%              cifras_emt.csv y seccion_informe_emt.md (va despues de 'tex'
%              porque recoge las cifras del pie de la tabla)
% Salidas en ronda2\resultados\T90\emt\ (o en 'out_dir'). Registro: r2_final_emt_log.txt
here = fileparts(mfilename('fullpath')); RAIZ = fileparts(here);
R2 = fullfile(RAIZ, '_revision_TPWRD', 'ronda2');
addpath(here, fullfile(RAIZ, 'src', 'fault_classifiers'), fullfile(RAIZ, 'src', 'utils'));
p.partes = {'base', 'ruido', 'resumen', 'tex'}; p.ids = {}; p.out_dir = fullfile(R2, 'resultados', 'T90', 'emt');
for k = 1:2:numel(varargin)
    if ~isfield(p, varargin{k}), error('r2_final_emt:opt', 'Opcion desconocida: %s', varargin{k}); end
    p.(varargin{k}) = varargin{k+1};
end
p.partes = cellstr(p.partes); p.ids = cellstr(p.ids); OUT = p.out_dir;
if ~exist(OUT, 'dir'), mkdir(OUT); end
logf = fullfile(OUT, 'r2_final_emt_log.txt');
if ismember('base', p.partes) && exist(logf, 'file'), delete(logf); end
diary(logf); cu = onCleanup(@() diary('off'));
ws = warning('off', 'all'); cw = onCleanup(@() warning(ws));
[o, info] = r2_final_opts();
fprintf('\n==== R2_FINAL_EMT  %s  MATLAB %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), version);
fprintf('Partes: %s\n', strjoin(p.partes, ', '));
fprintf('Opciones finales (r2_final_opts): %s\n', info.text);
fprintf('Opciones resueltas de la cadena:\n'); fn = fieldnames(info.resolved);
for k = 1:numel(fn)
    x = info.resolved.(fn{k});
    if isa(x, 'function_handle'), x = func2str(x); end
    if isnumeric(x) || islogical(x), x = mat2str(x); end
    fprintf('   %-14s = %s\n', fn{k}, char(x));
end
md5_ini = r2_t20_md5(); fprintf('Huella r2_t20_md5 al empezar: %s\n', md5_ini);
tic;

if ismember('base', p.partes), parte_base(R2, OUT, o, p.ids); end
if ismember('ruido', p.partes), parte_ruido(R2, OUT, o, p.ids); end
md5_fin = r2_t20_md5(); fprintf('Huella r2_t20_md5 tras la clasificacion: %s\n', md5_fin);
if ~strcmp(md5_ini, md5_fin)
    fprintf('AVISO: la huella ha cambiado durante la ejecucion. Repetir.\n');
else
    fprintf('La huella no ha cambiado durante la ejecucion.\n');
end
if ismember('tex', p.partes)
    cmd = sprintf('uv run --with pandas python -I "%s" "%s"', fullfile(here, 'r2_final_emt_tex.py'), OUT);
    [st, sal] = system(cmd); fprintf('%s\n%s\n', cmd, sal);
    if st ~= 0, fprintf('AVISO: r2_final_emt_tex.py ha fallado (estado %d). Ejecutarlo a mano.\n', st); end
end
if ismember('resumen', p.partes)
    r2_final_emt_resumen(OUT);
end
fprintf('Huella r2_t20_md5 al acabar: %s\n', r2_t20_md5());
fprintf('FIN r2_final_emt (%.0f s)\n', toc);
end

% =========================================================================
function parte_base(R2, OUT, o, ids)
% Una decision por caso y terminal con la cadena final.
I = leer_indice(R2, ids);
P21 = readtable(fullfile(R2, 'emt', 'T21_emt_cases_congelado.csv'), 'TextType', 'char');
TERM = {'Pex', 'DG1'}; ZBASE = 5.775;   % Vbase/Ibase = 231 V / 40 A (GridDesign.m)
rows = cell(1, 2*height(I)); nr = 0; n_incoh = 0;
for k = 1:height(I)
    S = load(fullfile(R2, 'emt', 'casos', [I.id{k} '.mat'])); c = S.caso;
    if ~strcmp(c.id, I.id{k}) || ~strcmp(c.tipo, I.tipo{k}) || abs(c.Rf - I.Rf_ohm(k)) > 1e-12 || abs(c.Rtierra - I.Rtierra_ohm(k)) > 1e-12 ...
            || abs(c.t_falta - I.t_falta_s(k)) > 2e-5 || abs(c.red_f - I.red_f_Hz(k)) > 1e-9
        n_incoh = n_incoh + 1; fprintf('INCOHERENCIA entre el archivo y el indice: %s\n', I.id{k});
    end
    tipo = c.tipo; f = c.red_f; t_on = c.t_falta; t_off = t_on + c.dur_falta;
    tp = upper(strrep(tipo, '-', '')); ph = find(ismember('ABC', tp));
    % angulo de incidencia medido en el bus de la falta (Pex), como en r2_emt_postprocess
    [Vp, tt] = extract_converter_data(S, 'Pex'); Vp = double(Vp); tt = double(tt(:));
    if numel(ph) == 1, vb = Vp(:, ph); else, vb = Vp(:, ph(1)) - Vp(:, ph(2)); end
    m = tt >= t_on - 5/f - 1e-9 & tt < t_on - 1e-9;
    Pb = 2*mean(vb(m).*exp(-1j*2*pi*f*(tt(m) - t_on)));
    ang_med = mod(rad2deg(angle(Pb)) + 90, 360);
    for q = 1:2
        [V, t] = extract_converter_data(S, TERM{q}); V = double(V); t = double(t(:));
        [V1p, V2p, V0p] = seq(V, t, t_on - 2/f, f);       % prefalta: ciclo que acaba un ciclo antes
        [V1f, V2f, V0f] = seq(V, t, t_on + 2/f, f);       % falta: tercer ciclo tras la incepcion
        rho = V2f/V1f; th_ley = mod(165 - rad2deg(angle(rho))/2, 180);
        out = r2_decision_chain(V, t, o, 'truth', tipo, 't_on', t_on, 't_off', t_off, 'f', f, ...
            't_range', [t_on - 1.0, min(t(end), t_off + 0.12)]);
        e = out.event; w = out.win;
        full = w.t_start >= t_on - 1e-9 & w.t_end <= t_off + 1e-9;
        th_emt = medcirc(w.angle_deg(full));
        r = struct();
        r.id = I.id{k}; r.conjunto = I.conjunto{k}; r.factor = I.factor{k}; r.terminal = TERM{q}; r.tipo = tipo;
        r.valida = double(I.(['valida_' TERM{q}])(k) == 1);
        r.Rf_ohm = c.Rf; r.Rf_pu = c.Rf/ZBASE; r.Rtierra_ohm = c.Rtierra; r.Lf_H = c.Lf; r.R_bucle_ohm = I.R_bucle_ohm(k);
        r.ang_inc_pedido_deg = I.ang_inc_pedido_deg(k); r.ang_inc_medido_deg = ang_med;
        r.red_Z_ohm = I.red_Z_ohm(k); r.red_XR = I.red_XR(k); r.red_f_Hz = f; r.Id_ref_pu = c.Id_ref_pu; r.dg1_dP_W = c.dg1_dP;
        r.dg1_on = double(c.dg1_on); r.dg2_on = double(c.dg2_on); r.t_on_s = t_on; r.t_off_s = t_off;
        r.V1_pre_Vpico = abs(V1p); r.V2V1_pre_pct = 100*abs(V2p)/abs(V1p); r.V0V1_pre_pct = 100*abs(V0p)/abs(V1p);
        r.V1_falta_pu = abs(V1f)/abs(V1p); r.V2V1_falta_pct = 100*abs(V2f)/abs(V1f); r.V0V1_falta_pct = 100*abs(V0f)/abs(V1f);
        r.argV2V1_falta_deg = rad2deg(angle(rho));
        r.theta_ley_deg = th_ley; r.theta_emt_deg = th_emt; r.err_ley_deg = abs(mod(th_emt - th_ley + 90, 180) - 90);
        r.T21_id = ''; r.T21_theta_pred_deg = NaN; r.T21_etiqueta_pred = ''; r.err_T21_deg = NaN;
        if strcmp(I.conjunto{k}, 'S1t')
            n21 = str2double(I.id{k}(end-2:end));
            r.T21_id = P21.id{n21}; r.T21_theta_pred_deg = P21.pred_theta_Pex_deg(n21); r.T21_etiqueta_pred = P21.pred_label_Pex{n21};
            r.err_T21_deg = abs(mod(th_emt - r.T21_theta_pred_deg + 90, 180) - 90);
        end
        r = campos_evento(r, e);
        r.full_modal_label = char(e.full_modal_label); r.full_ground_pct = e.full_ground_pct; r.full_healthy_pct = e.full_healthy_pct;
        r.retention_pct = e.retention_pct; r.max_run_correct = e.max_run_correct;
        r.full_v0stat_med = e.full_pinmax_med; r.full_v0stat_max = e.full_pinmax_max; r.full_roundness_med = e.full_roundness_med;
        r.pre_n_win = e.pre_n_win; r.n_full = e.n_full;
        nr = nr + 1; rows{nr} = r;
    end
    if mod(k, 25) == 0, fprintf('base: %d de %d casos (%.0f s)\n', k, height(I), toc); end
end
T = struct2table([rows{1:nr}], 'AsArray', true);
writetable(T, fullfile(OUT, 'eventos_todos.csv'));
V = T(T.valida == 1, :);
writetable(V(startsWith(V.conjunto, 'S1'), :), fullfile(OUT, 'S1_final.csv'));
writetable(V(startsWith(V.conjunto, 'S2'), :), fullfile(OUT, 'S2_final.csv'));
X = T(T.valida ~= 1, {'id', 'conjunto', 'factor', 'terminal', 'tipo', 'red_f_Hz', 'dg1_on', 'dg2_on', 'V2V1_pre_pct', 'V0V1_pre_pct'});
[~, loc] = ismember(X.id, I.id); X.aviso_indice = I.aviso(loc);
writetable(X, fullfile(OUT, 'excluidas.csv'));
fprintf('base: %d casos, %d eventos (caso y terminal), %d validos, %d excluidos; incoherencias archivo-indice: %d (%.0f s)\n', ...
    height(I), height(T), height(V), height(X), n_incoh, toc);
end

% =========================================================================
function parte_ruido(R2, OUT, o, ids)
% Barrido de Rf (S2a + S2d) con ruido de 40 dB y el peor caso diferencial de
% clase 0.5 de T52 (r2_v0gate_variants, parte 'combo': CLS 'cl0.5' = +1.0 % de
% ganancia y +40/60 grados de retardo en un mismo canal, los tres canales por
% turno). Tres condiciones: 'cl05' (solo el error de canal), 'snr40' (solo el
% ruido) y 'cl05_snr40' (las dos cosas; es la variante pedida).
I = leer_indice(R2, ids); I = I(ismember(I.conjunto, {'S2a', 'S2d'}), :);
[~, ordI] = ismember(I.id, leer_indice(R2, {}).id);   % posicion en el indice completo: fija la semilla
TERM = {'Pex', 'DG1'}; CH = 'ABC'; FASE = 40/60; GAN = 1.0; SNR = 40; NREP = 5;
C = struct('condicion', {}, 'canal', {}, 'rep', {}, 'fase', {}, 'gan', {}, 'snr', {});
for ch = 1:3, C(end+1) = struct('condicion', 'cl05', 'canal', ch, 'rep', 0, 'fase', FASE, 'gan', GAN, 'snr', Inf); end %#ok<AGROW>
for rep = 1:NREP, C(end+1) = struct('condicion', 'snr40', 'canal', 0, 'rep', rep, 'fase', 0, 'gan', 0, 'snr', SNR); end %#ok<AGROW>
for ch = 1:3
    for rep = 1:NREP, C(end+1) = struct('condicion', 'cl05_snr40', 'canal', ch, 'rep', rep, 'fase', FASE, 'gan', GAN, 'snr', SNR); end %#ok<AGROW>
end
fprintf('ruido: %d casos x 2 terminales x %d condiciones (clase 0.5: +%.4f grados de retardo y +%.1f %% de ganancia en un canal; SNR %d dB; %d semillas)\n', ...
    height(I), numel(C), FASE, GAN, SNR, NREP);
rows = cell(1, 2*height(I)*(numel(C) + 1)); nr = 0;
for k = 1:height(I)
    S = load(fullfile(R2, 'emt', 'casos', [I.id{k} '.mat'])); c = S.caso;
    tipo = c.tipo; f = c.red_f; t_on = c.t_falta; t_off = t_on + c.dur_falta;
    for q = 1:2
        [V, t] = extract_converter_data(S, TERM{q}); V = double(V); t = double(t(:));
        vpk = max(abs(V(t >= t_on - 0.6 & t < t_on - 0.2, :)), [], 'all');   % pico sano (T52: 2.4 a 2.8 s con t_on = 3 s)
        a = {'truth', tipo, 't_on', t_on, 't_off', t_off, 'f', f, 't_range', [t_on - 1.0, min(t(end), t_off + 0.12)]};
        ob = r2_decision_chain(V, t, o, a{:}); eb = ob.event;
        for i = 0:numel(C)
            r = struct(); r.id = I.id{k}; r.conjunto = I.conjunto{k}; r.terminal = TERM{q}; r.tipo = tipo;
            r.valida = double(I.(['valida_' TERM{q}])(k) == 1); r.Rf_ohm = c.Rf; r.Rf_pu = c.Rf/5.775;
            if i == 0
                r.condicion = 'base'; r.canal = ''; r.rep = 0; r.semilla = NaN; r.fase_deg = 0; r.ganancia_pct = 0; r.snr_dB = Inf; r.v_pico = vpk; e = eb;
            else
                x = C(i); pd = [0 0 0]; gg = [0 0 0]; canal = '';
                if x.canal > 0, pd(x.canal) = x.fase; gg(x.canal) = x.gan/100; canal = CH(x.canal); end
                if isfinite(x.snr)
                    seed = 5000000 + ordI(k)*1000 + q*100 + x.canal*10 + x.rep;
                    Ve = r2_apply_meas_error(V, t, 'phase_deg', pd, 'gain', gg, 'snr_db', x.snr, 'seed', seed, 'v_peak', vpk, 'f', f);
                else
                    seed = NaN;
                    Ve = r2_apply_meas_error(V, t, 'phase_deg', pd, 'gain', gg, 'v_peak', vpk, 'f', f);
                end
                oe = r2_decision_chain(Ve, t, o, a{:}); e = oe.event;
                r.condicion = x.condicion; r.canal = canal; r.rep = x.rep; r.semilla = seed; r.fase_deg = x.fase; r.ganancia_pct = x.gan; r.snr_dB = x.snr; r.v_pico = vpk;
            end
            r = campos_evento(r, e);
            r.max_run_correct = e.max_run_correct; r.retention_pct = e.retention_pct; r.full_ground_pct = e.full_ground_pct;
            r.pre_nonhealthy_pct = e.pre_nonhealthy_pct;
            r.base_outcome = char(eb.outcome); r.base_decl_class = char(eb.decl_class); r.base_decl_ms = eb.decl_ms;
            r.cambia = double(~strcmp(r.outcome, r.base_outcome) || ~strcmp(r.decl_class, r.base_decl_class));
            nr = nr + 1; rows{nr} = r;
        end
    end
    fprintf('ruido: %d de %d casos (%.0f s)\n', k, height(I), toc);
end
T = struct2table([rows{1:nr}], 'AsArray', true);
writetable(T, fullfile(OUT, 'S2_ruido_eventos.csv'));
fprintf('ruido: %d filas (%.0f s)\n', height(T), toc);
end

% =========================================================================
function r = campos_evento(r, e)
% Metricas por evento de r2_decision_chain que se guardan en todos los CSV.
r.outcome = char(e.outcome); r.decl_class = char(e.decl_class); r.decl_ms = e.decl_ms;
r.first_label_ms = e.first_label_ms; r.correct_decl_ms = e.correct_decl_ms;
r.n_decl_pre = e.n_decl_pre; r.pre_classes = char(e.pre_classes);
r.n_decl_correct = e.n_decl_correct; r.n_decl_wrong = e.n_decl_wrong; r.wrong_classes = char(e.wrong_classes);
r.n_decl_after = e.n_decl_after; r.after_classes = char(e.after_classes);
% resultado en cuatro clases: correcta, par_adyacente (falta entre fases declarada como
% otro par entre fases), erronea_otra y sin_declaracion
LL = {'A-B', 'B-C', 'C-A'};
if strcmp(r.outcome, 'erronea')
    if ismember(char(e.truth), LL) && ismember(r.decl_class, LL), r.resultado4 = 'par_adyacente'; else, r.resultado4 = 'erronea_otra'; end
else
    r.resultado4 = r.outcome;
end
end

function I = leer_indice(R2, ids)
opt = detectImportOptions(fullfile(R2, 'emt', 'indice_casos.csv'), 'TextType', 'char');
opt = setvartype(opt, {'id', 'conjunto', 'factor', 'estado', 'tipo', 'aviso'}, 'char');
I = readtable(fullfile(R2, 'emt', 'indice_casos.csv'), opt);
I = I(strcmp(I.estado, 'ok') & ~strcmp(I.conjunto, 'XC'), :);
if ~isempty(ids), I = I(ismember(I.id, ids), :); end
end

function [V1, V2, V0] = seq(V, t, t0, f)
% Fasores de un ciclo [t0, t0 + 1/f) por minimos cuadrados y secuencias (ref. fase a).
% Copia de la funcion local de r2_emt_postprocess.m (T23).
m = t >= t0 - 1e-9 & t < t0 + 1/f - 1e-9;
A = [cos(2*pi*f*t(m)), -sin(2*pi*f*t(m)), ones(nnz(m), 1)];
X = A \ V(m, :);
P = X(1, :) + 1j*X(2, :);
a = exp(2j*pi/3);
V0 = sum(P)/3; V1 = (P(1) + a*P(2) + a^2*P(3))/3; V2 = (P(1) + a^2*P(2) + a*P(3))/3;
end

function m = medcirc(th)
% Mediana de angulos con periodo 180 grados (copia de r2_emt_postprocess.m).
th = th(isfinite(th));
if isempty(th), m = NaN; return, end
c = angle(mean(exp(2j*deg2rad(th))))/2;
d = mod(deg2rad(th) - c + pi/2, pi) - pi/2;
m = mod(rad2deg(c + median(d)), 180);
end
