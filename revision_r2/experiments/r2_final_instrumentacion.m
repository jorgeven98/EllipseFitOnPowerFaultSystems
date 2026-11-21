function r2_final_instrumentacion(out_dir)
% R2_FINAL_INSTRUMENTACION  Errores de instrumentacion y ruido, por evento, sobre las
% 24 faltas del paper (12 archivos x Pex y DG1) con el clasificador FINAL (T90).
%
%   r2_final_instrumentacion          % salidas en ronda2\resultados\T90\instrumentacion
%
% Rejilla (subconjunto de la de T52, con sus mismas condiciones, signos y semillas,
% para poder comparar evento a evento con la columna G5_005 de T52):
%   fase en un canal (A, B y C por turno): -1, -0.5, -0.25, +0.25, +0.5, +1, +2 grados
%   ganancia en un canal: -3, -1, +1, +3 %
%   peor caso diferencial de clase 0.5 (retardo de 40/60 grados y +1 % en un canal) y
%   de clase 1.0 (80/60 grados y +2 %), sin ruido
%   ruido blanco de 60, 50, 40 y 30 dB, 20 semillas (las de T02)
%   clase 0.5 con ruido de 40 dB, 5 semillas
% El error se aplica a la senal completa (r2_apply_meas_error) antes de segmentar y la
% decision es la cadena completa con r2_final_opts (puerta adaptativa causal, M = 10,
% marcas al final de la ventana), registro t = 2.0 a 3.2 s (1 s de prefalta, de la que
% las 80 primeras ventanas inicializan el nivel ambiente).
% Falsa en sano = declaracion asegurada con t_decl <= t_on (n_decl_pre > 0).
% Salidas: instr_eventos.csv (una fila por evento), instr_resumen.csv (por condicion,
% terminal y familia), instr_vs_T52.csv (comparacion evento a evento), instr_log.txt
here = fileparts(mfilename('fullpath')); RAIZ = fileparts(here);
addpath(here, fullfile(RAIZ, 'src', 'fault_classifiers'), fullfile(RAIZ, 'src', 'utils'));
R2 = fullfile(RAIZ, '_revision_TPWRD', 'ronda2');
if nargin < 1 || isempty(out_dir), out_dir = fullfile(R2, 'resultados', 'T90', 'instrumentacion'); end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
logf = fullfile(out_dir, 'instr_log.txt'); if exist(logf, 'file'), delete(logf); end
diary(logf); cu = onCleanup(@() diary('off'));
ws = warning('off', 'all'); cw = onCleanup(@() warning(ws));
[FO, info0] = r2_final_opts();
fprintf('R2_FINAL_INSTRUMENTACION  %s  MATLAB %s\nOpciones finales: %s\nHuella: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), version, info0.text, info0.md5);
[T, info] = r2_sim_cases();
C = conds();
fprintf('%d condiciones x 24 eventos\n', numel(C));
rows = {}; tic;
for k = 1:height(T)
    d = load(fullfile(info.dir, T.file{k}));
    for c = 1:2
        cv = info.converters{c}; [V, t] = extract_converter_data(d, cv); V = double(V); t = double(t(:));
        vpk = max(abs(V(t >= 2.4 & t < 2.8, :)), [], 'all');
        for i = 1:numel(C)
            q = C(i);
            if isfinite(q.snr)
                if strcmp(q.parte, 'ruido'), seed = ((k - 1)*2 + c)*100000 + q.snr*100 + q.rep;      % semilla de T02 y T52
                else, seed = 7000000 + ((k - 1)*2 + c)*10000 + q.i52; end                         % semilla de T52 (parte combo)
                Ve = r2_apply_meas_error(V, t, 'phase_deg', q.phase_deg, 'gain', q.gain, 'snr_db', q.snr, 'seed', seed, 'v_peak', vpk);
            else
                seed = NaN; Ve = r2_apply_meas_error(V, t, 'phase_deg', q.phase_deg, 'gain', q.gain, 'v_peak', vpk);
            end
            o = r2_decision_chain(Ve, t, FO, 'truth', T.truth{k}, 't_on', info.t_on, 't_off', info.t_off, 't_range', info.t_range);
            e = o.event;
            fam = T.family{k}; if startsWith(fam, 'LL') && ~strcmp(fam, 'LLG'), famg = 'L-L'; elseif strcmp(fam, 'LG'), famg = 'L-G'; else, famg = 'L-L-G'; end
            rows(end+1, :) = {q.parte, q.id, q.grupo, q.nivel, q.snr, q.rep, seed, T.file{k}, T.truth{k}, famg, cv, ...
                e.outcome, e.decl_class, e.decl_ms, e.correct_decl_ms, e.first_label_ms, e.max_run_correct, e.retention_pct, ...
                e.n_decl_wrong, e.wrong_classes, e.n_decl_pre, e.pre_classes}; %#ok<AGROW>
        end
    end
    fprintf('archivo %d de %d (%.0f s)\n', k, height(T), toc);
end
TE = cell2table(rows, 'VariableNames', {'parte','condicion','grupo','nivel','SNR_dB','rep','semilla','archivo','truth','familia','terminal', ...
    'resultado','clase','decl_ms','correct_decl_ms','primera_etq_ms','racha_max','retencion_pct','n_decl_erroneas','clases_erroneas','n_decl_prefalta','clases_prefalta'});
writetable(TE, fullfile(out_dir, 'instr_eventos.csv'));

% ---- resumen por grupo y nivel, terminal y familia
TE.gk = strcat(TE.grupo, '|', string(TE.nivel));
G = unique(TE(:, {'parte','grupo','nivel'}), 'stable'); res = {};
for g = 1:height(G)
    sg = strcmp(TE.parte, G.parte{g}) & strcmp(TE.grupo, G.grupo{g}) & TE.nivel == G.nivel(g);
    for term = {'Pex', 'DG1', 'los dos'}
        st = sg; if ~strcmp(term{1}, 'los dos'), st = st & strcmp(TE.terminal, term{1}); end
        for fm = {'L-L', 'L-G', 'L-L-G', 'todas'}
            s = st; if ~strcmp(fm{1}, 'todas'), s = s & strcmp(TE.familia, fm{1}); end
            if ~any(s), continue, end
            ok = s & strcmp(TE.resultado, 'correcta'); lat = TE.decl_ms(ok);
            res(end+1, :) = {G.parte{g}, G.grupo{g}, G.nivel(g), term{1}, fm{1}, nnz(s), nnz(ok), nnz(s & strcmp(TE.resultado, 'erronea')), ...
                nnz(s & strcmp(TE.resultado, 'sin_declaracion')), sum(TE.n_decl_erroneas(s)), nnz(s & TE.n_decl_prefalta > 0), ...
                med(lat), mx(lat)}; %#ok<AGROW>
        end
    end
end
TR = cell2table(res, 'VariableNames', {'parte','grupo','nivel','terminal','familia','n_eventos','correctos','primera_decl_erronea','sin_declaracion', ...
    'n_decl_clase_erronea','registros_falsa_en_sano','latencia_mediana_ms','latencia_max_ms'});
writetable(TR, fullfile(out_dir, 'instr_resumen.csv'));
disp(TR(strcmp(TR.terminal, 'los dos') | strcmp(TR.familia, 'todas'), :));

% ---- comparacion con T52 (columna G5_005), evento a evento
T52 = fullfile(R2, 'resultados', 'T52'); cmp = {};
for p = {'fase', 'ganancia', 'ruido', 'combo'}
    f52 = fullfile(T52, sprintf('variants_eventos_%s.csv', p{1}));
    if ~isfile(f52), continue, end
    o = detectImportOptions(f52, 'TextType', 'char'); o = setvartype(o, intersect({'clase','clases_erroneas','clases_prefalta','clases_posfalta','etapa','resultado','condicion','grupo','variante','archivo','convertidor'}, o.VariableNames), 'char');
    A = readtable(f52, o); A = A(strcmp(A.variante, 'G5_005'), :);
    s = strcmp(TE.parte, p{1});
    for r = find(s)'
        m = find(strcmp(A.archivo, TE.archivo{r}) & strcmp(A.convertidor, TE.terminal{r}) & strcmp(A.condicion, TE.condicion{r}) & (A.rep == TE.rep(r)), 1);
        if isempty(m), continue, end
        cmp(end+1, :) = {p{1}, TE.condicion{r}, TE.rep(r), TE.archivo{r}, TE.terminal{r}, TE.resultado{r}, A.resultado{m}, TE.clase{r}, A.clase{m}, ...
            TE.decl_ms(r), A.decl_ms(m), TE.n_decl_prefalta(r), A.n_decl_prefalta(m), ...
            strcmp(TE.resultado{r}, A.resultado{m}) && strcmp(TE.clase{r}, A.clase{m}) && isequaln(round(TE.decl_ms(r), 6), round(A.decl_ms(m), 6)) && TE.n_decl_prefalta(r) == A.n_decl_prefalta(m)}; %#ok<AGROW>
    end
end
TC = cell2table(cmp, 'VariableNames', {'parte','condicion','rep','archivo','terminal','resultado_final','resultado_T52','clase_final','clase_T52','decl_ms_final','decl_ms_T52','prefalta_final','prefalta_T52','igual'});
writetable(TC, fullfile(out_dir, 'instr_vs_T52.csv'));
fprintf('Comparacion con T52 G5_005: %d eventos emparejados, %d identicos (resultado, clase, latencia y falsas en prefalta)\n', height(TC), nnz(TC.igual));
fprintf('Huella al terminar: %s\nFIN r2_final_instrumentacion (%.0f s)\n', r2_t20_md5(), toc);
end

function v = med(x), if isempty(x), v = NaN; else, v = median(x); end, end
function v = mx(x), if isempty(x), v = NaN; else, v = max(x); end, end

function C = conds()
% Mismas definiciones que conds() de r2_v0gate_variants.m (T52); i52 = indice de la
% condicion en la parte de T52 (fija la semilla de la parte combo).
C = struct('parte', {}, 'id', {}, 'phase_deg', {}, 'gain', {}, 'snr', {}, 'rep', {}, 'grupo', {}, 'nivel', {}, 'i52', {});
ch = 'ABC';
C(end+1) = mk('fase', 'sin_error', [0 0 0], [0 0 0], Inf, 0, 'base', 0, 1);
i = 1;
for c = 1:3
    for dg = [-1 -0.5 -0.25 0.25 0.5 1 2 4]
        i = i + 1; if dg == 4, continue, end
        pd = [0 0 0]; pd(c) = dg;
        C(end+1) = mk('fase', sprintf('fase_%c_%+g', ch(c), dg), pd, [0 0 0], Inf, 0, 'fase_1ch', abs(dg), i); %#ok<AGROW>
    end
end
i = 0;
for c = 1:3
    for g = [-3 -1 1 3 6]
        i = i + 1; if g == 6, continue, end
        gg = [0 0 0]; gg(c) = g/100;
        C(end+1) = mk('ganancia', sprintf('gan_%c_%+g', ch(c), g), [0 0 0], gg, Inf, 0, 'gan_1ch', abs(g), i); %#ok<AGROW>
    end
end
for s = [60 50 40 30]
    for rep = 1:20
        C(end+1) = mk('ruido', sprintf('snr_%d', s), [0 0 0], [0 0 0], s, rep, 'ruido', s, NaN); %#ok<AGROW>
    end
end
CLS = {'cl0.5', 1.0, 40/60; 'cl1.0', 2.0, 80/60; 'cl3P', 6.0, 4.0}; i = 0;
for q = 1:size(CLS, 1)
    for c = 1:3
        pd = [0 0 0]; pd(c) = CLS{q, 3}; gg = [0 0 0]; gg(c) = CLS{q, 2}/100;
        i = i + 1;
        if q <= 2, C(end+1) = mk('combo', sprintf('%s_%c', CLS{q, 1}, ch(c)), pd, gg, Inf, 0, CLS{q, 1}, Inf, i); end %#ok<AGROW>
        for s = [40 30]
            for rep = 1:5
                i = i + 1;
                if q == 1 && s == 40, C(end+1) = mk('combo', sprintf('%s_%c_snr%d', CLS{q, 1}, ch(c), s), pd, gg, s, rep, CLS{q, 1}, s, i); end %#ok<AGROW>
            end
        end
    end
end
end

function c = mk(parte, id, pd, gg, snr, rep, grupo, nivel, i52)
c = struct('parte', parte, 'id', id, 'phase_deg', pd, 'gain', gg, 'snr', snr, 'rep', rep, 'grupo', grupo, 'nivel', nivel, 'i52', i52);
end
