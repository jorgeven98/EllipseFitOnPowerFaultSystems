function r2_final_nominal(out_dir)
% R2_FINAL_NOMINAL  Tabla nominal de precision (tab:acc_comparison, columna GA) y
% precision por tipo de falta (fig:acc_per_fault) con el clasificador FINAL (T90).
%
%   r2_final_nominal                 % salidas en ronda2\resultados\T90\nominal
%
% Protocolo de la tabla publicada (compare_classifier_accuracy.m, reproducido en
% T02 por r2_verify_codex_nominal): 12 archivos, DG1 y Pex, ventanas de 80, 60, 40
% y 20 muestras con paso round(0.25 W), ventanas que empiezan en t >= 3.00 s y
% acaban en t <= 3.16 s, etiqueta cruda del clasificador (gate = 'none'), acierto
% por archivo y media sobre las nueve clases (los tres pares L-L se promedian
% entre su archivo RL y su archivo franco antes de la media de clases).
% Se calcula con las opciones finales (r2_final_opts) y con la puerta publicada
% (valores por defecto) para dejar constancia de que no cambia ninguna celda.
% Ademas: una fila por evento (24 faltas) con la cadena completa final
% (puerta adaptativa causal, M = 10, t_range [2.0 3.2]).
% Salidas: nominal_por_archivo.csv, nominal_tabla.csv, nominal_por_tipo.csv,
% nominal_eventos.csv, nominal_log.txt
here = fileparts(mfilename('fullpath')); RAIZ = fileparts(here);
addpath(here, fullfile(RAIZ, 'src', 'fault_classifiers'), fullfile(RAIZ, 'src', 'utils'));
if nargin < 1 || isempty(out_dir), out_dir = fullfile(RAIZ, '_revision_TPWRD', 'ronda2', 'resultados', 'T90', 'nominal'); end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
logf = fullfile(out_dir, 'nominal_log.txt'); if exist(logf, 'file'), delete(logf); end
diary(logf); cu = onCleanup(@() diary('off'));
ws = warning('off', 'all'); cw = onCleanup(@() warning(ws));
[FO, info] = r2_final_opts();
fprintf('R2_FINAL_NOMINAL  %s  MATLAB %s\nOpciones finales: %s\nHuella: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), version, info.text, info.md5);
[T, si] = r2_sim_cases();
Ws = [80 60 40 20]; WL = {'1', '3/4', '1/2', '1/4'}; convs = {'DG1', 'Pex'};
VAR = {'final', FO; 'publicada', struct()};
rows = {}; acc = nan(height(T), numel(Ws), 2, 2); ev = {};
for k = 1:height(T)
    d = load(fullfile(si.dir, T.file{k}));
    for c = 1:2
        [V, t] = extract_converter_data(d, convs{c});
        for w = 1:numel(Ws)
            for v = 1:2
                o = r2_decision_chain(V, t, VAR{v, 2}, 'gate', 'none', 'W', Ws(w), 'stride', max(1, round(Ws(w)*0.25)), 't_range', [3.0 3.16]);
                acc(k, w, c, v) = 100*mean(strcmp(o.win.label_raw, T.truth{k}));
                rows(end+1, :) = {T.file{k}, T.truth{k}, convs{c}, Ws(w), VAR{v, 1}, height(o.win), acc(k, w, c, v)}; %#ok<AGROW>
            end
        end
        o = r2_decision_chain(V, t, FO, 'truth', T.truth{k}, 't_on', si.t_on, 't_off', si.t_off, 't_range', si.t_range);
        e = o.event;
        ev(end+1, :) = {T.file{k}, T.truth{k}, T.family{k}, convs{c}, e.outcome, e.decl_class, e.decl_ms, e.first_label_ms, e.n_decl_wrong, e.n_decl_pre, e.max_run_correct, e.retention_pct}; %#ok<AGROW>
    end
end
TA = cell2table(rows, 'VariableNames', {'archivo','truth','convertidor','W','variante','n_ventanas','acierto_pct'});
writetable(TA, fullfile(out_dir, 'nominal_por_archivo.csv'));
ut = unique(T.truth, 'stable'); tab = {}; tipo = {};
for v = 1:2
    for c = 1:2
        agg = nan(numel(ut), numel(Ws));
        for u = 1:numel(ut), agg(u, :) = mean(acc(strcmp(T.truth, ut{u}), :, c, v), 1); end
        m9 = mean(agg, 1);
        for w = 1:numel(Ws)
            tab(end+1, :) = {VAR{v, 1}, convs{c}, Ws(w), WL{w}, m9(w), round(m9(w), 1)}; %#ok<AGROW>
            for u = 1:numel(ut)
                tipo(end+1, :) = {VAR{v, 1}, convs{c}, Ws(w), WL{w}, ut{u}, agg(u, w), round(agg(u, w), 1)}; %#ok<AGROW>
            end
        end
    end
end
TT = cell2table(tab, 'VariableNames', {'variante','convertidor','W','ciclos','acierto_media9_pct','acierto_1dec'});
writetable(TT, fullfile(out_dir, 'nominal_tabla.csv'));
TP = cell2table(tipo, 'VariableNames', {'variante','convertidor','W','ciclos','clase','acierto_pct','acierto_1dec'});
writetable(TP, fullfile(out_dir, 'nominal_por_tipo.csv'));
TE = cell2table(ev, 'VariableNames', {'archivo','truth','familia','convertidor','outcome','decl_class','decl_ms','first_label_ms','n_decl_wrong','n_decl_pre','max_run_correct','retention_pct'});
writetable(TE, fullfile(out_dir, 'nominal_eventos.csv'));
fprintf('\nTabla nominal (media de 9 clases):\n'); disp(TT);
dif = max(abs(acc(:, :, :, 1) - acc(:, :, :, 2)), [], 'all');
fprintf('Diferencia maxima final - publicada por archivo: %.4g puntos\n', dif);
fprintf('Eventos (cadena final): %d de %d correctos; erroneas %d; en prefalta %d; latencia mediana %.2f ms, max %.2f ms\n', ...
    nnz(strcmp(TE.outcome, 'correcta')), height(TE), sum(TE.n_decl_wrong), sum(TE.n_decl_pre), median(TE.decl_ms), max(TE.decl_ms));
fprintf('Huella al terminar: %s\nFIN r2_final_nominal\n', r2_t20_md5());
end
