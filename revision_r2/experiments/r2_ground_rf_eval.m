function r2_ground_rf_eval(que)
% R2_GROUND_RF_EVAL  (T23b) Clasificador publicado frente a la regla corregida
% de la rama de tierra, por evento, con la cadena de referencia de T03
% (r2_decision_chain con sus valores por defecto; solo cambia 'classifier').
%   r2_ground_rf_eval            % todo
%   r2_ground_rf_eval('orig')    % 12 registros originales x (Pex, DG1)
%   r2_ground_rf_eval('emt')     % casos EMT de S2 (85) y XC (2) x (Pex, DG1, DG2)
%   r2_ground_rf_eval('nominal') % tabla nominal de precision por ventana
%   r2_ground_rf_eval('ll')      % indicacion de confianza en la rama entre fases
% La linea base y las variantes se ejecutan en la misma pasada, con el mismo
% archivo de cadena (su SHA-256 queda en el log).
% Variantes: base (ga_fault_classifier_pure); n15 (regla nueva con sus valores
% por defecto: |V2/V1| >= 0.01, |V2/V0| >= 0.15, estimador 'conj'); n10
% (|V2/V0| >= 0.10); n15b (|V2/V1| >= 0.02); n00 (sin la condicion |V2/V0|);
% g15 (estimador 'geom', sin frecuencia).
% Salidas en ronda2\resultados\T23b\: eval_orig.csv, eval_emt.csv,
% eval_nominal.csv, eval_nominal_archivos.csv, ll_confianza.csv,
% ll_confianza_resumen.csv y eval_log_<que>.txt
if nargin < 1, que = 'todo'; end
raiz  = fileparts(fileparts(mfilename('fullpath')));
raiz2 = fullfile(raiz,'_revision_TPWRD','ronda2');
res   = fullfile(raiz2,'resultados','T23b');
addpath(fullfile(raiz,'experiments'), fullfile(raiz,'src','fault_classifiers'), fullfile(raiz,'src','utils'));
logf = fullfile(res, ['eval_log_' que '.txt']); if exist(logf,'file'), delete(logf); end
diary(logf); cu = onCleanup(@() diary('off'));
warning('off','all');
fprintf('R2_GROUND_RF_EVAL(%s)  %s  MATLAB %s\n', que, char(datetime('now')), version);
[~, h] = system(['certutil -hashfile "' fullfile(raiz,'experiments','r2_decision_chain.m') '" SHA256']);
hh = regexp(h, '[0-9a-f]{64}', 'match', 'once'); fprintf('SHA-256 de r2_decision_chain.m: %s\n', hh);
[~, h] = system(['certutil -hashfile "' fullfile(raiz,'src','fault_classifiers','ga_fault_classifier_pure.m') '" SHA256']);
fprintf('SHA-256 de ga_fault_classifier_pure.m: %s\n', regexp(h, '[0-9a-f]{64}', 'match', 'once'));
VAR = {'base', @ga_fault_classifier_pure
       'n15',  @(V,fs,f) r2_ga_classifier_ground(V, fs, f, struct('u2_min', 0.01, 'v2v0_min', 0.15))
       'n10',  @(V,fs,f) r2_ga_classifier_ground(V, fs, f, struct('u2_min', 0.01, 'v2v0_min', 0.10))
       'n15b', @(V,fs,f) r2_ga_classifier_ground(V, fs, f, struct('u2_min', 0.02, 'v2v0_min', 0.15))
       'n00',  @(V,fs,f) r2_ga_classifier_ground(V, fs, f, struct('u2_min', 0.01, 'v2v0_min', 0))
       'g15',  @(V,fs,f) r2_ga_classifier_ground(V, fs, f, struct('u2_min', 0.01, 'v2v0_min', 0.15, 'metodo', 'geom'))};
todo = strcmpi(que, 'todo');

%% 1. Registros originales
if todo || strcmpi(que, 'orig')
    [T, info] = r2_sim_cases(); rows = {};
    for k = 1:height(T)
        d = load(fullfile(info.dir, T.file{k}));
        for c = 1:numel(info.converters)
            [V, t] = extract_converter_data(d, info.converters{c});
            lab0 = {};
            for v = 1:size(VAR,1)
                o = r2_decision_chain(V, t, 'truth', T.truth{k}, 't_on', info.t_on, 't_off', info.t_off, ...
                    't_range', info.t_range, 'classifier', VAR{v,2});
                if v == 1, lab0 = o.win.label_raw; end
                rows(end+1,:) = fila_(T.file{k}, 'orig', info.converters{c}, T.truth{k}, NaN, NaN, NaN, VAR{v,1}, o, lab0); %#ok<AGROW>
            end
        end
        fprintf('orig %s\n', T.file{k});
    end
    E = tabla_(rows); writetable(E, fullfile(res,'eval_orig.csv'));
    resumen_(E, 'Registros originales (12 x Pex, DG1)');
end

%% 2. Casos EMT con y sin tierra de S2, y XC
if todo || strcmpi(que, 'emt')
    I = readtable(fullfile(raiz2,'emt','indice_casos.csv'), 'TextType','char');
    I = I(ismember(I.conjunto, {'S2a','S2d','S2b','S2c','XC'}), :);
    convs = {'Pex','DG1','DG2'}; rows = {};
    for k = 1:height(I)
        S = load(fullfile(raiz2,'emt','casos',[I.id{k} '.mat'])); c = S.caso;
        t_on = c.t_falta; t_off = t_on + c.dur_falta;
        for q = 1:numel(convs)
            [V, t] = extract_converter_data(S, convs{q});
            lab0 = {};
            for v = 1:size(VAR,1)
                o = r2_decision_chain(V, t, 'truth', c.tipo, 't_on', t_on, 't_off', t_off, 'f', c.red_f, ...
                    't_range', [t_on - 1.0, min(t(end), t_off + 0.12)], 'classifier', VAR{v,2});
                if v == 1, lab0 = o.win.label_raw; end
                rows(end+1,:) = fila_(c.id, I.conjunto{k}, convs{q}, c.tipo, c.Rf, c.Rtierra, I.ang_inc_pedido_deg(k), VAR{v,1}, o, lab0); %#ok<AGROW>
            end
        end
        if mod(k, 10) == 0, fprintf('emt %d/%d\n', k, height(I)); end
    end
    E = tabla_(rows); writetable(E, fullfile(res,'eval_emt.csv'));
    resumen_(E(~strcmp(E.conjunto,'XC'),:), 'EMT S2 (85 casos x Pex, DG1, DG2)');
end

%% 3. Tabla nominal de precision por ventana (metodo de r2_verify_codex_nominal)
if todo || strcmpi(que, 'nominal')
    [T, info] = r2_sim_cases(); Ws = [80 60 40 20]; convs = {'DG1','Pex'}; rowsF = {}; rowsS = {};
    acc = nan(height(T), numel(Ws), 2, size(VAR,1));
    for k = 1:height(T)
        d = load(fullfile(info.dir, T.file{k}));
        for c = 1:2
            [V, t] = extract_converter_data(d, convs{c});
            for w = 1:numel(Ws)
                for v = 1:size(VAR,1)
                    o = r2_decision_chain(V, t, 'gate', 'none', 'W', Ws(w), 'stride', max(1, round(Ws(w)*0.25)), ...
                        't_range', [3.0 3.16], 'classifier', VAR{v,2});
                    acc(k,w,c,v) = 100*mean(strcmp(o.win.label_raw, T.truth{k}));
                    rowsF(end+1,:) = {T.file{k}, T.truth{k}, convs{c}, Ws(w), VAR{v,1}, height(o.win), acc(k,w,c,v)}; %#ok<AGROW>
                end
            end
        end
        fprintf('nominal %s\n', T.file{k});
    end
    writetable(cell2table(rowsF, 'VariableNames', {'archivo','truth','convertidor','W','variante','n_ventanas','acierto_pct'}), fullfile(res,'eval_nominal_archivos.csv'));
    ut = unique(T.truth, 'stable');
    for v = 1:size(VAR,1)
        agg = nan(numel(ut), numel(Ws), 2);
        for u = 1:numel(ut), agg(u,:,:) = mean(acc(strcmp(T.truth, ut{u}),:,:,v), 1); end
        m9 = squeeze(mean(agg, 1));
        for c = 1:2, for w = 1:numel(Ws)
            rowsS(end+1,:) = {convs{c}, Ws(w), VAR{v,1}, m9(w,c), min(agg(:,w,c))}; %#ok<AGROW>
        end, end
    end
    TS = cell2table(rowsS, 'VariableNames', {'convertidor','W','variante','acierto_medio_9_clases_pct','peor_clase_pct'});
    writetable(TS, fullfile(res,'eval_nominal.csv'));
    disp(unstack(TS(:,1:4), 'acierto_medio_9_clases_pct', 'variante'));
end

%% 4. Rama entre fases: par asegurado por geometria (redondez <= 1/3)
if todo || strcmpi(que, 'll')
    rows = {};
    [T, info] = r2_sim_cases();
    for k = find(ismember(T.family, {'LL','LL_RL'}))'
        d = load(fullfile(info.dir, T.file{k}));
        for c = 1:numel(info.converters)
            [V, t] = extract_converter_data(d, info.converters{c});
            o = r2_decision_chain(V, t, 'truth', T.truth{k}, 't_on', info.t_on, 't_off', info.t_off, 't_range', info.t_range);
            rows(end+1,:) = fila_ll_(T.file{k}, 'orig', info.converters{c}, T.truth{k}, NaN, 1, o); %#ok<AGROW>
        end
    end
    I = readtable(fullfile(raiz2,'emt','indice_casos.csv'), 'TextType','char');
    I = I(ismember(I.tipo, {'A-B','B-C','C-A'}) & ~strcmp(I.conjunto,'XC'), :);
    convs = {'Pex','DG1','DG2'};
    for k = 1:height(I)
        S = load(fullfile(raiz2,'emt','casos',[I.id{k} '.mat'])); c = S.caso;
        t_on = c.t_falta; t_off = t_on + c.dur_falta;
        for q = 1:numel(convs)
            [V, t] = extract_converter_data(S, convs{q});
            o = r2_decision_chain(V, t, 'truth', c.tipo, 't_on', t_on, 't_off', t_off, 'f', c.red_f, ...
                't_range', [t_on - 1.0, min(t(end), t_off + 0.12)]);
            rows(end+1,:) = fila_ll_(c.id, I.conjunto{k}, convs{q}, c.tipo, 2*c.Rf, I.(['valida_' convs{q}])(k), o); %#ok<AGROW>
        end
        if mod(k, 20) == 0, fprintf('ll %d/%d\n', k, height(I)); end
    end
    L = cell2table(rows, 'VariableNames', {'id','conjunto','conv','tipo','Rf_total_ohm','medida_valida','outcome','decl_class','decl_ms', ...
        'redondez_max_racha','redondez_med_falta','par_asegurado','n_vent_falta','pct_vent_aseguradas'});
    writetable(L, fullfile(res,'ll_confianza.csv'));
    L = L(L.medida_valida == 1, :);
    L.origen = repmat({'EMT'}, height(L), 1); L.origen(strcmp(L.conjunto,'orig')) = {'orig'};
    G = groupsummary(L, {'origen','conv','outcome','par_asegurado'});
    writetable(G, fullfile(res,'ll_confianza_resumen.csv')); disp(G);
end
end

%% ------------------------------------------------------------------------
function r = fila_(id, conj, conv, tipo, Rf, Rg, ang, var, o, lab0)
e = o.event; w = o.win;
enf = w.t_start >= e.t_on - 1e-9 & w.t_end <= e.t_off + 1e-9;
camb = NaN; if ~isempty(lab0), camb = 100*mean(~strcmp(w.label_raw(enf), lab0(enf))); end
r = {id, conj, conv, tipo, Rf, Rg, ang, var, e.outcome, e.decl_class, e.decl_ms, e.correct_decl_ms, e.first_label_ms, ...
    e.n_decl_correct, e.n_decl_wrong, e.wrong_classes, e.n_decl_pre, e.pre_classes, e.n_decl_after, e.max_run_correct, ...
    e.retention_pct, e.full_modal_label, e.full_ground_pct, camb};
end

function E = tabla_(rows)
E = cell2table(rows, 'VariableNames', {'id','conjunto','conv','tipo','Rf_ohm','Rg_ohm','ang_inc_deg','variante','outcome','decl_class','decl_ms', ...
    'correct_decl_ms','first_label_ms','n_decl_correct','n_decl_wrong','wrong_classes','n_decl_pre','pre_classes','n_decl_after', ...
    'max_run_correct','retention_pct','modal_label','ground_pct','pct_vent_cambiadas'});
end

function resumen_(E, titulo)
fprintf('\n== %s\n', titulo);
vs = unique(E.variante, 'stable'); cs = unique(E.conv, 'stable'); tp = unique(E.tipo, 'stable');
for c = 1:numel(cs)
    for i = 1:numel(tp)
        for v = 1:numel(vs)
            k = strcmp(E.conv, cs{c}) & strcmp(E.tipo, tp{i}) & strcmp(E.variante, vs{v});
            fprintf('%4s %6s %5s: %3d eventos, correcta %3d, erronea %3d, sin declaracion %3d, declaraciones en prefalta %d, latencia correcta mediana %.1f ms (max %.1f)\n', ...
                cs{c}, tp{i}, vs{v}, nnz(k), nnz(k & strcmp(E.outcome,'correcta')), nnz(k & strcmp(E.outcome,'erronea')), ...
                nnz(k & strcmp(E.outcome,'sin_declaracion')), sum(E.n_decl_pre(k)), ...
                median(E.decl_ms(k & strcmp(E.outcome,'correcta')), 'omitnan'), max([E.decl_ms(k & strcmp(E.outcome,'correcta')); NaN]));
        end
    end
end
end

function r = fila_ll_(id, conj, conv, tipo, Rft, valida, o)
% Par asegurado: las M ventanas de la racha que declara tienen redondez <= 1/3
% (|V2/V1| >= 1/2), condicion bajo la que el giro no cruza la frontera de
% sector con ningun X/R (T21, apartado 5).
e = o.event; w = o.win; M = o.opts.M;
enf = w.t_start >= e.t_on - 1e-9 & w.t_end <= e.t_off + 1e-9;
LL = ismember(w.label, {'A-B','B-C','C-A'});
rmax = NaN; aseg = 0;
if e.declared
    j = o.decl.j(find(o.decl.t_decl > e.t_on & o.decl.t_decl <= e.t_off + 1e-9, 1));
    rmax = max(w.roundness(j-M+1:j)); aseg = double(rmax <= 1/3);
end
r = {id, conj, conv, tipo, Rft, valida, e.outcome, e.decl_class, e.decl_ms, rmax, median(w.roundness(enf & LL), 'omitnan'), aseg, ...
    nnz(enf), 100*mean(w.roundness(enf & LL) <= 1/3)};
end
