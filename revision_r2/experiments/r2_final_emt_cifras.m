function r2_final_emt_cifras(OUT)
% R2_FINAL_EMT_CIFRAS  Escribe cifras_emt.csv y seccion_informe_emt.md (T90, bloque EMT)
% a partir de los CSV de resumen de r2_final_emt_resumen. No calcula resultados:
% solo les da formato. Cada cifra lleva su CSV de origen y el valor que daba la
% puerta publicada (T23 y T23b), calculado por r2_final_emt_resumen con las mismas
% funciones sobre las columnas T03_* de T23 y la misma poblacion (Pex y DG1,
% medidas validas). La columna linea_main_rev2 sale de una busqueda en
% paper\main_rev2.tex hecha aqui (valor con su unidad en una linea que hable de
% resistencia de falta, incidencia o de la campana EMT).
here = fileparts(mfilename('fullpath')); RAIZ = fileparts(here); R2 = fullfile(RAIZ, '_revision_TPWRD', 'ronda2');
if nargin < 1 || isempty(OUT), OUT = fullfile(R2, 'resultados', 'T90', 'emt'); end
addpath(here);
L = @(f) r2_final_emt_leer(fullfile(OUT, f));
rel = strrep(strrep(OUT, [R2 filesep], ''), '\', '/');           % ruta relativa a ronda2
if strcmp(rel, strrep(OUT, '\', '/')), rel = 'resultados/T90/emt'; end
U = L('S2_umbrales.csv'); U0 = L('S2_umbrales_T23.csv'); RC = L('S2_recuentos.csv'); INC = L('S2_incidencia.csv'); TR = L('S2_tierra_resistiva.csv');
LAT = L('S2_latencias.csv'); SF = L('S1_por_factor.csv'); ST = L('S1_totales.csv'); ST0 = L('S1_totales_T23.csv'); SA = L('S1_angulo.csv');
C = L('controles.csv'); D = L('controles_diferencias.csv'); X = L('cambios_vs_T23.csv'); E = L('eventos_todos.csv'); W = L('S2_barrido_tabla.csv');
hayR = exist(fullfile(OUT, 'S2_ruido_resumen.csv'), 'file') == 2; NR = table(); NU = table(); NC = table();
if hayR, NR = L('S2_ruido_resumen.csv'); NU = L('S2_ruido_umbrales.csv'); NC = L('S2_ruido_cambios.csv'); end
main = fullfile(RAIZ, 'paper', 'main_rev2.tex'); ML = {};
if exist(main, 'file'), ML = cellstr(readlines(main, 'Encoding', 'UTF-8')); end
CTX = 'fault resistance|R_f|R_\{f\}|inception angle|EMT|electromagnetic transient|rf_sweep|5\.775';
ctx = find(~cellfun(@isempty, regexp(ML, CTX, 'once')))';
fprintf('main_rev2.tex: %d lineas; lineas con contexto de la campana EMT (%s): %s\n', numel(ML), CTX, mat2str(ctx));
SCR = 'experiments/r2_final_emt.m (r2_final_emt_resumen.m, r2_final_emt_cifras.m)';
M_G5 = 'puerta G5 0.05: con V0 por debajo del umbral la ventana va a la rama entre fases y la falta a tierra resistiva deja de declararse con clase de tierra';
NUEVA = 'cifra nueva de T90: no tiene equivalente en T23 ni en T23b';
Q = {};

    function add(id, desc, valor, unidad, csv, antiguo, motivo)
        if isempty(antiguo), cambia = 1; mot = NUEVA; if nargin > 6 && ~isempty(motivo), mot = [NUEVA '; ' motivo]; end
        elseif strcmp(valor, antiguo), cambia = 0; mot = 'se confirma con la puerta final';
        else, cambia = 1; mot = motivo; if isempty(mot), mot = M_G5; end
        end
        Q(end+1, :) = {id, desc, valor, unidad, [rel '/' csv], SCR, cambia, antiguo, grep_main_(ML, ctx, valor, unidad), mot}; %#ok<SETNU>
    end

%% poblacion
nv = @(t) nnz(E.valida == 1 & strcmp(E.terminal, t)); nx = @(t) nnz(E.valida ~= 1 & strcmp(E.terminal, t));
add('EMT_N_CASOS', 'Casos EMT procesados (sin los dos XC)', sprintf('%d', numel(unique(E.id))), 'casos', 'eventos_todos.csv', sprintf('%d', numel(unique(E.id))), '');
add('EMT_N_EVENTOS', 'Eventos (caso y terminal, Pex y DG1)', sprintf('%d', height(E)), 'eventos', 'eventos_todos.csv', sprintf('%d', height(E)), '');
add('EMT_N_EXCLUIDAS', 'Medidas no validas excluidas de los recuentos (Pex y DG1)', sprintf('%d', nnz(E.valida ~= 1)), 'eventos', 'excluidas.csv', sprintf('%d', nnz(E.valida ~= 1)), '');
add('EMT_N_EXCLUIDAS_PEX', 'Medidas no validas en Pex (prefalta no equilibrada, casos de frecuencia)', sprintf('%d', nx('Pex')), 'eventos', 'excluidas.csv', sprintf('%d', nx('Pex')), '');
add('EMT_N_EXCLUIDAS_DG1', 'Medidas no validas en DG1 (DG1 desconectada o prefalta no equilibrada)', sprintf('%d', nx('DG1')), 'eventos', 'excluidas.csv', sprintf('%d', nx('DG1')), '');
add('EMT_N_VALIDAS_PEX', 'Medidas validas en Pex (S1 y S2)', sprintf('%d', nv('Pex')), 'eventos', 'eventos_todos.csv', sprintf('%d', nv('Pex')), '');
add('EMT_N_VALIDAS_DG1', 'Medidas validas en DG1 (S1 y S2)', sprintf('%d', nv('DG1')), 'eventos', 'eventos_todos.csv', sprintf('%d', nv('DG1')), '');

%% S2: umbrales
for i = 1:height(U)
    tp = U.tipo{i}; tm = U.terminal{i}; tag = [strrep(tp, '-', '') '_' upper(tm)]; k = find(strcmp(U0.tipo, tp) & strcmp(U0.terminal, tm));
    pre = sprintf('Barrido de Rf (S2a+S2d), %s medida en %s: ', tp, tm);
    add(['EMT_S2_' tag '_OK_OHM'], [pre 'mayor Rf por fase con clasificacion correcta (todas las menores tambien correctas)'], fo(U.Rf_ok_contigua_ohm(i)), 'ohm', 'S2_umbrales.csv', fo(U0.Rf_ok_contigua_ohm(k)), '');
    add(['EMT_S2_' tag '_OK_PU'], [pre 'idem en pu de Zbase = 5.775 ohm'], fpu(U.Rf_ok_contigua_ohm(i)/5.775), 'pu', 'S2_umbrales.csv', fpu(U0.Rf_ok_contigua_ohm(k)/5.775), '');
    add(['EMT_S2_' tag '_FALLO_OHM'], [pre 'primera Rf del barrido sin clasificacion correcta'], fo(U.Rf_primer_fallo_ohm(i)), 'ohm', 'S2_umbrales.csv', fo(U0.Rf_primer_fallo_ohm(k)), '');
    add(['EMT_S2_' tag '_DECL_OHM'], [pre 'mayor Rf por fase con alguna declaracion, correcta o erronea'], fo(U.Rf_decl_max_ohm(i)), 'ohm', 'S2_umbrales.csv', fo(U0.Rf_decl_max_ohm(k)), '');
    add(['EMT_S2_' tag '_DECL_PU'], [pre 'idem en pu'], fpu(U.Rf_decl_max_ohm(i)/5.775), 'pu', 'S2_umbrales.csv', fpu(U0.Rf_decl_max_ohm(k)/5.775), '');
    add(['EMT_S2_' tag '_NODECL_OHM'], [pre 'menor Rf del barrido desde la que no se declara nada'], fo(U.Rf_nodecl_desde_ohm(i)), 'ohm', 'S2_umbrales.csv', fo(U0.Rf_nodecl_desde_ohm(k)), '');
    add(['EMT_S2_' tag '_NODECL_PU'], [pre 'idem en pu'], fpu(U.Rf_nodecl_desde_ohm(i)/5.775), 'pu', 'S2_umbrales.csv', fpu(U0.Rf_nodecl_desde_ohm(k)/5.775), '');
    add(['EMT_S2_' tag '_ERR'], [pre 'clase erronea declarada entre el umbral de clasificacion y el de declaracion'], vac(U.clases_erroneas{i}), 'texto', 'S2_umbrales.csv', vac(U0.clases_erroneas{k}), '');
    add(['EMT_S2_' tag '_MONOTONO'], [pre 'el resultado es monotono con Rf (1 si, 0 no)'], sprintf('%d', U.monotono_correcta(i) && U.monotono_declaracion(i)), 'logico', 'S2_umbrales.csv', sprintf('%d', U0.monotono_correcta(k) && U0.monotono_declaracion(k)), '');
    if any(tp == 'G'), col = 'V0V1'; nm = '|V0|/|V1|'; else, col = 'V2V1'; nm = '|V2|/|V1|'; end
    add(['EMT_S2_' tag '_FIRMA_OK'], [pre nm ' en la mayor Rf con clasificacion correcta'], fpc(U.([col '_en_ok_pct'])(i)), '%', 'S2_umbrales.csv', fpc(U0.([col '_en_ok_pct'])(k)), '');
    add(['EMT_S2_' tag '_FIRMA_DECL'], [pre nm ' en la mayor Rf con alguna declaracion'], fpc(U.([col '_en_decl_max_pct'])(i)), '%', 'S2_umbrales.csv', fpc(U0.([col '_en_decl_max_pct'])(k)), '');
    add(['EMT_S2_' tag '_FIRMA_NODECL'], [pre nm ' en la menor Rf desde la que no se declara nada'], fpc(U.([col '_en_nodecl_desde_pct'])(i)), '%', 'S2_umbrales.csv', fpc(U0.([col '_en_nodecl_desde_pct'])(k)), '');
end

%% S2: recuentos
GC = {'S2 completo', 'S2'; 'barrido (S2a+S2d)', 'BARRIDO'; 'incidencia (S2b)', 'S2B'; 'tierra resistiva (S2c)', 'S2C'};
for i = 1:height(RC)
    g = find(strcmp(GC(:, 1), RC.grupo{i})); tag = sprintf('%s_%s_%s', GC{g, 2}, upper(RC.terminal{i}), upper(strrep(RC.tipo{i}, '-', '')));
    add(['EMT_N_' tag], sprintf('%s, %s en %s: eventos con declaracion correcta, erronea y sin declaracion (de %d)', RC.grupo{i}, RC.tipo{i}, RC.terminal{i}, RC.n_eventos(i)), ...
        sprintf('%d/%d/%d', RC.correctas(i), RC.erroneas(i), RC.sin_declaracion(i)), 'correctas/erroneas/sin declaracion', 'S2_recuentos.csv', ...
        sprintf('%d/%d/%d', RC.T23_correctas(i), RC.T23_erroneas(i), RC.T23_sin_declaracion(i)), '');
end
q = strcmp(RC.grupo, 'S2 completo') & strcmp(RC.tipo, 'todos');
add('EMT_S2_PREFALTA', 'S2, Pex y DG1: eventos con alguna declaracion en prefalta', sprintf('%d de %d', sum(RC.eventos_con_decl_prefalta(q)), sum(RC.n_eventos(q))), 'eventos', 'S2_recuentos.csv', ...
    sprintf('%d de %d', sum(RC.T23_eventos_con_decl_prefalta(q)), sum(RC.n_eventos(q))), '');

%% S2b: incidencia
for i = 1:height(INC)
    tag = sprintf('%s_%s_R%s', upper(strrep(INC.tipo{i}, '-', '')), upper(INC.terminal{i}), strrep(sprintf('%g', INC.Rf_ohm(i)), '.', 'p'));
    add(['EMT_S2B_' tag], sprintf('Incidencia (S2b), %s en %s con Rf = %g ohm: resultado en los angulos de 0, 45, 90 y 135 grados y latencia de la declaracion correcta', INC.tipo{i}, INC.terminal{i}, INC.Rf_ohm(i)), ...
        inc_txt(INC.resultado_por_angulo{i}, INC.n_correcta(i), INC.decl_min_ms(i), INC.decl_max_ms(i)), 'texto (ms)', 'S2_incidencia.csv', ...
        inc_txt(INC.T23_resultado_por_angulo{i}, INC.T23_n_correcta(i), INC.T23_decl_min_ms(i), INC.T23_decl_max_ms(i)), '');
end
ig0 = 0; for i = 1:height(INC), ig0 = ig0 + (numel(unique(strsplit(strtrim(INC.T23_resultado_por_angulo{i}), ' '))) == 1); end
add('EMT_S2B_GRUPOS_IGUALES', 'Incidencia (S2b), Pex y DG1: grupos de tipo, terminal y Rf con el mismo resultado en los cuatro angulos', sprintf('%d de %d', nnz(INC.mismo_resultado_en_todos), height(INC)), 'grupos', 'S2_incidencia.csv', sprintf('%d de %d', ig0, height(INC)), '');
add('EMT_S2B_ANGULO_ERROR', 'Incidencia (S2b): diferencia maxima entre el angulo de incidencia medido y el pedido', sprintf('%.2f', max(INC.error_angulo_max_deg)), 'grados', 'S2_incidencia.csv', sprintf('%.2f', max(INC.error_angulo_max_deg)), '');

%% S2c: tierra resistiva
for i = 1:height(TR)
    tag = sprintf('%s_RG%g', upper(TR.terminal{i}), TR.Rtierra_ohm(i));
    add(['EMT_S2C_' tag], sprintf('A-B-G con fases francas y Rtierra = %g ohm (S2c), medida en %s (|V0|/|V1| = %.2f %%): resultado y latencia de la declaracion', TR.Rtierra_ohm(i), TR.terminal{i}, TR.V0V1_falta_pct(i)), ...
        ev_txt(TR.outcome{i}, TR.decl_class{i}, TR.decl_ms(i)), 'texto (ms)', 'S2_tierra_resistiva.csv', ev_txt(TR.T23_outcome{i}, TR.T23_decl_class{i}, TR.T23_decl_ms(i)), '');
end

%% S2: latencias
for i = 1:height(LAT)
    if ~ismember(LAT.grupo{i}, {'barrido (S2a+S2d)', 'incidencia (S2b)', 'S2 completo'}), continue, end
    g = find(strcmp(GC(:, 1), LAT.grupo{i})); tag = sprintf('%s_%s_%s', GC{g, 2}, upper(LAT.terminal{i}), upper(strrep(LAT.tipo{i}, '-', '')));
    pre = sprintf('%s, %s en %s, eventos con declaracion correcta (%d): ', LAT.grupo{i}, LAT.tipo{i}, LAT.terminal{i}, LAT.n_correctos(i));
    add(['EMT_LAT1_' tag], [pre 'primera etiqueta correcta, minimo a maximo'], rng_(LAT.primera_min_ms(i), LAT.primera_max_ms(i)), 'ms', 'S2_latencias.csv', rng_(LAT.T23_primera_min_ms(i), LAT.T23_primera_max_ms(i)), ...
        'cambia el conjunto de eventos correctos o las ventanas que la puerta manda a la rama de tierra');
    add(['EMT_LATD_' tag], [pre 'declaracion asegurada (M = 10), minimo a maximo y mediana'], rng3_(LAT.decl_min_ms(i), LAT.decl_max_ms(i), LAT.decl_mediana_ms(i)), 'ms', 'S2_latencias.csv', ...
        rng3_(LAT.T23_decl_min_ms(i), LAT.T23_decl_max_ms(i), LAT.T23_decl_mediana_ms(i)), 'cambia el conjunto de eventos correctos o las ventanas que la puerta manda a la rama de tierra');
end
q = strcmp(LAT.grupo, 'S2 completo') & strcmp(LAT.tipo, 'todos');
add('EMT_TAB_PRIMERA', 'Pie de T23_tabla.tex: primera etiqueta correcta en los eventos con declaracion correcta de S2 (Pex y DG1), minimo a maximo', ...
    sprintf('%.1f a %.1f', min(LAT.primera_min_ms(q)), max(LAT.primera_max_ms(q))), 'ms', 'S2_latencias.csv', sprintf('%.1f a %.1f', min(LAT.T23_primera_min_ms(q)), max(LAT.T23_primera_max_ms(q))), ...
    'la A-B-G con 100 ohm en el camino de tierra (primera etiqueta a 7.25 ms) ya no se declara como A-B-G');
if exist(fullfile(OUT, 'tabla_tex_cifras.csv'), 'file')
    TC = L('tabla_tex_cifras.csv'); v = @(n) TC.valor(strcmp(TC.cifra, n));
    if abs(v('primera_etiqueta_correcta_min_ms') - min(LAT.primera_min_ms(q))) > 1e-9 || abs(v('primera_etiqueta_correcta_max_ms') - max(LAT.primera_max_ms(q))) > 1e-9
        fprintf('AVISO: las cifras del pie de T23_tabla.tex (tabla_tex_cifras.csv) no coinciden con S2_latencias.csv. Regenerar el tex.\n');
    end
end

%% S2 con ruido de 40 dB y clase 0.5
if hayR
    CN = {'cl05_snr40', 'RUIDO', 'ruido de 40 dB mas peor caso diferencial de clase 0.5 (0.67 grados y 1 % en un canal; 3 canales x 5 semillas)'
          'snr40', 'SOLO40DB', 'solo ruido de 40 dB (5 semillas)'
          'cl05', 'SOLOCL05', 'solo el error de clase 0.5 (0.67 grados y 1 % en un canal; 3 canales)'};
    for c = 1:size(CN, 1)
        k = find(strcmp(NR.condicion, CN{c, 1}) & strcmp(NR.terminal, 'los dos') & strcmp(NR.tipo, 'todos'));
        pre = sprintf('Barrido de Rf (S2a+S2d) con %s, Pex y DG1: ', CN{c, 3});
        add(['EMT_' CN{c, 2} '_REALIZ_CAMBIAN'], [pre 'realizaciones (evento y perturbacion) cuyo resultado o clase cambia respecto al caso sin perturbar'], sprintf('%d de %d', NR.realizaciones_que_cambian(k), NR.n_realizaciones(k)), 'realizaciones', 'S2_ruido_resumen.csv', '', '');
        add(['EMT_' CN{c, 2} '_EVENTOS_CAMBIAN'], [pre 'eventos que cambian en alguna realizacion'], sprintf('%d de %d', NR.eventos_que_cambian_alguna_vez(k), NR.n_eventos(k)), 'eventos', 'S2_ruido_resumen.csv', '', '');
        add(['EMT_' CN{c, 2} '_CORRECTAS_PERDIDAS'], [pre 'realizaciones de eventos correctos sin perturbar que dejan de declararse bien'], sprintf('%d de %d', NR.correctas_en_base_que_dejan_de_serlo(k), NR.base_correctas(k)), 'realizaciones', 'S2_ruido_resumen.csv', '', '');
        add(['EMT_' CN{c, 2} '_PREFALTA'], [pre 'realizaciones con alguna declaracion en prefalta'], sprintf('%d de %d', NR.realizaciones_con_decl_prefalta(k), NR.n_realizaciones(k)), 'realizaciones', 'S2_ruido_resumen.csv', '', '');
        q = strcmp(NU.condicion, CN{c, 1});
        add(['EMT_' CN{c, 2} '_SERIES_UMBRAL_MOVIDO'], [pre 'series (tipo y terminal) en las que se mueve alguno de los tres umbrales en alguna realizacion'], ...
            sprintf('%d de %d%s', nnz(NU.se_mueve_algun_umbral(q)), nnz(q), lista_(strcat(NU.tipo(q & NU.se_mueve_algun_umbral == 1), {' en '}, NU.terminal(q & NU.se_mueve_algun_umbral == 1)))), 'series', 'S2_ruido_umbrales.csv', '', '');
    end
    for i = find(strcmp(NU.condicion, 'cl05_snr40'))'
        tp = NU.tipo{i}; tm = NU.terminal{i}; tag = [strrep(tp, '-', '') '_' upper(tm)];
        k = find(strcmp(NR.condicion, 'cl05_snr40') & strcmp(NR.terminal, tm) & strcmp(NR.tipo, tp));
        pre = sprintf('Barrido de Rf con ruido de 40 dB y clase 0.5, %s en %s: ', tp, tm);
        add(['EMT_RUIDO_' tag '_CAMBIAN'], [pre 'realizaciones que cambian de resultado o clase, y eventos afectados'], ...
            sprintf('%d de %d realizaciones; %d de %d eventos', NR.realizaciones_que_cambian(k), NR.n_realizaciones(k), NR.eventos_que_cambian_alguna_vez(k), NR.n_eventos(k)), 'texto', 'S2_ruido_resumen.csv', '', '');
        add(['EMT_RUIDO_' tag '_OK_OHM'], [pre 'mayor Rf con clasificacion correcta, minimo a maximo entre las 15 realizaciones (0 = ni la franca) y valor sin perturbar'], ...
            sprintf('%g a %g (sin perturbar %g; distinto en %d de %d)', NU.Rf_ok_contigua_min_ohm(i), NU.Rf_ok_contigua_max_ohm(i), NU.base_Rf_ok_contigua_ohm(i), NU.realizaciones_con_ok_distinto(i), NU.n_realizaciones(i)), 'ohm', 'S2_ruido_umbrales.csv', '', '');
        add(['EMT_RUIDO_' tag '_DECL_OHM'], [pre 'mayor Rf con alguna declaracion, minimo a maximo entre realizaciones y valor sin perturbar'], ...
            sprintf('%g a %g (sin perturbar %g; distinto en %d de %d)', NU.Rf_decl_max_min_ohm(i), NU.Rf_decl_max_max_ohm(i), NU.base_Rf_decl_max_ohm(i), NU.realizaciones_con_decl_distinto(i), NU.n_realizaciones(i)), 'ohm', 'S2_ruido_umbrales.csv', '', '');
        add(['EMT_RUIDO_' tag '_NODECL_OHM'], [pre 'menor Rf desde la que no se declara nada, minimo a maximo entre realizaciones y valor sin perturbar'], ...
            sprintf('%g a %g (sin perturbar %g; distinto en %d de %d)', NU.Rf_nodecl_desde_min_ohm(i), NU.Rf_nodecl_desde_max_ohm(i), NU.base_Rf_nodecl_desde_ohm(i), NU.realizaciones_con_nodecl_distinto(i), NU.n_realizaciones(i)), 'ohm', 'S2_ruido_umbrales.csv', '', '');
        add(['EMT_RUIDO_' tag '_DECL_MS'], [pre 'declaracion correcta, mediana y maximo, y mayor retraso respecto al caso sin perturbar'], ...
            sprintf('mediana %.2f, maximo %.2f, retraso maximo %.2f', NR.decl_correcta_mediana_ms(k), NR.decl_correcta_max_ms(k), NR.retraso_max_decl_correcta_ms(k)), 'ms', 'S2_ruido_resumen.csv', '', '');
    end
end

%% S1
for tm = {'Pex', 'DG1'}
    k = find(strcmp(SF.conjunto, 'S1 completo') & strcmp(SF.terminal, tm{1})); T = upper(tm{1});
    pre = sprintf('S1 (S1a+S1b+S1c+S1t), faltas entre fases medidas en %s: ', tm{1});
    add(['EMT_S1_' T '_VALIDAS'], [pre 'medidas validas'], sprintf('%d', SF.n_validas(k)), 'eventos', 'S1_por_factor.csv', sprintf('%d', SF.n_validas(k)), '');
    add(['EMT_S1_' T '_EXCLUIDAS'], [pre 'medidas no validas excluidas'], sprintf('%d', SF.n_excluidas(k)), 'eventos', 'S1_por_factor.csv', sprintf('%d', SF.n_excluidas(k)), '');
    add(['EMT_S1_' T '_CORRECTAS'], [pre 'declaracion correcta'], sprintf('%d', SF.correctas(k)), 'eventos', 'S1_por_factor.csv', sprintf('%d', SF.T23_correctas(k)), '');
    add(['EMT_S1_' T '_ADYACENTE'], [pre 'declaracion del par adyacente'], sprintf('%d', SF.par_adyacente(k)), 'eventos', 'S1_por_factor.csv', sprintf('%d', SF.T23_par_adyacente(k)), '');
    add(['EMT_S1_' T '_OTRA'], [pre 'declaracion de otra clase (de tierra)'], sprintf('%d', SF.erronea_otra(k)), 'eventos', 'S1_por_factor.csv', sprintf('%d', SF.T23_erronea_otra(k)), '');
    add(['EMT_S1_' T '_SINDECL'], [pre 'sin declaracion'], sprintf('%d', SF.sin_declaracion(k)), 'eventos', 'S1_por_factor.csv', sprintf('%d', SF.T23_sin_declaracion(k)), '');
    add(['EMT_S1_' T '_PREFALTA'], [pre 'eventos con alguna declaracion en prefalta'], sprintf('%d', SF.eventos_con_decl_prefalta(k)), 'eventos', 'S1_por_factor.csv', sprintf('%d', SF.eventos_con_decl_prefalta(k)), '');
    k1 = find(strcmp(ST.terminal, tm{1}) & strcmp(ST.tipo, 'todos')); k0 = find(strcmp(ST0.terminal, tm{1}) & strcmp(ST0.tipo, 'todos'));
    add(['EMT_S1_' T '_SENTIDO'], [pre 'sentido de las declaraciones erroneas'], vac(ST.sentido_del_error{k1}), 'texto', 'S1_totales.csv', vac(ST0.sentido_del_error{k0}), '');
end
for i = 1:height(SF)
    if strcmp(SF.conjunto{i}, 'S1 completo'), continue, end
    tag = regexprep(upper(sprintf('%s_%s_%s', SF.terminal{i}, SF.conjunto{i}, SF.factor{i})), '[^A-Z0-9]+', '_'); tag = regexprep(tag, '_$', '');
    add(['EMT_S1F_' tag], sprintf('S1, conjunto %s, factor %s, medida en %s: correctas, par adyacente y sin declaracion sobre las medidas validas (excluidas: %d)', SF.conjunto{i}, SF.factor{i}, SF.terminal{i}, SF.n_excluidas(i)), ...
        sprintf('%d/%d/%d de %d', SF.correctas(i), SF.par_adyacente(i), SF.sin_declaracion(i), SF.n_validas(i)), 'correctas/par adyacente/sin declaracion', 'S1_por_factor.csv', ...
        sprintf('%d/%d/%d de %d', SF.T23_correctas(i), SF.T23_par_adyacente(i), SF.T23_sin_declaracion(i), SF.n_validas(i)), '');
end
k = find(strcmp(SA.referencia, 'ley 165 - arg(V2/V1)/2') & strcmp(SA.terminal, 'Pex') & strcmp(SA.conjunto, 'S1 completo'));
add('EMT_S1_ANG_LEY_N', 'S1 en Pex: medidas validas en las que se compara el angulo de la elipse EMT con la ley 165 - arg(V2/V1)/2', sprintf('%d', SA.n(k)), 'eventos', 'S1_angulo.csv', sprintf('%d', SA.n(k)), '');
add('EMT_S1_ANG_LEY_MEDIANA', 'S1 en Pex: mediana de |theta_EMT - (165 - arg(V2/V1)/2)|', sprintf('%.1f', SA.error_mediana_deg(k)), 'grados', 'S1_angulo.csv', sprintf('%.1f', SA.T23_error_mediana_deg(k)), 'angulo calculado con la cadena final');
add('EMT_S1_ANG_LEY_MAX', sprintf('S1 en Pex: maximo de |theta_EMT - (165 - arg(V2/V1)/2)| (caso %s, V2/V1 = %.2f %%)', SA.id_del_maximo{k}, SA.V2V1_del_maximo_pct(k)), sprintf('%.1f', SA.error_max_deg(k)), 'grados', 'S1_angulo.csv', sprintf('%.1f', SA.T23_error_max_deg(k)), 'angulo calculado con la cadena final');
k = find(strcmp(SA.referencia, 'prediccion de T21') & startsWith(SA.poblacion, 'validas con'));
add('EMT_S1_ANG_T21_N', 'S1t en Pex: medidas validas con V2/V1 >= 2 % comparadas con la prediccion de T21', sprintf('%d', SA.n(k)), 'eventos', 'S1_angulo.csv', sprintf('%d', SA.n(k)), '');
add('EMT_S1_ANG_T21_MEDIANA', 'S1t en Pex: mediana de |theta_EMT - theta_T21| (validas con V2/V1 >= 2 %)', sprintf('%.1f', SA.error_mediana_deg(k)), 'grados', 'S1_angulo.csv', sprintf('%.1f', SA.T23_error_mediana_deg(k)), 'angulo calculado con la cadena final');
add('EMT_S1_ANG_T21_MAX', sprintf('S1t en Pex: maximo de |theta_EMT - theta_T21| (caso %s)', SA.id_del_maximo{k}), sprintf('%.1f', SA.error_max_deg(k)), 'grados', 'S1_angulo.csv', sprintf('%.1f', SA.T23_error_max_deg(k)), 'angulo calculado con la cadena final');

%% cambios por la puerta
gnd = contains(X.tipo, 'G');
ng = nnz(E.valida == 1 & contains(E.tipo, 'G')); nl = nnz(E.valida == 1 & ~contains(E.tipo, 'G'));
add('EMT_CAMBIO_TIERRA', 'Faltas a tierra de S2 (Pex y DG1): eventos cuyo resultado o clase cambia al pasar de la puerta publicada a la final', sprintf('%d de %d', nnz(gnd & X.cambia_resultado == 1), ng), 'eventos', 'cambios_vs_T23.csv', '', '');
add('EMT_CAMBIO_TIERRA_LAT', 'Faltas a tierra de S2 (Pex y DG1): eventos con el mismo resultado y otra latencia de declaracion', sprintf('%d de %d', nnz(gnd & X.cambia_resultado == 0), ng), 'eventos', 'cambios_vs_T23.csv', '', '');
add('EMT_CAMBIO_LL', 'Faltas entre fases de S1 y S2 (Pex y DG1): eventos cuyo resultado o clase cambia al pasar de la puerta publicada a la final', sprintf('%d de %d', nnz(~gnd & X.cambia_resultado == 1), nl), 'eventos', 'cambios_vs_T23.csv', '', '');
add('EMT_CAMBIO_LL_LAT', 'Faltas entre fases de S1 y S2 (Pex y DG1): eventos con el mismo resultado y otra latencia de declaracion', sprintf('%d de %d', nnz(~gnd & X.cambia_resultado == 0), nl), 'eventos', 'cambios_vs_T23.csv', '', '');

CF = cell2table(Q, 'VariableNames', {'id', 'descripcion', 'valor', 'unidad', 'csv_origen', 'script', 'cambia', 'valor_antiguo', 'linea_main_rev2', 'motivo'});
if numel(unique(CF.id)) ~= height(CF), error('r2_final_emt_cifras:id', 'Hay identificadores de cifra repetidos.'); end
writetable(CF, fullfile(OUT, 'cifras_emt.csv'), 'QuoteStrings', true);
fprintf('cifras_emt.csv: %d cifras; cambian %d (de ellas nuevas %d); se confirman %d; con linea en main_rev2.tex: %d\n', height(CF), nnz(CF.cambia), nnz(cellfun(@isempty, CF.valor_antiguo)), nnz(~CF.cambia), nnz(~cellfun(@isempty, CF.linea_main_rev2)));

if exist('r2_final_emt_informe', 'file') == 2   % T90: el redactor de la seccion no llego a escribirse; la seccion esta en informes/T90.md
    r2_final_emt_informe(OUT, CF, ctx);
else
    fprintf('r2_final_emt_informe no existe: seccion_informe_emt.md no se escribe (ver informes/T90.md, apartado 5).\n');
end
end

% =========================================================================
function s = grep_main_(ML, ctx, valor, unidad)
% Lineas de main_rev2.tex que hablan de la campana EMT (ctx) y contienen alguno de los
% numeros de la cifra con su unidad. Vacio si no hay ninguna.
s = '';
if isempty(ctx), return, end
num = regexp(valor, '(?<![\w.])\d+(\.\d+)?', 'match'); hit = [];
switch unidad
    case 'ohm', suf = '\s*(~|\\,|\\;)?\s*(\$\\Omega\$|\\Omega|ohm)';
    case 'pu', suf = '\s*(~|\\,|\\;)?\s*p\.?u';
    case {'ms', 'texto (ms)'}, suf = '\s*(~|\\,|\\;)?\s*ms';
    case 'grados', suf = '\s*(\^\{?\\circ|~?degrees)';
    case '%', suf = '\s*\\?%';
    otherwise, suf = '';
end
if isempty(suf) || isempty(num), return, end
for k = ctx
    for j = 1:numel(num)
        if ~isempty(regexp(ML{k}, ['(?<![\d.])' regexptranslate('escape', num{j}) '0*' suf], 'once')), hit(end+1) = k; break, end %#ok<AGROW>
    end
end
if ~isempty(hit), s = strjoin(arrayfun(@(x) sprintf('%d', x), hit, 'UniformOutput', false), ';'); end
end

function s = fo(x)
if isnan(x), s = 'ninguna'; else, s = sprintf('%g', x); end
end

function s = fpu(x)
if isnan(x), s = 'ninguna'; elseif x < 0.01, s = sprintf('%.4f', x); elseif x < 1, s = sprintf('%.3f', x); elseif x < 10, s = sprintf('%.2f', x); else, s = sprintf('%.1f', x); end
end

function s = fpc(x)
if isnan(x), s = 'no aplica'; else, s = sprintf('%.2f', x); end
end

function s = vac(x)
if isempty(x), s = 'ninguna'; else, s = x; end
end

function s = rng_(a, b)
if isnan(a), s = 'sin eventos correctos'; else, s = sprintf('%.2f a %.2f', a, b); end
end

function s = rng3_(a, b, m)
if isnan(a), s = 'sin eventos correctos'; else, s = sprintf('%.2f a %.2f (mediana %.2f)', a, b, m); end
end

function s = ev_txt(outcome, cls, ms)
switch outcome
    case 'correcta', s = sprintf('correcta (%s) a %.2f ms', cls, ms);
    case 'erronea', s = sprintf('erronea: declara %s a %.2f ms', cls, ms);
    otherwise, s = 'sin declaracion';
end
end

function s = inc_txt(res, nok, dmin, dmax)
tk = strsplit(strtrim(res), ' '); u = unique(tk, 'stable'); p = cell(1, numel(u));
for k = 1:numel(u), p{k} = sprintf('%s x%d', u{k}, nnz(strcmp(tk, u{k}))); end
s = strjoin(p, ', ');
if nok > 0, s = sprintf('%s; declaracion correcta de %.2f a %.2f ms', s, dmin, dmax); end
end

function s = lista_(c)
if isempty(c), s = ''; else, s = [' (' strjoin(c(:)', ', ') ')']; end
end
