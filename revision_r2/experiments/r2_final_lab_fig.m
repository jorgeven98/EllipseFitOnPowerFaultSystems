function R = r2_final_lab_fig(varargin)
% R2_FINAL_LAB_FIG  (T90, bloque de laboratorio) Datos de las figuras de laboratorio del
% paper con la senal adoptada en la ronda 2: conv1 referida al lado de red (Dy11).
%
%   r2_final_lab_fig()                      % escribe en ronda2\resultados\T90\figuras\
%   r2_final_lab_fig('W', 150, 'stamp', 'end', 'zthr', 0.9, 'outdir', carpeta)
%
% Figuras de paper\main_rev2.tex y CSV que leen (paper\data\, que NO se toca):
%   fig:FLL_Outkontactor (activa)  Falla_A-B_ellip.csv, Falla_B-C_ellip.csv, Falla_A-C_ellip.csv
%                                  columnas t, R_v, r_v, a_v, R_i, r_i, a_i
%   fig:F_ABC (aparcada, \iffalse) Falla_A-B-C.csv (t, va, vb, vc, ia, ib, ic) y
%                                  Falla_A-B-C_ellip.csv (R_v, r_v, R_i, r_i)
%   Falla_A-B.csv existe en paper\data pero ninguna figura lo lee.
%
% Lo que se ha podido deducir de los CSV antiguos (el script que los genero no esta):
%   - *_ellip.csv: 557 filas; el eje t es exactamente imresize(linspace(0,0.45,4501)',[557 1])
%     (reduccion bicubica de una serie de 4501 valores, uno por muestra de un tramo de 0.45 s
%     a 10 kHz); semiejes en unidades tales que el circulo de prefalta vale 1.18 (el teorico
%     es sqrt(3/2) = 1.2247); angulos en radianes en [0, pi), puestos a 0 cuando la curva es
%     casi circular; tension y corriente normalizadas por separado.
%   - Falla_A-B.csv y Falla_A-B-C.csv: 226 filas, una muestra de cada 20 (2 ms) de un tramo de
%     0.45 s; tensiones de fase de conv1 SIN referir, en pu de 326.6 V, y corrientes de conv1 en
%     pu de 90*sqrt(2) A, las dos con el orden de canales [1 3 2] (no el [3 2 1] de la ronda 2).
%     Este script localiza el registro y la muestra inicial por busqueda exhaustiva.
%   - La ventana, el suavizado y la normalizacion exactos de los *_ellip.csv no se han podido
%     reproducir (PENDIENTE). La senal de tension coincide con meas_box2 o con conv1 referida
%     al lado de red, no con conv1 sin referir.
%
% Como se generan los CSV nuevos (mismo nombre, mismas columnas, mismo numero de filas):
%   - Tramo de 4501 muestras (0.45 s). A-B y A-B-C: el mismo tramo que los CSV de ondas
%     antiguos (identificado de forma exacta). B-C y A-C: registro y muestra inicial estimados
%     ajustando la serie antigua (se informa del error de ajuste de cada candidato).
%   - Tension: r.sig.conv1_ref_red (r2_lab_signals), orden [3 2 1]. Corriente: corriente de conv1
%     referida al lado de red con la misma relacion (r.Mv' * i), orden [3 2 1].
%   - Normalizacion: amplitud de la secuencia directa de prefalta (cinco primeros ciclos del
%     tramo), por separado para tension y corriente. Un sistema equilibrado da un circulo de
%     radio sqrt(3/2) = 1.2247, el valor de la tabla teorica del paper.
%   - Elipse: ga_fault_classifier_pure sobre una ventana de W = 150 muestras (0.75 ciclos, lo
%     que declara el paper para las visualizaciones) que TERMINA en cada muestra del tramo
%     (convencion de la ronda 2: marcas al final de la ventana).
%   - Angulo en radianes, en [0, pi); 0 cuando redondez b/a >= zthr (0.9): con la curva casi
%     circular el eje no esta definido.
%   - Reduccion a 557 filas: media de bloques de 8 o 9 muestras consecutivas (media axial para
%     el angulo); t es el tiempo medio del bloque, relativo al inicio del tramo.
%   - Falla_A-B.csv y Falla_A-B-C.csv: una muestra de cada 20; tension en pu de 326.6 V y
%     corriente en pu de 90*sqrt(2) A, con las senales referidas al lado de red.
% Salidas: los seis CSV, figuras_comparacion.csv (antiguo frente a nuevo, por columna),
% figuras_identificacion.csv, figuras_angulos_sim_lab.csv, figuras_comparacion.png y
% r2_final_lab_fig_log.txt. Sin azar. Solo lectura de datos y de paper\data.
%
% Ejecucion (R2024b):
%   "C:\Program Files\MATLAB\R2024b\bin\matlab.exe" -batch "addpath('<raiz del proyecto>\experiments'); r2_final_lab_fig"
p = struct('W', 150, 'stamp', 'end', 'zthr', 0.9, 'outdir', '');
for k = 1:2:numel(varargin), p.(varargin{k}) = varargin{k+1}; end
here = fileparts(mfilename('fullpath')); addpath(here); addpath(fullfile(here, '..', 'src', 'fault_classifiers'));
RAIZ = fileparts(here); R2 = fullfile(RAIZ, '_revision_TPWRD', 'ronda2'); OLD = fullfile(RAIZ, 'paper', 'data');
outdir = p.outdir; if isempty(outdir), outdir = fullfile(R2, 'resultados', 'T90', 'figuras'); end
if ~exist(outdir, 'dir'), mkdir(outdir); end
logf = fopen(fullfile(outdir, 'r2_final_lab_fig_log.txt'), 'w'); P = @(varargin) plog(logf, varargin{:});
P('r2_final_lab_fig %s MATLAB %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), version);
P('Opciones: W = %d muestras, marca de tiempo = %s de la ventana, angulo a 0 si b/a >= %.2f\n', p.W, p.stamp, p.zthr);
fs = 1e4; f = 50; cyc = fs/f; NS = 4501; NR = 557; IB = 90*sqrt(2);
sets = {'Outkontaktor', 'Autokontaktor'};
% incepcion por registro: la de r2_lab_autokontaktor (lab_severidad.csv), ya regenerada por r2_final_lab
sv = fullfile(R2, 'resultados', 'T90', 'lab', 'lab_severidad.csv'); if ~isfile(sv), sv = fullfile(R2, 'resultados', 'T20', 'G5_005', 'lab_severidad.csv'); end
S0 = readtable(sv, 'TextType', 'string'); P('Incepciones: %s\n', sv);
% rejilla antigua: imresize(linspace(0,0.45,4501)',[557 1]) coincide con ((k-0.5)*4501/557-0.5)*1e-4 salvo en las dos filas de cada extremo
ug = ((1:NR)' - 0.5)*NS/NR + 0.5;                                   % posicion (muestra del tramo, base 1) del centro de cada fila

%% =========================================================== 1. CSV de ondas antiguos: registro y tramo
P('\n1. Registro y tramo de los CSV de ondas antiguos (busqueda exhaustiva: 2 senales x 6 ordenes x 2 signos x 20 desfases)\n');
WV = {'Falla_A-B', 'ab'; 'Falla_A-B-C', 'abc'}; pp = perms(1:3); wave = struct();
for q = 1:2
    Wt = readtable(fullfile(OLD, [WV{q,1} '.csv'])); wv = [Wt.va Wt.vb Wt.vc]; wi = [Wt.ia Wt.ib Wt.ic]; n = height(Wt);
    best = struct('e', Inf);
    for is = 1:2
        fl = dir(fullfile(RAIZ, 'data', 'lab', ['Fault_emulator_' sets{is}], [WV{q,2} '1*.mat']));
        for k = 1:numel(fl)
            r = r2_lab_signals(fullfile(fl(k).folder, fl(k).name));
            sig = {'conv1 fase (filas 59:61, pu)', r.Vcpu; 'meas_box2 (filas 183:185)/326.6 V', r.Vb/r.VB};
            for s = 1:2
                for ip = 1:6
                    for sg = [1 -1]
                        for off = 0:19
                            X = sg*sig{s,2}(off+1:20:end, pp(ip,:)); m = size(X, 1); if m < n, continue, end
                            e = zeros(m - n + 1, 1);
                            for c = 1:3
                                xx = movsum(X(:,c).^2, [0 n-1]); e = e + xx(1:m-n+1) - 2*conv(X(:,c), flipud(wv(:,c)), 'valid') + sum(wv(:,c).^2);
                            end
                            [em, j] = min(e);
                            if em < best.e, best = struct('e', em, 'rec', fl(k).name(1:end-4), 'set', sets{is}, 'sig', sig{s,1}, 'perm', pp(ip,:), 'sg', sg, 'k0', off + 1 + 20*(j-1)); end
                        end
                    end
                end
            end
        end
    end
    r = r2_lab_signals(fullfile(RAIZ, 'data', 'lab', ['Fault_emulator_' best.set], [best.rec '.mat'])); idx = best.k0 + 20*(0:n-1)';
    Ip = reshape(r.Ic(idx, best.perm), [], 1); g = Ip \ wi(:); best.i_res = norm(wi(:) - g*Ip)/norm(wi(:)); best.i_base = 1/g;
    best.rms = sqrt(max(best.e, 0)/(3*n)); js = find(S0.registro == best.rec & S0.conjunto == best.set, 1); best.k_on = S0.k_on(js);
    P('  %s.csv = registro %s (%s), %s, orden de canales %s, signo %+d, muestras %d + 20k (t = %.4f s), error rms %.2g pu;\n     corrientes = conv1 (filas 54:56) con el mismo orden / %.2f A (residuo %.1g); incepcion a %.4f s del inicio del tramo\n', ...
        WV{q,1}, best.rec, best.set, best.sig, mat2str(best.perm), best.sg, best.k0, (best.k0-1)/fs, best.rms, best.i_base, best.i_res, (best.k_on - best.k0)/fs);
    wave.(WV{q,2}) = best;
end

%% =========================================================== 2. B-C y A-C: registro y tramo estimados
P('\n2. Registro y tramo de los *_ellip.csv antiguos por ajuste de la serie (ventana centrada de %d muestras, senal adoptada)\n', p.W);
P('   error = rms(r_v/R_v antiguo - b/a) + rms(R_i antiguo/prefalta - a_i/prefalta); desplazamiento d = muestras entre el inicio del tramo y la incepcion\n');
FG = {'Falla_A-B_ellip', 'ab', 'A-B'; 'Falla_B-C_ellip', 'bc', 'B-C'; 'Falla_A-C_ellip', 'ac', 'C-A'; 'Falla_A-B-C_ellip', 'abc', 'A-B-C'};
idrows = {}; seg = struct();
for q = 1:4
    O = readtable(fullfile(OLD, [FG{q,1} '.csv'])); ro = O.r_v./O.R_v; io = O.R_i/median(O.R_i(O.t < 0.10));
    cand = S0(S0.tipo == FG{q,2}, :); bestq = struct('e', Inf);
    for k = 1:height(cand)
        r = r2_lab_signals(fullfile(RAIZ, 'data', 'lab', ['Fault_emulator_' char(cand.conjunto(k))], [char(cand.registro(k)) '.mat']));
        [V, I] = signals(r); k_on = cand.k_on(k); N = size(V, 1);
        lo = max(p.W, k_on - 1900); hi = min(N - p.W, k_on + 3400); ii = (lo:4:hi)';
        [av, bv] = slide(V, ii, p.W, 'center', fs, f); ai = slide(I, ii, p.W, 'center', fs, f);
        dd = (1100:2:1900)'; ev = nan(size(dd)); ei = ev;
        for u = 1:numel(dd)
            x = k_on - dd(u) + ug - 1; if x(1) < ii(1) || x(end) > ii(end), continue, end
            rn = interp1(ii, bv./av, x); an = interp1(ii, ai, x); an = an/median(an(O.t < 0.10));
            ev(u) = sqrt(mean((ro - rn).^2)); ei(u) = sqrt(mean((io - an).^2));
        end
        [em, u] = min(ev + ei);
        w = struct('archivo', FG{q,1}, 'registro', char(cand.registro(k)), 'conjunto', char(cand.conjunto(k)), 'P_pre_kW', cand.P_pre_kW(k), 'k_on', k_on, ...
            'd_muestras', dd(u), 'k0', k_on - dd(u), 'error_redondez', ev(u), 'error_corriente', ei(u), 'error_total', em, 'elegido', false, 'metodo', 'ajuste de la serie');
        idrows{end+1} = w; %#ok<AGROW>
        if em < bestq.e, bestq = struct('e', em, 'i', numel(idrows)); end
    end
    idrows{bestq.i}.elegido = true; w = idrows{bestq.i};
    P('  %-18s mejor candidato: %s (%s), P = %.1f kW, d = %d, error %.3f (redondez %.3f, corriente %.3f)\n', FG{q,1}, w.registro, w.conjunto, w.P_pre_kW, w.d_muestras, w.error_total, w.error_redondez, w.error_corriente);
    seg.(FG{q,2}) = struct('rec', w.registro, 'set', w.conjunto, 'k0', w.k0, 'k_on', w.k_on, 'metodo', 'estimado por ajuste de la serie antigua', 'err', w.error_total);
end
ID = struct2table([idrows{:}]');
for c = string(ID.Properties.VariableNames), if iscell(ID.(c)), ID.(c) = string(ID.(c)); elseif ischar(ID.(c)), ID.(c) = string(cellstr(ID.(c))); end, end
% donde hay CSV de ondas, manda la identificacion exacta
for tp = {'ab', 'abc'}
    b = wave.(tp{1}); s = seg.(tp{1});
    P('  %s: identificacion exacta por las ondas = %s (%s), muestra inicial %d; el ajuste de la serie daba %s (%s), muestra inicial %d (diferencia %d muestras)\n', ...
        tp{1}, b.rec, b.set, b.k0, s.rec, s.set, s.k0, s.k0 - b.k0);
    ID.metodo(ID.archivo == FG{strcmp(FG(:,2), tp{1}), 1} & ID.elegido) = "ajuste de la serie (comprobacion; se usa el tramo exacto de las ondas)";
    seg.(tp{1}) = struct('rec', b.rec, 'set', b.set, 'k0', b.k0, 'k_on', b.k_on, 'metodo', 'exacto (coincide con el CSV de ondas antiguo)', 'err', s.err);
    ID.ajuste_coincide_con_ondas(ID.archivo == FG{strcmp(FG(:,2), tp{1}), 1} & ID.elegido) = strcmp(s.rec, b.rec) && strcmp(s.set, b.set);
    ID.k0_exacto_ondas(ID.archivo == FG{strcmp(FG(:,2), tp{1}), 1} & ID.elegido) = b.k0;
end
writetable(ID, fullfile(outdir, 'figuras_identificacion.csv'));

%% =========================================================== 3. CSV nuevos
P('\n3. CSV nuevos en %s\n', outdir);
cmp = {}; NEW = struct(); OLDT = struct();
for q = 1:4
    s = seg.(FG{q,2}); r = r2_lab_signals(fullfile(RAIZ, 'data', 'lab', ['Fault_emulator_' s.set], [s.rec '.mat']));
    [V, I, Vpu, Ipu] = signals(r); ix = (s.k0 : s.k0 + NS - 1)'; tt = (0:NS-1)'/fs;
    % normalizacion: amplitud de secuencia directa en los cinco primeros ciclos del tramo
    nv = seq1(V(ix(1:5*cyc),:), tt(1:5*cyc), f); ni = seq1(I(ix(1:5*cyc),:), tt(1:5*cyc), f);
    [av, bv, thv] = slide(V/nv, ix, p.W, p.stamp, fs, f); [ai, bi, thi] = slide(I/ni, ix, p.W, p.stamp, fs, f);
    blk = min(NR, floor(((1:NS)' - 0.5)*NR/NS) + 1);
    T = table(accumarray(blk, tt, [NR 1], @mean), accumarray(blk, av, [NR 1], @mean), accumarray(blk, bv, [NR 1], @mean), blockang(blk, thv, bv./av, p.zthr, NR), ...
        accumarray(blk, ai, [NR 1], @mean), accumarray(blk, bi, [NR 1], @mean), blockang(blk, thi, bi./ai, p.zthr, NR), 'VariableNames', {'t','R_v','r_v','a_v','R_i','r_i','a_i'});
    writetable(T, fullfile(outdir, [FG{q,1} '.csv'])); NEW.(FG{q,2}) = T;
    O = readtable(fullfile(OLD, [FG{q,1} '.csv'])); OLDT.(FG{q,2}) = O;
    P('  %-22s registro %s (%s), muestras %d a %d (t = %.4f a %.4f s), incepcion a %.4f s del inicio; tramo %s\n     base de tension %.4f pu de 326.6 V; base de corriente %.2f A de pico (%.2f A rms)\n', ...
        [FG{q,1} '.csv'], s.rec, s.set, ix(1), ix(end), (ix(1)-1)/fs, (ix(end)-1)/fs, (s.k_on - s.k0)/fs, s.metodo, nv, ni*IB, ni*IB/sqrt(2));
    pre = O.t < 0.10; fl = O.t > 0.17 & O.t < 0.26;
    for c = {'R_v','r_v','a_v','R_i','r_i','a_i'}
        o = O.(c{1}); x = T.(c{1}); isang = c{1}(1) == 'a';
        if isang, d = mod(x - o + pi/2, pi) - pi/2; mo = @(v, m) axmed(v(m)); else, d = x - o; mo = @(v, m) median(v(m)); end
        w = struct('archivo', [FG{q,1} '.csv'], 'figura', tern(q < 4, 'fig:FLL_Outkontactor (activa)', 'fig:F_ABC (aparcada)'), 'columna', c{1}, 'usada_en_figura', q < 4 || ~isang, 'filas', NR, 'dif_t_max_ms', 1e3*max(abs(T.t - O.t)), ...
            'antiguo_prefalta', mo(o, pre), 'nuevo_prefalta', mo(x, pre), 'antiguo_falta', mo(o, fl), 'nuevo_falta', mo(x, fl), 'dif_falta', mo(x, fl) - mo(o, fl), ...
            'dif_falta_rel_pct', 100*(mo(x, fl) - mo(o, fl))/max(abs(mo(o, fl)), eps), 'dif_abs_max', max(abs(d)), 'dif_rms', sqrt(mean(d.^2)), 'dif_max_rel_pct', 100*max(abs(d))/max(abs(o)), ...
            'dif_abs_max_regimen', max(abs(d(fl))), 'dif_max_rel_regimen_pct', 100*max(abs(d(fl)))/max(abs(o(fl))), 'registro_nuevo', s.rec, 'k0_nuevo', s.k0, 'tramo', s.metodo, 'sustituir', q < 4);
        cmp{end+1} = w; %#ok<AGROW>
    end
end
% ondas (una muestra de cada 20), senales referidas al lado de red
for q = 1:2
    s = seg.(WV{q,2}); r = r2_lab_signals(fullfile(RAIZ, 'data', 'lab', ['Fault_emulator_' s.set], [s.rec '.mat'])); [~, ~, Vpu, Ipu] = signals(r);
    idx = s.k0 + 20*(0:225)'; T = table((0:225)'*0.002, Vpu(idx,1), Vpu(idx,2), Vpu(idx,3), Ipu(idx,1), Ipu(idx,2), Ipu(idx,3), 'VariableNames', {'t','va','vb','vc','ia','ib','ic'});
    writetable(T, fullfile(outdir, [WV{q,1} '.csv'])); O = readtable(fullfile(OLD, [WV{q,1} '.csv']));
    P('  %-22s registro %s (%s), muestras %d + 20k; tension conv1 referida al lado de red en pu de 326.6 V, corriente referida en pu de %.2f A\n', [WV{q,1} '.csv'], s.rec, s.set, s.k0, IB);
    for c = {'va','vb','vc','ia','ib','ic'}
        o = O.(c{1}); x = T.(c{1}); d = x - o;
        w = struct('archivo', [WV{q,1} '.csv'], 'figura', tern(q == 2, 'fig:F_ABC (aparcada)', 'ninguna'), 'columna', c{1}, 'usada_en_figura', q == 2, 'filas', height(T), 'dif_t_max_ms', 1e3*max(abs(T.t - O.t)), ...
            'antiguo_prefalta', max(abs(o(O.t < 0.10))), 'nuevo_prefalta', max(abs(x(O.t < 0.10))), 'antiguo_falta', max(abs(o(O.t > 0.17 & O.t < 0.26))), 'nuevo_falta', max(abs(x(O.t > 0.17 & O.t < 0.26))), 'dif_falta', NaN, ...
            'dif_falta_rel_pct', NaN, 'dif_abs_max', max(abs(d)), 'dif_rms', sqrt(mean(d.^2)), 'dif_max_rel_pct', 100*max(abs(d))/max(abs(o)), 'dif_abs_max_regimen', NaN, 'dif_max_rel_regimen_pct', NaN, ...
            'registro_nuevo', s.rec, 'k0_nuevo', s.k0, 'tramo', s.metodo, 'sustituir', false);
        cmp{end+1} = w; %#ok<AGROW>
    end
end
C = struct2table([cmp{:}]');
for c = string(C.Properties.VariableNames), if iscell(C.(c)), C.(c) = string(C.(c)); elseif ischar(C.(c)), C.(c) = string(cellstr(C.(c))); end, end
writetable(C, fullfile(outdir, 'figuras_comparacion.csv'));
P('\n4. Antiguo frente a nuevo (prefalta t < 0.10 s; falta 0.17 < t < 0.26 s; angulos en rad, diferencias modulo pi; en las ondas, pico en vez de mediana)\n');
P('  %-22s %-4s | prefalta ant. / nuevo | falta ant. / nuevo (dif. rel.) | dif. max (rel. al max. antiguo) | dif. rms | dif. max en regimen\n', 'archivo', 'col');
for i = 1:height(C)
    P('  %-22s %-4s | %8.4f / %8.4f | %8.4f / %8.4f (%+6.1f %%) | %8.4f (%6.1f %%) | %8.4f | %8.4f\n', C.archivo(i), C.columna(i), C.antiguo_prefalta(i), C.nuevo_prefalta(i), C.antiguo_falta(i), C.nuevo_falta(i), ...
        C.dif_falta_rel_pct(i), C.dif_abs_max(i), C.dif_max_rel_pct(i), C.dif_rms(i), C.dif_abs_max_regimen(i));
end

%% =========================================================== 5. Angulos: simulacion RL (paper), laboratorio antiguo y nuevo
ang = {}; SM = {'FallaRL_A-B_Pex_ellip', 'ab', 'A-B'; 'FallaRL_B-C_Pex_ellip', 'bc', 'B-C'; 'FallaRL_A-C_Pex_ellip', 'ac', 'C-A'};
for q = 1:3
    w = struct('falta', SM{q,3}); fsim = fullfile(OLD, [SM{q,1} '.csv']);
    if isfile(fsim)
        Sx = readtable(fsim); m = Sx.r_v./Sx.R_v < 0.9; j1 = find(m, 1); j2 = find(m, 1, 'last'); m = false(height(Sx), 1); m(round(j1 + 0.25*(j2-j1)) : round(j1 + 0.75*(j2-j1))) = true;
        w.sim_RL_a_v_deg = axmed(Sx.a_v(m))*180/pi; w.sim_RL_a_i_deg = axmed(Sx.a_i(m))*180/pi;
    else, w.sim_RL_a_v_deg = NaN; w.sim_RL_a_i_deg = NaN;
    end
    O = OLDT.(SM{q,2}); T = NEW.(SM{q,2}); fl = O.t > 0.17 & O.t < 0.26;
    w.lab_antiguo_a_v_deg = axmed(O.a_v(fl))*180/pi; w.lab_antiguo_a_i_deg = axmed(O.a_i(fl))*180/pi; w.lab_nuevo_a_v_deg = axmed(T.a_v(fl))*180/pi; w.lab_nuevo_a_i_deg = axmed(T.a_i(fl))*180/pi;
    w.lab_antiguo_b_sobre_a_v = median(O.r_v(fl)./O.R_v(fl)); w.lab_nuevo_b_sobre_a_v = median(T.r_v(fl)./T.R_v(fl));
    ang{end+1} = w; %#ok<AGROW>
end
A = struct2table([ang{:}]'); A.falta = string(cellstr(A.falta)); writetable(A, fullfile(outdir, 'figuras_angulos_sim_lab.csv'));
P('\n5. Angulo del eje en la falta, en grados (simulacion RL de Pex del paper; laboratorio antiguo; laboratorio nuevo):\n');
for i = 1:height(A)
    P('  %s: tension %.1f / %.1f / %.1f; corriente %.1f / %.1f / %.1f\n', A.falta(i), A.sim_RL_a_v_deg(i), A.lab_antiguo_a_v_deg(i), A.lab_nuevo_a_v_deg(i), A.sim_RL_a_i_deg(i), A.lab_antiguo_a_i_deg(i), A.lab_nuevo_a_i_deg(i));
end

%% =========================================================== 6. Figura de control (antiguo a trazos, nuevo continuo)
try
    fg = figure('Visible', 'off', 'Color', [0.99 0.98 0.95], 'Position', [50 50 1100 800]); col = [0.10 0.20 0.45; 0.90 0.45 0.10; 0.25 0.25 0.25]; tl = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    pan = {'R_v', 'r_v', 'Semiejes de tension'; 'R_i', 'r_i', 'Semiejes de corriente'; 'a_v', '', 'Angulo de tension (rad)'; 'a_i', '', 'Angulo de corriente (rad)'}; tp = {'ab', 'bc', 'ac'};
    for u = 1:4
        ax = nexttile(tl); hold(ax, 'on'); set(ax, 'Color', [0.99 0.98 0.95], 'XLim', [0 0.45], 'GridColor', [0.6 0.7 0.85], 'XGrid', 'on', 'YGrid', 'on');
        for q = 1:3
            O = OLDT.(tp{q}); T = NEW.(tp{q});
            plot(ax, O.t, O.(pan{u,1}), '--', 'Color', col(q,:), 'LineWidth', 0.8); plot(ax, T.t, T.(pan{u,1}), '-', 'Color', col(q,:), 'LineWidth', 1.2);
            if ~isempty(pan{u,2}), plot(ax, O.t, O.(pan{u,2}), ':', 'Color', col(q,:), 'LineWidth', 0.8); plot(ax, T.t, T.(pan{u,2}), '-', 'Color', col(q,:), 'LineWidth', 0.6); end
        end
        ylabel(ax, pan{u,3}); if u == 4, xlabel(ax, 'Tiempo (s)'); end
    end
    title(tl, 'Figura de laboratorio: CSV antiguos (a trazos) y nuevos con conv1 referida al lado de red (continuo). Azul A-B, naranja B-C, gris C-A');
    exportgraphics(fg, fullfile(outdir, 'figuras_comparacion.png'), 'Resolution', 130); close(fg);
catch ME
    P('No se pudo generar figuras_comparacion.png: %s\n', ME.message);
end
fclose(logf);
R = struct('comparacion', C, 'identificacion', ID, 'angulos', A, 'tramos', seg, 'ondas', wave);
end

%% ======================================================================
function [V, I, Vpu, Ipu] = signals(r)
% Senal adoptada: tension de conv1 referida al lado de red (pu de 326.6 V, orden [3 2 1]) y
% corriente de conv1 referida al lado de red con la misma relacion (i_red = Mv' * i_conv).
Vpu = r.sig.conv1_ref_red; x = (r.Mv'*r.Ic')'/(90*sqrt(2)); Ipu = x(:, r.PERM); V = Vpu; I = Ipu;
end
function a1 = seq1(v, t, f)
% Amplitud (pico) de la secuencia directa de 50 Hz por minimos cuadrados
X = [cos(2*pi*f*t), sin(2*pi*f*t), ones(numel(t), 1)]; c = X \ v; ph = (c(1,:) - 1j*c(2,:)).';
a = exp(2j*pi/3); a1 = abs([1 a a^2]*ph/3);
end
function [a, b, th] = slide(X, idx, W, stamp, fs, f)
% Parametros de la elipse (los del clasificador) en una ventana de W muestras por cada indice:
% 'end' = la ventana termina en la muestra; 'center' = centrada en ella. Angulo en radianes.
n = numel(idx); a = nan(n, 1); b = a; th = a; ws = warning('off', 'all'); cl = onCleanup(@() warning(ws));
if strcmp(stamp, 'end'), o1 = -(W - 1); else, o1 = -floor(W/2); end
for j = 1:n
    i1 = idx(j) + o1; if i1 < 1 || i1 + W - 1 > size(X, 1), continue, end
    try
        [~, q] = ga_fault_classifier_pure(X(i1:i1+W-1,:), fs, f);
        if ~isempty(q.x_T), a(j) = q.semi_major; b(j) = q.semi_minor; th(j) = q.angle_deg*pi/180; end
    catch
    end
end
a = fillmissing(a, 'nearest'); b = fillmissing(b, 'nearest'); th = fillmissing(th, 'nearest');
end
function m = blockang(blk, th, rr, zthr, NR)
% Angulo por bloque: media axial de las ventanas no circulares; 0 si la mayoria son circulares
m = zeros(NR, 1);
for k = 1:NR
    s = blk == k; ok = s & rr < zthr;
    if nnz(ok) > nnz(s)/2, m(k) = mod(angle(mean(exp(2j*th(ok))))/2, pi); end
end
end
function m = axmed(th), th = th(isfinite(th)); m = mod(angle(mean(exp(2j*th)))/2, pi); end
function v = tern(c, a, b), if c, v = a; else, v = b; end, end
function plog(fid, varargin), fprintf(varargin{:}); fprintf(fid, varargin{:}); end
