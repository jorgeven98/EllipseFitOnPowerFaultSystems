%% R2_DATASET_INVENTORY  (T72, peticion R4.2c del revisor 4)
% Inventario medido de los 18 .mat de data/simulation_curated: tipo de falta
% (por caida de tension, no por el nombre), instante de incepcion y de despeje,
% angulo de incidencia, tension, corriente y potencia prefalta de Pex, DG1 y
% DG2, y numero de ventanas del clasificador. Complementa a
% r2_dataset_inventory_slx.py (valores de los bloques de falta del modelo).
%
% Todo sale de los datos. Lo que no se puede medir (Rf, Rg, Lf, SwitchTimes
% por archivo) NO se escribe aqui: se toma del modelo y es DEDUCIDO.
%
% Ejecutar desde cualquier carpeta:
%   "C:\Program Files\MATLAB\R2024b\bin\matlab.exe" -batch "run('<RAIZ>\experiments\r2_dataset_inventory.m')"
% Salidas (resultados\T72): dataset.csv (la columna used_in da el uso en los scripts), inventario_log.txt
%
% Convenciones:
%  - va = V*cos(theta); theta en grados en (-180,180]; theta = 0 es un pico
%    positivo y +-90 un cruce por cero. Bus de medida: tensiones fase-neutro en
%    bornes del convertidor indicado (Pex = nodo 4 del esquema).
%  - Resolucion temporal 0.25 ms = 4.5 grados electricos a 50 Hz.
%  - Marcas temporales de la incepcion: primera muestra con residuo > umbral
%    frente a la senoide prefalta ajustada (t_ini). La incepcion real esta en
%    (t_ini - Ts, t_ini]. t_sw es el SwitchTimes del modelo mas cercano.

clear; clc;
RAIZ = fileparts(fileparts(mfilename('fullpath')));
OUT  = fullfile(RAIZ, '_revision_TPWRD', 'ronda2', 'resultados', 'T72');
if ~exist(OUT, 'dir'), mkdir(OUT); end
DATA = fullfile(RAIZ, 'data', 'simulation_curated');
if exist(fullfile(OUT, 'inventario_log.txt'), 'file'), delete(fullfile(OUT, 'inventario_log.txt')); end
diary(fullfile(OUT, 'inventario_log.txt')); diary on;
fprintf('r2_dataset_inventory %s  MATLAB %s\n', datestr(now, 31), version);

F0 = 50; W = 20; STRIDE = 5;
% nombre, etiqueta esperada por el nombre (convencion del paper), SwitchTimes
% candidatos del modelo (DEDUCIDO, de slx_blocks.csv)
files = {
 'DataGridRFAB.mat',   'A-B';   'DataGridRFBC.mat',   'B-C';   'DataGridRFCA.mat',   'C-A';
 'DataGridRL_AB.mat',  'A-B';   'DataGridRL_BC.mat',  'B-C';   'DataGridRL_AC.mat',  'C-A';
 'DataGridRL_ABC.mat', 'A-B-C';
 'DataGrid_AB.mat',    'A-B';   'DataGrid_BC.mat',    'B-C';   'DataGrid_AC.mat',    'C-A';
 'DataGrid_ABC.mat',   'A-B-C';
 'DataGrid_AG.mat',    'A-G';   'DataGrid_BG.mat',    'B-G';   'DataGrid_CG.mat',    'C-G';
 'DataGrid_ABG.mat',   'A-B-G'; 'DataGrid_BCG.mat',   'B-C-G'; 'DataGrid_ACG.mat',   'C-A-G';
 'DataGrid_ABCG.mat',  'A-B-C-G'};
SW_CAND = [3.000 3.023 2.023];

% ---- uso de cada archivo en los scripts de experimentos (busqueda textual)
expdir = fullfile(RAIZ, 'experiments');
sc = dir(fullfile(expdir, '*.m'));
uso = strings(size(files,1), 1);
for i = 1:size(files,1)
    u = {};
    for j = 1:numel(sc)
        if strncmp(sc(j).name, 'r2_', 3), continue; end
        txt = fileread(fullfile(expdir, sc(j).name));
        if contains(txt, ['''' files{i,1} ''''])
            u{end+1} = sc(j).name; %#ok<AGROW>
        end
    end
    uso(i) = strjoin(string(u), ' ; ');
end

conv = {'Pex', 'DG1', 'DG2'};
R = struct([]); SWP = struct([]);
for i = 1:size(files, 1)
    fn = files{i,1}; fprintf('\n[%02d/%02d] %s (nombre: %s)\n', i, size(files,1), fn, files{i,2});
    d = load(fullfile(DATA, fn));
    G = d.GridConvs;
    [vP, tP] = getV(G, 'Pex');
    Ts = median(diff(tP)); fs = 1/Ts; N = numel(tP);
    spc = round(fs/F0);   % muestras por ciclo (80)
    assert(abs(spc - fs/F0) < 1e-9, 'fs no es multiplo de 50 Hz');
    w = 2*pi*F0;

    % --- 1. incepcion aproximada por diferencia de periodo (v(k)-v(k-1ciclo))
    kref = (tP >= 1.7 & tP <= 1.9);
    Vpk = max(abs(vP(kref, :)), [], 1);                 % pico prefalta por fase
    dper = [zeros(spc,3); abs(vP(spc+1:end,:) - vP(1:end-spc,:))] ./ Vpk;
    kr = find(max(dper,[],2) > 0.05 & tP >= 1.9, 1, 'first');
    % --- 2. ajuste senoidal prefalta (4 ciclos que acaban 2 ciclos antes) y residuo
    k2 = kr - 2*spc; k1 = k2 - 4*spc + 1;
    tt = tP(k1:k2);
    A = [cos(w*tt), sin(w*tt), ones(size(tt))];
    cf = zeros(3,3); sig = zeros(1,3);
    for p = 1:3
        cf(:,p) = A \ vP(k1:k2, p);
        sig(p) = std(vP(k1:k2,p) - A*cf(:,p));
    end
    amp = hypot(cf(1,:), cf(2,:));                       % amplitud ajustada
    phi = atan2(-cf(2,:), cf(1,:));                      % v = amp*cos(w t + phi)
    kk = (kr-2*spc):(kr+2*spc);
    Af = [cos(w*tP(kk)), sin(w*tP(kk)), ones(numel(kk),1)];
    res = abs(vP(kk,:) - Af*cf) ./ amp;                  % residuo en fraccion del pico
    thr = max(0.02, 6*max(sig./amp));
    hit = max(res, [], 2) > thr;
    j = find(hit(1:end-2) & hit(2:end-1) & hit(3:end), 1, 'first');   % 3 muestras seguidas
    kini = kk(j); t_ini = tP(kini);
    [~, ic] = min(abs(SW_CAND - t_ini)); t_sw = SW_CAND(ic);
    fprintf('  t_ini medido = %.4f s (umbral residuo %.3f pico), SwitchTimes candidato %.3f s (dif %.2f ms)\n', ...
        t_ini, thr, t_sw, 1e3*(t_ini - t_sw));

    % --- 3. despeje: segunda vez que la diferencia de periodo supera el umbral
    % Cada polo del bloque de falta abre en su paso por cero de corriente, asi que
    % las fases se despejan en instantes distintos. Se toma el despeje de la
    % ultima fase en falta: fases con caida (rms del 2o ciclo de falta < 0.8 del
    % prefalta) cuya envolvente (maximo de |v| en las 20 muestras previas) vuelve
    % a superar 0.85 de su pico prefalta. Resolucion del orden de un cuarto de ciclo.
    seg2 = kini+spc : kini+2*spc-1;
    rr2 = rms(vP(seg2,:),1) ./ rms(vP(k1:k2,:),1);
    fph = find(rr2 < 0.8);
    env = movmax(abs(vP), [19 0], 1);
    kc = kini + 2*spc;
    for p_ = fph
        kk_ = find((1:N)' > kini + 2*spc & env(:,p_) > 0.85*Vpk(p_), 1, 'first');
        kc = max(kc, kk_);
    end
    how = 'envolvente de las fases en falta';
    kfin = kc; t_fin = tP(kfin);
    fprintf('  despeje (%s) = %.4f s; duracion %.1f ms\n', how, t_fin, 1e3*(t_fin - t_ini));

    % --- 4. tipo de falta por caida de tension (Pex) en la falta estacionaria
    kf1 = kini + spc; kf2 = kini + 2*spc - 1;   % 2o ciclo de falta
    Vf = vP(kf1:kf2,:);  Vp = vP(k1:k2,:);
    rms_f = rms(Vf,1); rms_p = rms(Vp,1);
    ratio = rms_f ./ rms_p;
    v0f = rms(mean(Vf,2)); v0p = rms(mean(Vp,2));
    v0ratio = v0f ./ mean(rms_p);
    inferred = infer_type(ratio, v0ratio);
    ok = strcmp(inferred, files{i,2});
    fprintf('  Pex rms falta/pre = [%.3f %.3f %.3f], V0/Vfase_pre = %.3f -> %s (nombre %s) %s\n', ...
        ratio, v0ratio, inferred, files{i,2}, ternary(ok, 'OK', '*** NO COINCIDE ***'));

    % --- 5. angulo de incidencia: fase de cada tension en t_sw y en t_ini
    th_sw  = wrapdeg(rad2deg(w*t_sw  + phi));
    th_ini = wrapdeg(rad2deg(w*t_ini + phi));
    vpu_ini = cos(deg2rad(th_ini));
    % lazo en falta: para L-L, angulo de la tension compuesta del par que cae
    pair = pair_of(inferred);
    th_loop_sw = NaN;
    if numel(pair) == 2
        ph = [phi(pair(1)) amp(pair(1)); phi(pair(2)) amp(pair(2))];
        z = amp(pair(1))*exp(1j*(w*t_sw + phi(pair(1)))) - amp(pair(2))*exp(1j*(w*t_sw + phi(pair(2))));
        th_loop_sw = wrapdeg(rad2deg(angle(z)));
    end
    % angulo de referencia de la falta: lazo (L-L), fase en falta (con tierra), fase A (trifasica)
    if ~isnan(th_loop_sw), th_ref = th_loop_sw; th_what = 'loop';
    elseif any(contains(inferred, 'G')) && ~strcmp(inferred, 'A-B-C-G')
        pf = find(ratio < 0.3, 1, 'first'); th_ref = th_sw(pf); nm3 = 'ABC'; th_what = ['phase ' nm3(pf)];
    else, th_ref = th_sw(1); th_what = 'phase A'; end
    fprintf('  theta(fase a,b,c) en t_sw = [%.1f %.1f %.1f] grados; en t_ini = [%.1f %.1f %.1f]\n', th_sw, th_ini);

    % --- 6. punto de operacion prefalta (2 ciclos justo antes de t_ini) por convertidor
    kp2 = kini - 1; kp1 = kp2 - 2*spc + 1;
    op = struct();
    for c = 1:3
        [vv, tv] = getV(G, conv{c}); [ii, ~] = getI(G, conv{c});
        assert(abs(tv(1)-tP(1))<1e-12 && numel(tv)==N, 'ejes de tiempo distintos');
        tw = tP(kp1:kp2); E = exp(-1j*w*tw);
        Vph = 2/numel(tw) * (E.' * vv(kp1:kp2,:));       % fasor pico, fila
        Iph = 2/numel(tw) * (E.' * ii(kp1:kp2,:));
        S = 0.5 * sum(Vph .* conj(Iph));                  % P + jQ (signo de Io segun el modelo)
        op.(conv{c}).Vrms = mean(abs(Vph))/sqrt(2);
        op.(conv{c}).Irms = mean(abs(Iph))/sqrt(2);
        op.(conv{c}).P = real(S); op.(conv{c}).Q = imag(S);
        op.(conv{c}).Vunb = (max(abs(Vph)) - min(abs(Vph))) / mean(abs(Vph));
        % ratio de tension de falta (minimo rms/prefalta) de cada convertidor
        Vf_c = vv(kf1:kf2,:); Vp_c = vv(k1:k2,:);
        op.(conv{c}).ratio = rms(Vf_c,1) ./ rms(Vp_c,1);
        op.(conv{c}).Ifault_pk = max(abs(ii(kf1:kf2,:)), [], 'all');
        fprintf('  %s prefalta: V=%.1f V rms, I=%.2f A rms, P=%.0f W, Q=%.0f var; caida rms [%.2f %.2f %.2f]; I pico en falta %.1f A\n', ...
            conv{c}, op.(conv{c}).Vrms, op.(conv{c}).Irms, op.(conv{c}).P, op.(conv{c}).Q, ...
            op.(conv{c}).ratio, op.(conv{c}).Ifault_pk);
    end
    % Vdc de Pex
    vdc_pre = mean(G.Pex.Vdc.Data(kp1:kp2)); vdc_min = min(G.Pex.Vdc.Data(kf1:kf2));

    % --- 7. ventanas del clasificador (W = 20, paso 5)
    nfault_s = kfin - kini;
    nwin_fault_all = floor((nfault_s - W)/STRIDE) + 1;      % ventanas totalmente dentro de la falta medida
    idx0 = find(tP >= 3.0, 1, 'first'); idx1 = find(tP <= 3.16, 1, 'last');
    if isempty(idx0) || isempty(idx1) || idx1 < idx0
        nwin_acc = 0;
    else
        nwin_acc = floor((idx1 - idx0 + 1 - W)/STRIDE) + 1;
    end
    nwin_pre_total = floor((kini - 1 - W)/STRIDE) + 1;

    % --- 8. guardar fila
    r = struct();
    r.file = string(fn); r.name_label = string(files{i,2}); r.inferred_label = string(inferred);
    r.name_matches_data = ok;
    r.N = N; r.fs_Hz = fs; r.duration_s = tP(end) - tP(1);
    r.t_ini_s = t_ini; r.t_sw_candidate_s = t_sw; r.t_ini_minus_sw_ms = 1e3*(t_ini - t_sw);
    r.t_clear_s = t_fin; r.fault_duration_ms = 1e3*(t_fin - t_ini);
    r.pex_rms_ratio_a = ratio(1); r.pex_rms_ratio_b = ratio(2); r.pex_rms_ratio_c = ratio(3);
    r.pex_v0_over_vph = v0ratio;
    r.theta_a_sw_deg = th_sw(1); r.theta_b_sw_deg = th_sw(2); r.theta_c_sw_deg = th_sw(3);
    r.theta_loop_sw_deg = th_loop_sw; r.theta_ref_deg = th_ref; r.theta_ref_what = string(th_what);
    r.va_pu_at_ini = vpu_ini(1); r.vb_pu_at_ini = vpu_ini(2); r.vc_pu_at_ini = vpu_ini(3);
    for c = 1:3
        n = conv{c};
        r.([n '_Vrms_pre_V']) = op.(n).Vrms; r.([n '_Irms_pre_A']) = op.(n).Irms;
        r.([n '_P_pre_W']) = op.(n).P; r.([n '_Q_pre_var']) = op.(n).Q;
        r.([n '_Vmin_ratio']) = min(op.(n).ratio); r.([n '_Ifault_peak_A']) = op.(n).Ifault_pk;
    end
    r.Pex_Vdc_pre_V = vdc_pre; r.Pex_Vdc_min_fault_V = vdc_min;
    r.n_win_W20_in_measured_fault = nwin_fault_all;
    r.n_win_W20_in_analysis_region_3p00_3p16 = nwin_acc;
    r.n_win_W20_prefault = nwin_pre_total;
    r.used_in = uso(i);
    if isempty(R), R = r; else, R(end+1) = r; end %#ok<SAGROW>
end
T = struct2table(R);
writetable(T, fullfile(OUT, 'dataset.csv'));
fprintf('\nEscrito %s (%d filas)\n', fullfile(OUT,'dataset.csv'), height(T));

% ---- resumen de rangos (los que pide el revisor)
fprintf('\n== RESUMEN ==\n');
fprintf('Incoherencias nombre/datos: %d de %d\n', sum(~T.name_matches_data), height(T));
fprintf('t_ini - SwitchTimes candidato (ms): min %.2f max %.2f\n', min(T.t_ini_minus_sw_ms), max(T.t_ini_minus_sw_ms));
fprintf('Duracion de falta medida (ms): min %.1f max %.1f\n', min(T.fault_duration_ms), max(T.fault_duration_ms));
fprintf('theta_a en t_sw (grados): ');  fprintf('%.1f ', unique(round(T.theta_a_sw_deg,1))); fprintf('\n');
diary off;

%% ---------------- funciones locales ----------------
function [v, t] = getV(G, name)
    switch name
        case {'DG1','DG2'}, ts = G.(name).Vac;
        case 'Pex',         ts = G.Pex.Vabc;
    end
    v = ts.Data; t = ts.Time;
end
function [i, t] = getI(G, name)
    ts = G.(name).Io; i = ts.Data; t = ts.Time;
end
function a = wrapdeg(a)
    a = mod(a + 180, 360) - 180;
end
function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end
function lab = infer_type(ratio, v0r)
    % ratio: rms posfalta / prefalta por fase; v0r: V0 rms / tension de fase prefalta
    names = 'ABC';
    low = ratio < 0.3;
    ground = v0r > 0.10;
    if ground
        switch sum(low)
            case 1, lab = [names(low) '-G'];
            case 2, p = find(low); lab = pairname(p, true);
            case 3, lab = 'A-B-C-G';
            otherwise, lab = 'G?';
        end
    else
        if sum(low) == 3
            lab = 'A-B-C';
        else
            [~, ord] = sort(ratio, 'ascend');
            p = sort(ord(1:2));
            if ratio(ord(3)) > 0.8 && ratio(ord(2)) < 0.8, lab = pairname(p, false);
            else, lab = 'L-L?'; end
        end
    end
end
function s = pairname(p, g)
    if isequal(p, [1 2]), s = 'A-B'; elseif isequal(p, [2 3]), s = 'B-C'; else, s = 'C-A'; end
    if g, s = [s '-G']; end
end
function p = pair_of(lab)
    p = [];
    if endsWith(lab, '-G') || strcmp(lab, 'A-B-C') || strcmp(lab, 'A-B-C-G'), return; end
    switch lab
        case 'A-B', p = [1 2]; case 'B-C', p = [2 3]; case 'C-A', p = [3 1];
    end
end
