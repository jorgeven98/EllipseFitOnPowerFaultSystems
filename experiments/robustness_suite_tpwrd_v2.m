%% ROBUSTNESS_SUITE_TPWRD_V2 — defensas corregidas + control
% Cambios respecto a v1:
%  - Control sin perturbación (tasa base de falsa operación).
%  - Defensa DC: estimación LS del DC sobre el CICLO COMPLETO anterior
%    (80 muestras, base [cos,sin,1] bien condicionada), restada de la
%    ventana de 1/4 de ciclo actual. Emula el historial disponible en un IED.
%  - Defensa transitorios: persistencia M-de-M (M ventanas consecutivas con
%    la MISMA etiqueta de falta antes de declarar). Se mide falsa operación
%    tras persistencia y el coste en retardo sobre la falta real A-B.

clear; clc;
addpath(fullfile('..', 'classifiers'));
addpath(fullfile('..', 'utils'));

fs = 4000; f = 50;
outdir = fullfile('..', 'results', 'robustness');
if ~exist(outdir, 'dir'), mkdir(outdir); end

base = load(fullfile('..', 'data', 'DataGridRL_AB.mat'));
[Vpex, tpex] = extract_converter_data(base, 'Pex');
seg = tpex >= 2.4 & tpex < 2.8;
Vh = Vpex(seg, :); th = tpex(seg); th = th - th(1);
Npeak = max(abs(Vh(:)));
segF = tpex >= 3.00 & tpex < 3.14;
Vf = Vpex(segF, :); tf_ = tpex(segF); tf_ = tf_ - tf_(1);
% instante real de incepción dentro del segmento de falta
t_incep = 3.0 - tpex(find(segF,1));  % = 0 por construcción (segF empieza en 3.00)

W = 20; STRIDE = 5; M_PERS = 8;   % persistencia: 8 ventanas = 8*1.25 ms = 10 ms extra

%% CONTROL: señal sana sin perturbación
lab0 = cls(Vh, W, STRIDE, fs, f);
fprintf('CONTROL sano: falsa operación cruda %.1f%% (%d ventanas)\n', 100*mean(~strcmp(lab0,'A-B-C')), numel(lab0));
fprintf('CONTROL sano con persistencia M=%d: %.1f%%\n', M_PERS, 100*persist_rate(lab0, M_PERS));

%% EXP A2: DC decayente con estimación de DC sobre ciclo anterior
fprintf('\n=== EXP A2: DC decayente, defensa = LS sobre ciclo anterior (80 m.) ===\n');
dc_levels = [0.05 0.10 0.20 0.30]; tau = 0.030; resA = [];
for k = 1:numel(dc_levels)
    dc = dc_levels(k) * Npeak * exp(-th/tau);
    Vdc = Vh; Vdc(:,1) = Vdc(:,1) + dc;
    lab_raw = cls(Vdc, W, STRIDE, fs, f);
    lab_def = cls_dcrem_hist(Vdc, W, STRIDE, fs, f);
    fo_raw  = mean(~strcmp(lab_raw, 'A-B-C'));
    fo_def  = mean(~strcmp(lab_def, 'A-B-C'));
    fo_pers = persist_rate(lab_raw, M_PERS);
    % retención de la falta real con la defensa DC
    dcf = dc_levels(k) * Npeak * exp(-tf_/tau);
    Vfdc = Vf; Vfdc(:,1) = Vfdc(:,1) + dcf;
    lab_f = cls_dcrem_hist(Vfdc, W, STRIDE, fs, f);
    ret = mean(strcmp(lab_f(ceil(numel(lab_f)/4):end), 'A-B')); % ventanas ya dentro de falta
    resA = [resA; dc_levels(k), 100*fo_raw, 100*fo_def, 100*fo_pers, 100*ret];
    fprintf('I_dc=%.2f pu: cruda %5.1f%% | LS-hist %5.1f%% | persistencia %5.1f%% | retención A-B %5.1f%%\n', ...
        dc_levels(k), 100*fo_raw, 100*fo_def, 100*fo_pers, 100*ret);
end
writematrix([["Idc_pu","raw_pct","LShist_pct","persist_pct","retention_pct"]; resA], ...
    fullfile(outdir, 'expA2_dc_defensas.csv'));

%% EXP B2: cap-switching con persistencia
fprintf('\n=== EXP B2: cap-switching, defensa = persistencia M=%d ===\n', M_PERS);
fosc_list = [400 600 900]; amp_list = [0.10 0.20 0.30]; tau_osc = 0.010; t0 = 0.05;
resB = [];
for fo_ = fosc_list
    for am = amp_list
        Vcs = Vh;
        osc = am * Npeak * exp(-max(th - t0, 0)/tau_osc) .* sin(2*pi*fo_*(th - t0)) .* (th >= t0);
        [~, i0_] = min(abs(th - t0)); w0 = Vh(i0_, :) / Npeak;
        for ph = 1:3, Vcs(:,ph) = Vcs(:,ph) + osc .* w0(ph); end
        lab = cls(Vcs, W, STRIDE, fs, f);
        fo_raw  = mean(~strcmp(lab, 'A-B-C'));
        fo_pers = persist_rate(lab, M_PERS);
        resB = [resB; fo_, am, 100*fo_raw, 100*fo_pers];
        fprintf('f=%3d Hz, A=%.2f pu: cruda %5.1f%% | con persistencia %5.1f%%\n', fo_, am, 100*fo_raw, 100*fo_pers);
    end
end
writematrix([["fosc_Hz","amp_pu","raw_pct","persist_pct"]; resB], ...
    fullfile(outdir, 'expB2_cap_persistencia.csv'));

%% Coste de la persistencia: retardo de declaración en la falta real A-B
lab_f = cls(Vf, W, STRIDE, fs, f);
% primera ventana (tras incepción) con M etiquetas iguales consecutivas de falta
starts = 1:STRIDE:(size(Vf,1)-W+1);
tcent = (starts + (W-1)/2)/fs;                 % s desde inicio del segmento (incepción en 0)
det_raw = NaN; det_pers = NaN;
for i = 1:numel(lab_f)
    if strcmp(lab_f{i}, 'A-B') && isnan(det_raw), det_raw = tcent(i); end
    if i >= M_PERS && all(strcmp(lab_f(i-M_PERS+1:i), 'A-B')) && isnan(det_pers)
        det_pers = tcent(i); break
    end
end
fprintf('\nDetección falta real A-B: primera ventana correcta a %.2f ms; con persistencia M=%d a %.2f ms\n', ...
    1e3*det_raw, M_PERS, 1e3*det_pers);

%% EXP D2: secuencia negativa — persistencia no defiende (V2 sostenida); medir umbral del gate
fprintf('\n=== EXP D2: V2 ambiente, umbral del gate de circularidad ===\n');
v2_levels = [0.001 0.002 0.005 0.01 0.02]; resD = [];
for k2 = v2_levels
    v2 = k2*Npeak*[cos(2*pi*f*th), cos(2*pi*f*th + 2*pi/3), cos(2*pi*f*th - 2*pi/3)];
    lab = cls(Vh + v2, W, STRIDE, fs, f);
    fo_raw  = mean(~strcmp(lab, 'A-B-C'));
    fo_pers = persist_rate(lab, M_PERS);
    resD = [resD; k2, 100*fo_raw, 100*fo_pers];
    fprintf('V2=%.1f%%: cruda %5.1f%% | persistencia %5.1f%%\n', 100*k2, 100*fo_raw, 100*fo_pers);
end
writematrix([["V2_pu","raw_pct","persist_pct"]; resD], ...
    fullfile(outdir, 'expD2_v2_umbral.csv'));

fprintf('\nSuite v2 completada.\n');

%% ---------- funciones ----------
function labels = cls(V, W, stride, fs, f)
    n = size(V,1); labels = {};
    for s = 1:stride:(n - W + 1)
        try, ft = ga_fault_classifier(V(s:s+W-1,:), fs, f); catch, ft = 'Error'; end
        labels{end+1} = ft; %#ok<AGROW>
    end
end

function labels = cls_dcrem_hist(V, W, stride, fs, f)
    % DC estimado por LS sobre las 80 muestras ANTERIORES al final de la ventana
    % (un ciclo completo de historia, base bien condicionada), restado por fase.
    n = size(V,1); labels = {};
    H = 80; tH = (0:H-1)'/fs;
    X = [cos(2*pi*f*tH), sin(2*pi*f*tH), ones(H,1)];
    pinvX = pinv(X);
    for s = 1:stride:(n - W + 1)
        e = s + W - 1;
        hs = max(1, e - H + 1);
        Vw = V(s:e,:);
        if e - hs + 1 == H
            hist = V(hs:e,:);
            for ph = 1:3
                coef = pinvX * hist(:,ph);
                Vw(:,ph) = Vw(:,ph) - coef(3);
            end
        end
        try, ft = ga_fault_classifier(Vw, fs, f); catch, ft = 'Error'; end
        labels{end+1} = ft; %#ok<AGROW>
    end
end

function r = persist_rate(labels, M)
    % fracción de ventanas que dispararían tras exigir M etiquetas de falta
    % IDÉNTICAS consecutivas (etiqueta ~= 'A-B-C')
    n = numel(labels); fire = false(1,n);
    for i = M:n
        w = labels(i-M+1:i);
        if all(~strcmp(w, 'A-B-C')) && numel(unique(w)) == 1
            fire(i) = true;
        end
    end
    r = mean(fire);
end
