%% ROBUSTNESS_SUITE_TPWRD
% Suite de robustez no gaussiana para la revisión TPWRD-00232-2026 (P3).
% Cuatro experimentos con el clasificador ORIGINAL (ga_fault_classifier, GA-FuL):
%   A) DC decayente superpuesto (spoofing + defensa por LS DC-removal)   [R3.C1]
%   B) Transitorio de conmutación de banco de condensadores 400-900 Hz   [R3.C1]
%   C) SLG en sistema NO aterrado (desplazamiento de neutro, v0/i0)      [R1.4]
%   D) Inyección de secuencia negativa (grid codes VDE-AR-N 4120/IEEE2800) [R3.C2]
%
% Señales base: datos reales de simulación (data/simulation_curated, Pex/DG1)
% para A, B y D; generación analítica en forma cerrada para C.
% Salidas: _revision_TPWRD/robustness/*.csv + resumen en consola.
%
% Autor: Jorge (con Claude) - Revisión TPWRD, jul-2026

clear; clc;
addpath(fullfile('..', 'classifiers'));
addpath(fullfile('..', 'utils'));

fs = 4000; f = 50; Tsamp = 1/fs;
outdir = fullfile('..', 'results', 'robustness');
if ~exist(outdir, 'dir'), mkdir(outdir); end

% --- Señal base sana (real, Pex pre-falta) ---
base = load(fullfile('..', 'data', 'DataGridRL_AB.mat'));
[Vpex, tpex] = extract_converter_data(base, 'Pex');
[Vdg1, ~   ] = extract_converter_data(base, 'DG1');
seg = tpex >= 2.4 & tpex < 2.8;                 % 0.4 s sanos, régimen establecido
Vh_pex = Vpex(seg, :); Vh_dg1 = Vdg1(seg, :);
th = tpex(seg); th = th - th(1);
Npeak = max(abs(Vh_pex(:)));                    % amplitud nominal de pico (~1 pu)
fprintf('Señal sana: %d muestras, pico=%.3f\n', size(Vh_pex,1), Npeak);

% --- Ventanas de falta reales (para retención de clasificación) ---
segF = tpex >= 3.02 & tpex < 3.14;              % dentro de la falta A-B
Vf_pex = Vpex(segF, :); tf_ = tpex(segF); tf_ = tf_ - tf_(1);

% --- utilidades ---
W = 20; STRIDE = 5;                             % 1/4 ciclo, solape 75%
classify_all = @(V) classify_windows(V, W, STRIDE, fs, f);

%% ========================================================================
%  EXPERIMENTO A: DC DECAYENTE
%  ========================================================================
fprintf('\n=== EXP A: DC decayente (tau=30 ms, fase A) ===\n');
dc_levels = [0.05 0.10 0.20 0.30];
tau = 0.030;
resA = [];
for k = 1:numel(dc_levels)
    dc = dc_levels(k) * Npeak * exp(-th/tau);
    Vdc = Vh_pex; Vdc(:,1) = Vdc(:,1) + dc;     % DC en fase A
    % sin defensa
    lab_raw = classify_all(Vdc);
    fo_raw = mean(~strcmp(lab_raw, 'A-B-C'));
    % con defensa: LS DC-removal por ventana (base [cos,sin,1] a f nominal)
    lab_def = classify_windows_dcrem(Vdc, W, STRIDE, fs, f);
    fo_def = mean(~strcmp(lab_def, 'A-B-C'));
    % retención en falta real A-B con DC
    dcf = dc_levels(k) * Npeak * exp(-tf_/tau);
    Vfdc = Vf_pex; Vfdc(:,1) = Vfdc(:,1) + dcf;
    lab_f = classify_windows_dcrem(Vfdc, W, STRIDE, fs, f);
    ret = mean(strcmp(lab_f, 'A-B'));
    resA = [resA; dc_levels(k), 100*fo_raw, 100*fo_def, 100*ret]; %#ok<*AGROW>
    fprintf('I_dc=%.2f pu: falsa op. sin defensa %5.1f%%, con LS-DC %5.1f%% | retención falta A-B %5.1f%%\n', ...
        dc_levels(k), 100*fo_raw, 100*fo_def, 100*ret);
end
writematrix([["Idc_pu","falseop_raw_pct","falseop_LSdc_pct","fault_retention_pct"]; resA], ...
    fullfile(outdir, 'expA_dc_decayente.csv'));

%% ========================================================================
%  EXPERIMENTO B: CONMUTACION DE BANCO DE CONDENSADORES
%  ========================================================================
fprintf('\n=== EXP B: cap-switching (oscilación amortiguada, tau=10 ms) ===\n');
fosc_list = [400 600 900]; amp_list = [0.10 0.20 0.30]; tau_osc = 0.010;
t0 = 0.05;                                       % instante de conmutación
resB = [];
for fo_ = fosc_list
    for am = amp_list
        Vcs = Vh_pex;
        osc = am * Npeak * exp(-max(th - t0, 0)/tau_osc) .* sin(2*pi*fo_*(th - t0)) .* (th >= t0);
        % amplitud por fase proporcional a la tensión instantánea en t0 (física del cierre)
        [~, i0_] = min(abs(th - t0));
        w0 = Vh_pex(i0_, :) / Npeak;
        for ph = 1:3, Vcs(:,ph) = Vcs(:,ph) + osc .* w0(ph); end
        lab = classify_all(Vcs);
        % solo ventanas que tocan el transitorio (t0 a t0+3*tau)
        idx_tr = window_centers(size(Vcs,1), W, STRIDE)/fs;
        touch = idx_tr >= t0 - W/fs & idx_tr <= t0 + 3*tau_osc;
        fo_rate = mean(~strcmp(lab(touch), 'A-B-C'));
        resB = [resB; fo_, am, 100*fo_rate, sum(touch)];
        fprintf('f=%3d Hz, A=%.2f pu: falsa operación %5.1f%% (%d ventanas afectadas)\n', ...
            fo_, am, 100*fo_rate, sum(touch));
    end
end
writematrix([["fosc_Hz","amp_pu","falseop_pct","n_windows"]; resB], ...
    fullfile(outdir, 'expB_cap_switching.csv'));

%% ========================================================================
%  EXPERIMENTO C: SLG EN SISTEMA NO ATERRADO (analítico)
%  ========================================================================
fprintf('\n=== EXP C: SLG no aterrado, desplazamiento de neutro V_N ===\n');
% Sistema aislado: falta A-G con severidad k = |V_N|/E (k=1 <=> falta franca).
% Tensiones fase-tierra: v_ph = e_ph + v_N, con v_N = -k*e_a (en fase con e_a).
tC = (0:Tsamp:0.4-Tsamp)';
ea = cos(2*pi*f*tC); eb = cos(2*pi*f*tC - 2*pi/3); ec = cos(2*pi*f*tC + 2*pi/3);
k_list = [0.25 0.5 0.75 1.0];
C0 = 1e-6;                                       % capacidad homopolar por fase (típica distribución)
resC = [];
for k = k_list
    vN = -k * ea;
    Vslg = [ea + vN, eb + vN, ec + vN];
    lab = classify_all(Vslg);
    accAG = mean(strcmp(lab, 'A-G'));
    v0 = mean(Vslg, 2);                          % tensión homopolar instantánea
    i0max = 3 * 2*pi*f * C0 * max(abs(v0));      % corriente residual capacitiva (A/pu·Ohm base)
    resC = [resC; k, 100*accAG, max(abs(v0)), i0max];
    fprintf('k=%.2f: clasif. A-G %5.1f%% | v0_max=%.3f pu | i0 capacitiva ~%.2e pu\n', ...
        k, 100*accAG, max(abs(v0)), i0max);
end
writematrix([["severidad_k","acc_AG_pct","v0max_pu","i0_cap_pu"]; resC], ...
    fullfile(outdir, 'expC_slg_no_aterrado.csv'));
% CSV para figura v0/i0 del paper (caso k=1)
vN = -1.0 * ea; Vslg = [ea + vN, eb + vN, ec + vN];
v0 = mean(Vslg, 2); i0 = 3*C0*gradient(v0, Tsamp);   % i0 = 3*C0*dv0/dt
Tfig = table(tC, Vslg(:,1), Vslg(:,2), Vslg(:,3), v0, i0, ...
    'VariableNames', {'t','va','vb','vc','v0','i0'});
writetable(Tfig(1:2:800,:), fullfile(outdir, 'expC_slg_v0i0_waveforms.csv'));

%% ========================================================================
%  EXPERIMENTO D: INYECCION DE SECUENCIA NEGATIVA
%  ========================================================================
fprintf('\n=== EXP D: secuencia negativa V2 superpuesta ===\n');
v2_levels = [0.005 0.01 0.02 0.05 0.10];
resD = [];
for k2 = v2_levels
    % V2: sistema de secuencia inversa (a,c,b) en pu de la amplitud nominal
    v2a = k2*Npeak*cos(2*pi*f*th); v2b = k2*Npeak*cos(2*pi*f*th + 2*pi/3); v2c = k2*Npeak*cos(2*pi*f*th - 2*pi/3);
    % (a) ambiente: sano + V2 (¿falsa operación?)
    Vamb = Vh_pex + [v2a, v2b, v2c];
    lab_a = classify_all(Vamb);
    fo_amb = mean(~strcmp(lab_a, 'A-B-C'));
    % (b) falta real A-B + V2 (¿retención?)
    v2af = k2*Npeak*cos(2*pi*f*tf_); v2bf = k2*Npeak*cos(2*pi*f*tf_ + 2*pi/3); v2cf = k2*Npeak*cos(2*pi*f*tf_ - 2*pi/3);
    Vfv2 = Vf_pex + [v2af, v2bf, v2cf];
    lab_b = classify_all(Vfv2);
    ret = mean(strcmp(lab_b, 'A-B'));
    resD = [resD; k2, 100*fo_amb, 100*ret];
    fprintf('V2=%.1f%%: falsa op. ambiente %5.1f%% | retención falta A-B %5.1f%%\n', ...
        100*k2, 100*fo_amb, 100*ret);
end
writematrix([["V2_pu","falseop_ambiente_pct","fault_retention_pct"]; resD], ...
    fullfile(outdir, 'expD_secuencia_negativa.csv'));

fprintf('\nSuite completada. CSVs en %s\n', outdir);

%% ========================================================================
%  FUNCIONES LOCALES
%  ========================================================================
function labels = classify_windows(V, W, stride, fs, f)
    n = size(V,1); labels = {};
    for s = 1:stride:(n - W + 1)
        try
            ft = ga_fault_classifier(V(s:s+W-1,:), fs, f);
        catch
            ft = 'Error';
        end
        labels{end+1} = ft;
    end
end

function labels = classify_windows_dcrem(V, W, stride, fs, f)
    % Igual que classify_windows pero con estimación LS de DC por ventana:
    % por fase, regresión sobre [cos(2*pi*f*t), sin(2*pi*f*t), 1] y resta del término constante.
    n = size(V,1); labels = {};
    t = (0:W-1)'/fs;
    X = [cos(2*pi*f*t), sin(2*pi*f*t), ones(W,1)];
    pinvX = pinv(X);
    for s = 1:stride:(n - W + 1)
        Vw = V(s:s+W-1,:);
        for ph = 1:3
            coef = pinvX * Vw(:,ph);
            Vw(:,ph) = Vw(:,ph) - coef(3);
        end
        try
            ft = ga_fault_classifier(Vw, fs, f);
        catch
            ft = 'Error';
        end
        labels{end+1} = ft;
    end
end

function centers = window_centers(n, W, stride)
    starts = 1:stride:(n - W + 1);
    centers = starts + (W-1)/2;
end
