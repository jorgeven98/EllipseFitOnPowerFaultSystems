%% ROBUSTNESS_SENSOR_ERRORS_R12 — errores de ganancia/fase de sensor (R1.2)
% Peor caso realista: error en UNA fase (el error común a las 3 no cambia la
% geometría). Un error de ganancia eps en la fase A genera componentes V0 y V2
% de ~eps/3 cada una — es la vía por la que el error de sensor ataca al método.
% Métricas: falsa operación en sano (cruda y con gate M=10 de etiquetas
% idénticas) y retención de la falta real A-B con el error aplicado.
clear; clc;
addpath(fullfile('..', 'src', 'fault_classifiers'));
addpath(fullfile('..', 'src', 'utils'));
fs = 4000; f = 50; W = 20; STRIDE = 5; M = 10;
outdir = fullfile('..', '_revision_TPWRD', 'robustness');

base = load(fullfile('..', 'data', 'simulation_curated', 'DataGridRL_AB.mat'));
[Vpex, tpex] = extract_converter_data(base, 'Pex');
seg = tpex >= 2.4 & tpex < 2.8; Vh = Vpex(seg,:); th = tpex(seg); th = th - th(1);
segF = tpex >= 3.00 & tpex < 3.14; Vf = Vpex(segF,:);

fprintf('=== Error de GANANCIA en fase A ===\n');
resG = [];
for eps_ = [0.005 0.01 0.03]
    Vg = Vh; Vg(:,1) = Vg(:,1) * (1 + eps_);
    lab = cls(Vg, W, STRIDE, fs, f);
    fo_raw = mean(~strcmp(lab, 'A-B-C'));
    fo_sec = prate_strict(lab, M);
    Vgf = Vf; Vgf(:,1) = Vgf(:,1) * (1 + eps_);
    labf = cls(Vgf, W, STRIDE, fs, f);
    ret = mean(strcmp(labf, 'A-B'));
    resG = [resG; 100*eps_, 100*fo_raw, 100*fo_sec, 100*ret]; %#ok<*AGROW>
    fprintf('gain +%.1f%%: falsa op cruda %5.1f%% | asegurada %5.1f%% | retención A-B %5.1f%%\n', ...
        100*eps_, 100*fo_raw, 100*fo_sec, 100*ret);
end
writematrix([["gain_pct","raw_pct","secured_pct","retention_pct"]; resG], ...
    fullfile(outdir, 'r12_gain_error.csv'));

fprintf('\n=== Error de FASE en fase A (retardo temporal interpolado) ===\n');
resP = [];
for deg = [0.5 1.0]
    dt = deg/360/f;                          % s de retardo
    tq = th - dt;                            % muestrear la fase A retardada
    Vp = Vh;
    Vp(:,1) = interp1(th, Vh(:,1), tq, 'spline', 'extrap');
    lab = cls(Vp, W, STRIDE, fs, f);
    fo_raw = mean(~strcmp(lab, 'A-B-C'));
    fo_sec = prate_strict(lab, M);
    tf_ = (0:size(Vf,1)-1)'/fs;
    Vpf = Vf; Vpf(:,1) = interp1(tf_, Vf(:,1), tf_ - dt, 'spline', 'extrap');
    labf = cls(Vpf, W, STRIDE, fs, f);
    ret = mean(strcmp(labf, 'A-B'));
    resP = [resP; deg, 100*fo_raw, 100*fo_sec, 100*ret];
    fprintf('fase %.1f deg: falsa op cruda %5.1f%% | asegurada %5.1f%% | retención A-B %5.1f%%\n', ...
        deg, 100*fo_raw, 100*fo_sec, 100*ret);
end
writematrix([["phase_deg","raw_pct","secured_pct","retention_pct"]; resP], ...
    fullfile(outdir, 'r12_phase_error.csv'));

function labels = cls(V, W, stride, fs, f)
    n = size(V,1); labels = {};
    for s = 1:stride:(n - W + 1)
        try, ft = ga_fault_classifier(V(s:s+W-1,:), fs, f); catch, ft='Error'; end
        labels{end+1} = ft;
    end
end
function r = prate_strict(labels, M)
    n = numel(labels); fire = false(1,n);
    for i = M:n
        w = labels(i-M+1:i);
        if all(~strcmp(w,'A-B-C')) && numel(unique(w))==1, fire(i)=true; end
    end
    r = mean(fire);
end
