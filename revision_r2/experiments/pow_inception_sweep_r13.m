%% POW_INCEPTION_SWEEP_R13 — latencia vs ángulo de incepción (point-on-wave, R1.3)
% Falta A-B sintética con el modelo de severidad del paper (v_a-v_b escalado
% por (1-k), que reproduce la tabla: r_min = 1.2247*(1-k)), conmutada en 12
% ángulos de incepción de la fase A (0:30:330). Se mide la latencia hasta la
% primera ventana correcta y hasta la declaración asegurada (M=10 idénticas),
% ambas respecto al instante de conmutación (fin de ventana).
clear; clc;
addpath(fullfile('..', 'src', 'fault_classifiers'));
fs = 4000; f = 50; W = 20; STRIDE = 5; M = 10;
outdir = fullfile('..', '_revision_TPWRD', 'robustness');

T = 1.0; t = (0:1/fs:T-1/fs)';
ea = cos(2*pi*f*t); eb = cos(2*pi*f*t - 2*pi/3); ec = cos(2*pi*f*t + 2*pi/3);

res = [];
for k = [0.3 0.6]
    for phi = 0:30:330
        % instante de conmutación: fase A en ángulo phi, tras >=0.5 s de régimen
        t0 = 0.5 + (phi/360)/f;
        on = t >= t0;
        d = (ea - eb)/2;
        va = ea - k*d.*on; vb = eb + k*d.*on; vc = ec;
        V = [va, vb, vc];
        % clasificar desde 1 ciclo antes de t0
        i0 = find(t >= t0, 1) - round(fs/f);
        lab = {}; tend = [];
        for s = i0:STRIDE:(numel(t)-W+1)
            try, ft = ga_fault_classifier(V(s:s+W-1,:), fs, f); catch, ft='Error'; end
            lab{end+1} = ft; tend(end+1) = t(s+W-1);
        end
        post = tend > t0; d1 = NaN; dM = NaN;
        for i = find(post)
            if strcmp(lab{i}, 'A-B') && isnan(d1), d1 = tend(i) - t0; end
            if i >= M && all(strcmp(lab(i-M+1:i), 'A-B')) && isnan(dM), dM = tend(i) - t0; break, end
        end
        res = [res; k, phi, 1e3*d1, 1e3*dM]; %#ok<AGROW>
    end
    rk = res(res(:,1)==k, :);
    fprintf('k=%.1f: 1a correcta [min %.2f, mediana %.2f, max %.2f] ms | asegurada [min %.2f, mediana %.2f, max %.2f] ms | fallos: %d/12\n', ...
        k, min(rk(:,3)), median(rk(:,3),'omitnan'), max(rk(:,3)), ...
        min(rk(:,4)), median(rk(:,4),'omitnan'), max(rk(:,4)), sum(isnan(rk(:,4))));
end
writematrix([["k","phi_deg","first_ms","secured_ms"]; res], ...
    fullfile(outdir, 'r13_pow_sweep.csv'));
