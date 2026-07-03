%% SEVERITY_ESTIMATOR_R17 — estimador de severidad no lineal (reviewer R1.7)
% Genera la relación teórica severidad k -> semieje menor r para faltas L-G
% (modelo de severidad del paper: v_a=(1-k)e_a, v_b=e_b, v_c=e_c) con el
% clasificador original a ventana de 1 ciclo, y compara estimadores inversos:
% lineal (el del paper), cuadrático y lineal a tramos (2 segmentos).
clear; clc;
addpath(fullfile('..', 'classifiers'));
fs = 4000; f = 50; W = 80;
t = (0:W-1)'/fs;
ea = cos(2*pi*f*t); eb = cos(2*pi*f*t - 2*pi/3); ec = cos(2*pi*f*t + 2*pi/3);

kk = (0.02:0.01:0.95)';
rr = zeros(size(kk));
for i = 1:numel(kk)
    V = [(1-kk(i))*ea, eb, ec];
    [~, p] = ga_fault_classifier(V, fs, f);
    rr(i) = p.semi_minor;
end
% comprobación contra la tabla del paper (k=0.1 -> 1.1446, k=0.9 -> 0.7141)
i1 = find(abs(kk-0.1)<1e-9); i9 = find(abs(kk-0.9)<1e-9);
fprintf('Verificación tabla: r(0.1)=%.4f (paper 1.1446), r(0.9)=%.4f (paper 0.7141)\n', rr(i1), rr(i9));

% --- estimadores inversos k_hat(r) ajustados sobre la curva densa ---
% 1) lineal
c1 = polyfit(rr, kk, 1);  k_lin = polyval(c1, rr);
% 2) cuadrático
c2 = polyfit(rr, kk, 2);  k_qua = polyval(c2, rr);
% 3) lineal a tramos (nudo en la mediana de r)
rknot = median(rr);
lo = rr <= rknot; hi = ~lo;
cl = polyfit(rr(lo), kk(lo), 1); ch = polyfit(rr(hi), kk(hi), 1);
k_pw = zeros(size(kk)); k_pw(lo) = polyval(cl, rr(lo)); k_pw(hi) = polyval(ch, rr(hi));

mae = @(kh) 100*mean(abs(kh - kk));
fprintf('MAE severidad L-G: lineal %.2f%% | cuadrático %.2f%% | lineal a tramos %.2f%%\n', ...
    mae(k_lin), mae(k_qua), mae(k_pw));
fprintf('Coefs cuadrático: k = %.4f r^2 + %.4f r + %.4f\n', c2(1), c2(2), c2(3));

out = table(kk, rr, k_lin, k_qua, k_pw, 'VariableNames', {'k','r_minor','k_lineal','k_cuadratico','k_tramos'});
outdir = fullfile('..', 'results', 'robustness'); if ~exist(outdir, 'dir'), mkdir(outdir); end
writetable(out, fullfile(outdir, 'r17_severidad_LG.csv'));
fprintf('Guardado r17_severidad_LG.csv\n');
