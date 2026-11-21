function [k_ini, k_end, info] = r2_lab_inception(V, t, fs, f)
% R2_LAB_INCEPTION  Estima el inicio y el final de la falta en un registro de laboratorio.
% Mismo criterio que experiments/lab_latency_r3m1.m (no lo modifica): region de falta por
% desequilibrio del valor eficaz por ciclo y, dentro de ella, primera muestra en que la
% tension se aparta mas del 5 % del pico de la extrapolacion sinusoidal de prefalta,
% sostenido 5 muestras. El final se estima como el ultimo ciclo con desequilibrio.
%   V [N x 3] tensiones (cualquier unidad), t [N x 1], fs, f.
%   k_ini, k_end: indices de muestra (NaN si no hay falta detectable).
% El emulador no registra el mando: son estimaciones a partir de la tension.
cyc = round(fs/f); nw = floor(numel(t)/cyc);
dev = zeros(nw,1); lev = zeros(nw,1);
for i = 1:nw
    seg = (i-1)*cyc+1 : i*cyc;
    r3 = sqrt(mean(V(seg,:).^2, 1));
    dev(i) = max(r3) - min(r3); lev(i) = mean(r3);
end
base = median(dev); lev0 = median(lev);
% desequilibrio (faltas asimetricas) o hundimiento comun (trifasicas)
anomv = dev > 5*base + 0.02*lev0 | lev < 0.8*lev0;
anom = find(anomv, 1);
k_ini = NaN; k_end = NaN; info = struct('base_dev', base, 'level', lev0);
if isempty(anom) || anom < 5, return; end
i_region = (anom-2)*cyc;
idx_fit = i_region-2*cyc : i_region-1;
tt = t(idx_fit); X = [cos(2*pi*f*tt), sin(2*pi*f*tt), ones(numel(tt),1)];
coefs = X \ V(idx_fit,:);
idx_test = i_region : min(i_region+3*cyc, numel(t));
Xt = [cos(2*pi*f*t(idx_test)), sin(2*pi*f*t(idx_test)), ones(numel(idx_test),1)];
resid = abs(V(idx_test,:) - Xt*coefs);
Vpk = max(abs(V(idx_fit,:)), [], 'all');
hit = all(movmin(max(resid,[],2), [0 4]) > 0.05*Vpk, 2);
ii = find(hit, 1);
if isempty(ii), return; end
k_ini = idx_test(ii);
last = find(anomv, 1, 'last');
k_end = min(last*cyc, numel(t));
info.n_cycles_anom = nnz(anomv);
end
