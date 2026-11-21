function S = r2_sector_scan(V, t, W, STRIDE, fs, f)
% R2_SECTOR_SCAN  Pasa ga_fault_classifier_pure por ventanas deslizantes (T21).
% Devuelve por ventana: etiqueta, angulo, redondez, rama de linea y tiempos de
% inicio y fin (las marcas temporales se toman al FINAL de la ventana).
n = size(V,1); idx = 1:STRIDE:(n-W+1); m = numel(idx);
S.label = cell(1,m); S.theta = zeros(1,m); S.round = zeros(1,m);
S.isline = false(1,m); S.tstart = zeros(1,m); S.tend = zeros(1,m);
for k = 1:m
    s = idx(k);
    [ft, P] = ga_fault_classifier_pure(V(s:s+W-1,:), fs, f);
    S.label{k} = char(ft); S.theta(k) = P.angle_deg; S.round(k) = P.roundness;
    S.isline(k) = P.is_line; S.tstart(k) = t(s); S.tend(k) = t(s+W-1);
end
end
