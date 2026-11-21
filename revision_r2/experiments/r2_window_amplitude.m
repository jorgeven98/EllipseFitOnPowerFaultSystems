function [a, b] = r2_window_amplitude(V, i_start, W, fs, f)
% R2_WINDOW_AMPLITUDE  Semiejes de la trayectoria de cada ventana por ajuste sinusoidal.
%
% Para cada ventana ajusta por minimos cuadrados v(t) = c cos(wt) + s sin(wt)
% en las tres fases (sin termino continuo: con un cuarto de ciclo no se puede
% separar) y devuelve los valores singulares de la matriz [c s] (3 x 2), que
% son los semiejes mayor y menor de la elipse 3D a la frecuencia f, en las
% unidades de V. En regimen equilibrado de pico Vp: a = b = sqrt(3/2) Vp.
% Sirve para medir la longitud real de la linea en la rama degenerada, donde
% el clasificador fuerza semi_major = sqrt(3)/sqrt(2) sin medirlo (T41).
m = numel(i_start); a = nan(m, 1); b = nan(m, 1);
tw = (0:W-1)'/fs; X = [cos(2*pi*f*tw), sin(2*pi*f*tw)]; Xp = pinv(X);
for j = 1:m
    if ~isfinite(i_start(j)), continue, end
    C = Xp * V(i_start(j):i_start(j)+W-1, :);      % 2 x 3
    s = svd(C);
    a(j) = s(1); b(j) = s(2);
end
end
