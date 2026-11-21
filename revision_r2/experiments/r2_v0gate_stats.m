function st = r2_v0gate_stats(Vw, fs, f, vpk)
% R2_V0GATE_STATS  Estadisticos candidatos para la puerta de secuencia cero
% (T51 y T52, ronda 2 TPWRD). No clasifica: solo mide la ventana.
%
%   st = r2_v0gate_stats(Vw, fs, f, vpk)
%
% Vw  [N x 3] ventana de tensiones fase-tierra (cualquier unidad).
% vpk pico de fase nominal en las mismas unidades (solo para los estadisticos
%     '_nom'); si se omite o es NaN esos campos salen NaN.
%
% Notacion: k^ = (1,1,1)/sqrt(3); p_i = x_i . k^ = sqrt(3) v0_i, con
% v0 = (va + vb + vc)/3. Todos los estadisticos 'rms_*' y 'ph_*' estan
% escalados para valer V0/Vpico en un sistema sano equilibrado con una
% secuencia cero de amplitud V0 (valor eficaz de ciclo completo), de modo que
% su umbral se lee como «V0 en por unidad».
%
% Campos:
%   s3       G0 (puerta publicada): max |x^_i . k^| en las tres muestras idx1,
%            idx2, idx3 de ga_fault_classifier_pure (normaliza muestra a muestra)
%   tilt     G1: seno del angulo entre la normal del plano ajustado por SVD
%            (todas las muestras, sin centrar) y k^
%   tilt2    G1b: lo mismo con el bivector de dos muestras x1 ^ x2 que ya
%            calcula el clasificador (angulo entre bivector y plano de Kirchhoff)
%   rms_win  G5: ||X k^||_F / ||X||_F = sqrt(3) rms(v0) / rms(||x||)
%            (misma magnitud que G0 pero agregada ANTES de normalizar)
%   rms_nom  G5n: sqrt(2) rms(v0) / vpk (escala nominal fija)
%   rms_max  G5m: sqrt(3) rms(v0) / max_i ||x_i|| (escala: modulo maximo de la
%            ventana, cota inferior del semieje mayor)
%   ph_win   G5p: amplitud del fasor fundamental de v0 por minimos cuadrados
%            (cos, sin; sin termino continuo) sobre rms(||x||)
%   ph_nom   G5pn: la misma amplitud sobre vpk
%   sv       [1x3] valores singulares de la ventana
%   n_svd    [1x3] normal del plano por SVD (signo: n . k^ >= 0)
%   d_svd    [1x3] direccion principal (primer vector singular derecho)
%   gnd_svd  etiqueta de tierra que daria classify_by_bivector con n_svd
%   gnd_line etiqueta de tierra a partir de la direccion principal cuando la
%            trayectoria es una recta (fase sana dominante -> L-L-G de las otras
%            dos); '' si la direccion no senala una sola fase
if nargin < 4 || isempty(vpk), vpk = NaN; end
N = size(Vw, 1);
q = round(fs/(4*f)); if q < 1, q = 1; end
if N < q, i1 = 1; i2 = N; i3 = 1;
else, c = round(N/2); i1 = max(1, c - q); i2 = min(N, i1 + q); i3 = round(abs(i1 - i2)/2);
end
k = [1 1 1]/sqrt(3);
xs = Vw([i1 i2 i3], :); xh = xs ./ sqrt(sum(xs.^2, 2));
st.s3 = max(abs(xh*k'));
n2 = cross(xh(1, :), xh(2, :)); bn = norm(n2);
st.tilt2 = sqrt(max(0, 1 - (n2*k'/bn)^2));
[~, Sv, Vr] = svd(Vw, 'econ'); s = diag(Sv)';
n = Vr(:, 3)'; if n*k' < 0, n = -n; end
st.tilt = sqrt(max(0, 1 - (n*k')^2));
st.sv = s; st.n_svd = n; st.d_svd = Vr(:, 1)';
p = Vw*k';
rp = sqrt(mean(p.^2)); nx = sqrt(sum(Vw.^2, 2)); rx = sqrt(mean(nx.^2));
st.rms_win = rp/rx;
st.rms_nom = rp/(sqrt(1.5)*vpk);
st.rms_max = rp/max(nx);
tt = (0:N-1)'/fs; X = [cos(2*pi*f*tt), sin(2*pi*f*tt)];
cf = X\p; A = hypot(cf(1), cf(2));
st.ph_win = A/(sqrt(2)*rx);
st.ph_nom = A/(sqrt(3)*vpk);
st.gnd_svd = by_bivector_([n(3), n(1), n(2)], 0.3);
d = abs(Vr(:, 1)'); [ds, id] = sort(d, 'descend');
if ds(2) < 0.5*ds(1)
    lab = {'B-C-G', 'C-A-G', 'A-B-G'}; st.gnd_line = lab{id(1)};
else
    st.gnd_line = '';
end
end

function fault_type = by_bivector_(biv, dif_thr)
% Copia literal de classify_by_bivector de ga_fault_classifier_pure
[biv_sort, idx] = sort(abs(biv), "descend");
dif1 = abs(biv_sort(1) - biv_sort(2));
dif2 = abs(biv_sort(2) - biv_sort(3));
if dif1 < dif2 || dif2 > dif_thr
    switch idx(3)
        case 1, fault_type = 'A-B-G';
        case 2, fault_type = 'B-C-G';
        case 3, fault_type = 'C-A-G';
    end
else
    switch idx(1)
        case 1, fault_type = 'C-G';
        case 2, fault_type = 'A-G';
        case 3, fault_type = 'B-G';
    end
end
end
