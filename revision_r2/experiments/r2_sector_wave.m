function V = r2_sector_wave(t, V1, V2, f, V0)
% R2_SECTOR_WAVE  Tensiones de fase a partir de fasores de secuencia (T21).
% Convenio: v_k(t) = Re{V_k e^{j 2 pi f t}}, secuencia abc,
%   Va = V0 + V1 + V2,  Vb = V0 + a^2 V1 + a V2,  Vc = V0 + a V1 + a^2 V2.
% V1, V2, V0 pueden ser escalares o vectores por muestra (cambio prefalta/falta).
if nargin < 5, V0 = 0; end
a = exp(1j*2*pi/3);
t = t(:); e = exp(1j*2*pi*f*t);
V1 = V1(:); V2 = V2(:); V0 = V0(:);
V = real([(V0+V1+V2).*e, (V0+a^2*V1+a*V2).*e, (V0+a*V1+a^2*V2).*e]);
end
