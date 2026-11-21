function [u0, u2, aux] = r2_ground_rf_seqlaw(tipo, Z0, Z2, Rf, Rg)
% R2_GROUND_RF_SEQLAW  Leyes de secuencia en el punto de falta para faltas a
% tierra resistivas (T23b). Devuelve u0 = V0/V1 y u2 = V2/V1 referidos a la
% fase a (secuencia abc), para:
%   'A-G','B-G','C-G'      falta monofasica con resistencia Rf (Rg se suma a Rf)
%   'B-C-G','C-A-G','A-B-G' falta bifasica a tierra con Rf en cada fase y Rg
%                           comun entre el punto estrella de la falta y tierra
% Z0, Z2: impedancias Thevenin homopolar e inversa vistas desde el punto de
% falta (la red directa no interviene en los cocientes). Rf, Rg escalares o
% matrices del mismo tamano.
%
% L-G (fase singular = fase en falta):
%   u2 = -Z2/(Z0+Z2+3R),  u0 = -Z0/(Z0+Z2+3R),  R = Rf+Rg   => V0/V2 = Z0/Z2
% L-L-G (fase singular = fase sana), Z2' = Z2+Rf, Z0' = Z0+Rf+3Rg:
%   D  = Rf (Z2'+Z0') + Z2' Z0'
%   u2 = Z2 Z0'/D,  u0 = Z0 Z2'/D          => V0/V2 = (Z0/Z2)(Z2+Rf)/(Z0+Rf+3Rg)
% Cambio de fase singular k (a->0, b->1, c->2): u2 -> u2 a^k, u0 -> u0 a^(2k)
%   (V1 = a^2 V1b, V2 = a V2b para k = 1).
a = exp(2j*pi/3);
tp = upper(strrep(char(tipo), '-', ''));
switch tp
    case 'AG',  fam = 1; k = 0;
    case 'BG',  fam = 1; k = 1;
    case 'CG',  fam = 1; k = 2;
    case 'BCG', fam = 2; k = 0;
    case {'CAG','ACG'}, fam = 2; k = 1;
    case 'ABG', fam = 2; k = 2;
    otherwise, error('tipo no admitido: %s', tipo);
end
if fam == 1
    D = Z0 + Z2 + 3*(Rf + Rg);
    u2s = -Z2 ./ D; u0s = -Z0 ./ D;
else
    Z2p = Z2 + Rf; Z0p = Z0 + Rf + 3*Rg;
    D = Rf.*(Z2p + Z0p) + Z2p.*Z0p;
    u2s = Z2.*Z0p ./ D; u0s = Z0.*Z2p ./ D;
end
% de la referencia en la fase singular a la referencia en la fase a
u2 = u2s * a^k; u0 = u0s * a^(2*k);
aux.fam = fam; aux.k = k; aux.u2s = u2s; aux.u0s = u0s;
% inclinacion compleja del plano: n ~ k + rho (cos(d), cos(d-120), cos(d+120))
aux.tau = 2*(conj(u2).*u0 - conj(u0)) ./ (1 - abs(u2).^2);
end
