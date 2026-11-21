function R = r2_sector_seqsolve(P)
% R2_SECTOR_SEQSOLVE  Redes de secuencia para una falta B-C con Rf (T21).
% Sistema:  convertidor -- M --[Z_MF]-- F --[Z_FG]-- red (Eg tras Z_FG)
%   M  punto de medida (bornes del convertidor); F punto de falta.
%   Z_FG incluye el tramo de linea F-G y la impedancia del equivalente de red.
% Entradas (struct P, en pu o en ohmios coherentes):
%   Eg     f.e.m. de la red (fasor)
%   Z_FG   impedancia F-red (se supone Z1 = Z2 en los elementos pasivos)
%   Z_MF   impedancia M-F
%   Rf     impedancia TOTAL de falta entre las fases b y c (puede ser compleja)
%   conv   'gfl' | 'gfm'
%   Z2c    impedancia de secuencia inversa del convertidor vista desde M
%          (Inf: corriente inversa suprimida)
%   gfl:   Imax (modulo de I1), gamma (rad, angulo de I1 respecto de V1 en M;
%          gamma = 0 corriente activa, gamma = -pi/2 reactiva inyectada)
%   gfm:   Ec (f.e.m. interna, fasor), Zc (impedancia interna en secuencia directa)
% Salidas: V1F, V2F, V1M, V2M, qF, qM, Z2th, y los factores de transferencia
%   k1 = V1M/V1F, k2 = V2M/V2F, de modo que qM = qF k2/k1.
if isinf(P.Z2c), Z2th = P.Z_FG; k2 = 1;
else, Z2th = 1/(1/P.Z_FG + 1/(P.Z_MF + P.Z2c)); k2 = P.Z2c/(P.Z2c + P.Z_MF); end
switch lower(P.conv)
    case 'gfm'
        Zb = P.Z_MF + P.Zc;
        V1th = (P.Eg/P.Z_FG + P.Ec/Zb)/(1/P.Z_FG + 1/Zb);
        Z1th = 1/(1/P.Z_FG + 1/Zb);
        I1F = V1th/(Z1th + Z2th + P.Rf);
        V1F = V1th - Z1th*I1F;
        V1M = V1F + P.Z_MF*(P.Ec - V1F)/Zb;
        I1c = (P.Ec - V1F)/Zb;
    case 'gfl'
        % fuente de corriente de secuencia directa; su angulo se refiere a V1M
        I1c = P.Imax*exp(1j*P.gamma);  V1M = 1;
        for it = 1:200
            V1th = P.Eg + P.Z_FG*I1c;  Z1th = P.Z_FG;
            I1F = V1th/(Z1th + Z2th + P.Rf);
            V1F = V1th - Z1th*I1F;
            V1Mn = V1F + P.Z_MF*I1c;
            I1n = P.Imax*exp(1j*(angle(V1Mn) + P.gamma));
            if abs(I1n - I1c) < 1e-13, V1M = V1Mn; I1c = I1n; break; end
            I1c = 0.5*I1c + 0.5*I1n; V1M = V1Mn;
        end
    otherwise, error('conv desconocido');
end
V2F = Z2th*I1F;  V2M = k2*V2F;
R.V1F = V1F; R.V2F = V2F; R.V1M = V1M; R.V2M = V2M; R.I1c = I1c;
R.qF = V2F/V1F; R.qM = V2M/V1M; R.Z2th = Z2th; R.k1 = V1M/V1F; R.k2 = k2;
end
