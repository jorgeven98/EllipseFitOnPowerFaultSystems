function [V, t, u0, u2, V1] = r2_ground_rf_synth(tipo, Z0, Z1, Z2, Rf, Rg, E, dur, fs, f, toff)
% R2_GROUND_RF_SYNTH  Onda sintetica en regimen permanente de una falta a
% tierra resistiva en el punto de falta, desde las redes de secuencia (T23b).
% E: f.e.m. Thevenin de pico (fase a, fase 0). toff: desplazamiento del origen
% de tiempos en muestras (0.37 por defecto, como T21, para no muestrear el
% paso por cero comun de una falta franca).
if nargin < 11, toff = 0.37; end
[u0, u2, ax] = r2_ground_rf_seqlaw(tipo, Z0, Z2, Rf, Rg);
if ax.fam == 1
    R = Rf + Rg; V1 = E*(Z0 + Z2 + 3*R)/(Z1 + Z2 + Z0 + 3*R);
else
    Z2p = Z2 + Rf; Z0p = Z0 + Rf + 3*Rg; Zp = Z2p*Z0p/(Z2p + Z0p);
    V1 = E*(Rf + Zp)/(Z1 + Rf + Zp);
end
t = ((0:round(dur*fs)-1)' + toff)/fs;
V = r2_sector_wave(t, V1, u2*V1, f, u0*V1);
end
