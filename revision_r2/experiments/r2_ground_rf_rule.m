function [lab, info] = r2_ground_rf_rule(u0, u2)
% R2_GROUND_RF_RULE  Regla de fase de la rama de tierra basada en la relacion
% entre las componentes homopolar e inversa (T23b). Entradas: u0 = V0/V1 y
% u2 = V2/V1 (complejos, referidos a la fase a, secuencia abc).
%
% Fundamento (r2_ground_rf_seqlaw). Fase singular k: la fase en falta en L-G,
% la fase sana en L-L-G; k = 0, 1, 2 para a, b, c.
%  - beta = arg(V2/V1) vale 120 k + 180 + psi en L-G y 120 k + psi en L-L-G.
%    El giro psi va de 0 (falta franca) a arg Z2 (resistencia alta): las seis
%    clases tienen bases cada 60 grados (B-C-G 0, C-G 60, C-A-G 120, A-G 180,
%    A-B-G 240, B-G 300) y giran en el mismo sentido.
%  - alfa = arg(V0/V2) vale 120 k + arg(Z0/Z2) en L-G, sin depender de Rf. En
%    L-L-G vale 120 k + arg(Z0/Z2) + arg(Z2+Rf) - arg Z2 + arg Z0 - arg(Z0+Rf+3Rg).
% Regla:
%  A. Si |u2| >= 1/sqrt(2) (redondez del lugar en el plano de Kirchhoff
%     <= 0.172) el giro de beta es menor de 45 grados con cualquier X/R: clase
%     por el sextante de beta, [base - 15, base + 45). No usa u0. Cubre la
%     L-L-G con las dos fases unidas y resistencia comun a tierra (V2 = V1).
%  B. Si no: eje k por el sector de alfa (centros 0, 120, 240; +-60) y
%     familia por el semiplano de beta_k = beta - 120 k: (-45, 135) -> L-L-G
%     con fase sana k; el resto -> L-G en la fase k.
% Los cortes (1/sqrt(2), 45, 60, 135) salen de las cotas del giro, no de un
% ajuste a faltas francas.
LG  = {'A-G','B-G','C-G'};            % fase en falta k = 0, 1, 2
LLG = {'B-C-G','C-A-G','A-B-G'};      % fase sana  k = 0, 1, 2
SEXT = {'B-C-G','C-G','C-A-G','A-G','A-B-G','B-G'};   % bases 0, 60, ..., 300
wr = @(x) mod(x + 180, 360) - 180;
alfa = rad2deg(angle(u0/u2)); beta = rad2deg(angle(u2));
info = struct('alfa_deg', alfa, 'beta_deg', beta, 'k', NaN, 'alfa_k', NaN, 'beta_k', NaN, 'paso', '');
lab = '';
if ~isfinite(beta), return, end
if abs(u2) >= 1/sqrt(2)
    q = mod(floor((beta + 15)/60), 6);
    lab = SEXT{q+1}; info.paso = 'A'; info.beta_k = wr(beta - 60*q); return
end
if ~isfinite(alfa), return, end
k = mod(round(alfa/120), 3); bk = wr(beta - 120*k);
info.k = k; info.alfa_k = wr(alfa - 120*k); info.beta_k = bk; info.paso = 'B';
if bk > -45 && bk < 135, lab = LLG{k+1}; else, lab = LG{k+1}; end
end
