function r = r2_lab_signals(matfile)
% R2_LAB_SIGNALS  Carga un registro del laboratorio (formato OPAL de 250 filas) y
% construye las tres senales de tension que se clasifican en la ronda 2.
%
% Las expresiones son las de la parte 2 de experiments\r2_lab_dy_check.m (T12), que
% es un script y no se puede llamar; r2_lab_autokontaktor comprueba que con ellas se
% reproduce resultados\T12\lab_dy_classification.csv en los 12 registros L-L.
%
%   r.sig.conv1_ref_red  conv1 referida al lado de red a traves del grupo Dy11
%                        (tensiones compuestas / sqrt(3): -T_Dy1 * v), en pu. PRINCIPAL
%   r.sig.box2           meas_box2 (lado de red, medida directa junto a la falta), pu
%   r.sig.conv1_fase     conv1 fase-estrella en pu (filas 59:61), lo publicado
% Las tres con la permutacion de canales [3 2 1] (orden a-b-c del clasificador).
% Tambien devuelve los canales en bruto: Vc, Vb (V), Ic, Ib (A), Psi, Vdc, ICH (Egston).
d = load(matfile); R = d.results; L = d.legendText;
has = @(i, s) contains(string(L{i}), s);
assert(has(62, 'From_conv1') && has(62, 'V_abc_meas_SI') && has(183, 'From_meas_box2') && has(183, 'V_abc_meas_SI') && ...
       has(54, 'From_conv1') && has(54, 'I_abc_meas_SI') && has(59, 'V_abc_meas_pu') && has(50, 'P_si') && has(58, 'Vdc_meas_SI') && ...
       has(179, 'From_meas_box2') && has(179, 'I_abc_meas_SI') && has(44, 'CDA_current.ICH2') && has(45, 'CDA_current.ICH3'), ...
       'canales inesperados en %s', matfile);
r.fs = 1e4; r.f = 50; r.PERM = [3 2 1];
r.VB = 400*sqrt(2)/sqrt(3);                     % base pu de tension del registro (pico fase-neutro)
C = [0 1 0; 0 0 1; 1 0 0]; T1 = (eye(3) - C')/sqrt(3);
r.t = R(1,:)'; r.t = r.t - r.t(1);
r.Vc = R(62:64,:)'; r.Vcpu = R(59:61,:)'; r.Vb = R(183:185,:)'; r.Ic = R(54:56,:)'; r.Ib = R(179:181,:)';
r.Psi = R(50,:)'; r.Vdc = R(58,:)'; r.ICH = R(43:48,:)'; r.Irated = R(103,1); r.Urated = R(102,1);
r.sig.conv1_fase = r.Vcpu(:, r.PERM);
x = ((-T1)*r.Vc')'/r.VB; r.sig.conv1_ref_red = x(:, r.PERM);
x = r.Vb/r.VB; r.sig.box2 = x(:, r.PERM);
% matriz de tensiones conv1 = Mv * box2 en el orden registrado (T12): se usa para
% referir las corrientes de conv1 al lado de red (i_red = Mv' * i_conv, sin perdidas)
r.Mv = [-1 1 0; 0 -1 1; 1 0 -1]/sqrt(3);
end
