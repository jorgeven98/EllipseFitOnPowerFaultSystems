%% R2_SECTOR_ANALYTIC  (T21; peticiones R1.2 y X4 de la ronda 2 de TPWRD)
% Barrido analitico de los centros de sector de las faltas L-L con los ocho
% factores del revisor 1 (X/R de la fuente, resistencia de falta, fortaleza de
% la fuente, angulo de potencia prefalta, frecuencia, control del convertidor,
% aportacion remota y errores de fase desiguales entre canales).
%
% Todas las comprobaciones "con el clasificador" llaman a
% src/fault_classifiers/ga_fault_classifier_pure.m sin modificarlo, ventana a
% ventana (W = 20 muestras a 4 kHz, paso 5). La decision por evento usa la
% regla de persistencia M = 10 sobre esas etiquetas (r2_sector_event.m).
%
% Requiere haber ejecutado antes r2_sector_paper_network.m (genera
% paper_network_z2.csv en la carpeta de trabajo).
%
% Ejecutar desde ronda2\resultados\T21 con R2024b:
%   matlab -batch "addpath('<RAIZ>\experiments'); r2_sector_analytic"
% Salidas: CSV en la carpeta actual y ronda2\textos\T21_tabla.tex.
%
% Convenio: v_k(t) = Re{V_k e^{jwt}}, secuencia abc. Entonces
%   v_alfa + j v_beta = V1 e^{jwt} + conj(V2) e^{-jwt}
% y la orientacion del eje mayor es (arg V1 - arg V2)/2 (mod 180). En el marco
% del clasificador:  theta = 165 - (1/2) arg(V2/V1)  (mod 180).
% Rf es SIEMPRE la resistencia total entre las dos fases en falta.

RAIZ = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(RAIZ,'src','fault_classifiers'));
addpath(fullfile(RAIZ,'src','utils'));
warning('off','all'); rng(21);
G.fs = 4000; G.f0 = 50; G.W = 20; G.STRIDE = 5; G.M = 10;
TEXDIR = fullfile(RAIZ,'_revision_TPWRD','ronda2','textos');
QMIN_FIX = (1-0.99)/(1+0.99);          % puerta fija de redondez 0.99  -> |V2/V1| minimo
QMIN_AMB = (1-0.97)/(1+0.97);          % ejemplo adaptativo: ambiente con V2 = 1 % (r_amb = 0.98, Delta = 0.01)
SUM = {};                              % filas de la tabla resumen de los 8 factores
tic;

%% ===================== A. CONVENCION ANGULAR ===============================
fprintf('=== A. Convencion angular ===\n');
% A1. rejilla de |q| y arg q, con un origen de tiempos arbitrario (arg V1 = 37 grados)
absq = [0.01 0.02 0.05 0.1 0.2 0.4 0.6 0.8 0.95 1];
argq = -175:10:175;      % evita caer exactamente en una frontera de sector
rows = {};
for m = absq
    for ag = argq
        V1 = exp(1j*deg2rad(37)); V2 = m*exp(1j*deg2rad(ag))*V1;
        s = statwin(V1, V2, G.f0, G);
        tp = thpred(V2/V1);
        tsum = mod(165 + 0.5*rad2deg(angle(V1)+angle(V2)), 180);   % formula con la SUMA (linea 573), constante ajustada a la B-C franca con arg V1 = 0
        rows(end+1,:) = {m, ag, tp, s.theta, wrap90(s.theta-tp), s.thmin, s.thmax, wrap90(tsum-s.theta), ...
            (1-m)/(1+m), s.round, s.label, sector(tp)}; %#ok<SAGROW>
    end
end
TA1 = cell2table(rows,'VariableNames',{'absq','argq_deg','theta_pred_deg','theta_cls_deg','err_deg','theta_cls_min_deg', ...
    'theta_cls_max_deg','err_sum_formula_deg','roundness_pred','roundness_cls','label_cls','label_pred'});
writetable(TA1,'conv_grid.csv');
ok = TA1.absq < 1;      % en |q| = 1 la trayectoria es una recta (rama de linea), se da aparte
fprintf('A1  %d casos: max |theta_cls - theta_pred| = %.3g grados (|q|<1), %.3g (|q|=1); etiquetas distintas: %d\n', ...
    height(TA1), max(abs(TA1.err_deg(ok))), max(abs(TA1.err_deg(~ok))), sum(~strcmp(TA1.label_cls,TA1.label_pred)));
fprintf('    formula con la suma: error max %.1f grados\n', max(abs(TA1.err_sum_formula_deg)));
A1max = max(abs(TA1.err_deg));

% A2. desplazamiento comun de fase de V1 y V2 (cambio del origen de tiempos)
rows = {};
for al = [0 30 70 110 150 200 290]
    V1 = exp(1j*deg2rad(al)); V2 = 0.5*exp(1j*deg2rad(-120))*V1;      % A-B, |q| = 0.5
    s = statwin(V1, V2, G.f0, G);
    rows(end+1,:) = {al, s.theta, thpred(V2/V1), mod(165 + 0.5*rad2deg(angle(V1)+angle(V2)),180), s.label}; %#ok<SAGROW>
end
TA2 = cell2table(rows,'VariableNames',{'common_shift_deg','theta_cls_deg','theta_difference_formula_deg','theta_sum_formula_deg','label_cls'});
writetable(TA2,'conv_shift.csv'); disp(TA2);

% A3. los tres pares: saltos de 120 grados en arg(V2/V1) dan 60 grados en theta
rows = {};
PAIRS = {'B-C',0; 'A-B',-120; 'C-A',120};
for p = 1:3
    for m = [1 0.5 0.1]
        V1 = exp(1j*0.3); V2 = m*exp(1j*deg2rad(PAIRS{p,2}))*V1;
        s = statwin(V1, V2, G.f0, G);
        rows(end+1,:) = {PAIRS{p,1}, PAIRS{p,2}, m, thpred(V2/V1), s.theta, s.label}; %#ok<SAGROW>
    end
end
TA3 = cell2table(rows,'VariableNames',{'pair','argq_bolted_deg','absq','theta_pred_deg','theta_cls_deg','label_cls'});
writetable(TA3,'conv_pairs.csv'); disp(TA3);

% A4. mismo caso en pu y en voltios (el clasificador tiene un umbral absoluto b < 0.01)
rows = {};
for sc = [1 231 20e3*sqrt(2/3)]
    for m = [0.02 0.42 0.98 1]
        V1 = sc*exp(1j*0.3); V2 = m*V1*exp(1j*deg2rad(40));
        s = statwin(V1, V2, G.f0, G);
        rows(end+1,:) = {sc, m, thpred(V2/V1), s.theta, wrap90(s.theta-thpred(V2/V1)), s.round, s.label, s.pline}; %#ok<SAGROW>
    end
end
TA4 = cell2table(rows,'VariableNames',{'scale_V_peak','absq','theta_pred_deg','theta_cls_deg','err_deg','roundness_cls','label_cls','pct_line_branch'});
writetable(TA4,'conv_units.csv'); disp(TA4);

% A5. tercer armonico de secuencia directa (lo generan los convertidores con
%     control en dq ante desequilibrio): la orientacion del lugar geometrico la
%     fija  q_ef = q + (V3/V1) e^{-j 2 arg V1}  a primer orden
rows = {};
tt = (0:199)'/G.fs;
for c = 1:6
    switch c
        case 1, V1 = 0.989*exp(1j*deg2rad(-71.8)); q = 0.0139*exp(1j*deg2rad(-67.1)); h3 = 0.0213*exp(1j*deg2rad(-3.9))/V1; % DG1, RL_AB medido
        case 2, V1 = 1; q = 0.02; h3 = 0.02*exp(1j*deg2rad(90));
        case 3, V1 = exp(1j*1); q = 0.05*exp(1j*deg2rad(50)); h3 = 0.02*exp(1j*deg2rad(-140));
        case 4, V1 = 1; q = 0.2*exp(1j*deg2rad(-100)); h3 = 0.02*exp(1j*deg2rad(30));
        case 5, V1 = 1; q = 0; h3 = 0.03*exp(1j*deg2rad(200));
        case 6, V1 = exp(-1j*0.7); q = 0.4*exp(1j*deg2rad(60)); h3 = 0.05*exp(1j*deg2rad(10));
    end
    V = r2_sector_wave(tt, V1, q*V1, G.f0) + r2_sector_wave(tt, h3*V1, 0, 3*G.f0);
    S = r2_sector_scan(V, tt, G.W, G.STRIDE, G.fs, G.f0);
    qef = q + h3*exp(-1j*2*angle(V1));
    rows(end+1,:) = {abs(q), rad2deg(angle(q)), abs(h3), abs(qef), rad2deg(angle(qef)), thpred(q), thpred(qef), cmed(S.theta), ...
        min(wrap90(S.theta-thpred(qef))), max(wrap90(S.theta-thpred(qef)))}; %#ok<SAGROW>
end
TA5 = cell2table(rows,'VariableNames',{'absq','argq_deg','absV3_over_V1','absq_eff','argq_eff_deg','theta_from_q_deg','theta_from_q_eff_deg', ...
    'theta_cls_median_deg','win_dev_min_deg','win_dev_max_deg'});
writetable(TA5,'conv_h3.csv'); disp(TA5);

%% ===================== B. PUNTO DE FALTA ===================================
fprintf('=== B. Punto de falta ===\n');
% B1. contraejemplo de Codex: Z1 = Z2, sin carga, medida en el punto de falta
PH = [60 70 75 85]; RHO = [0 1 2 2.1 3 3.5 5 6 10];
CODEX85 = [1 0 165 165; 0.67817 42.5 143.75 143.75; 0.43239 59.485 135.26 135.26; 0.41609 60.512 134.74 134.74; ...
           0.30827 67.116 131.44 131.44; 0.26861 69.479 130.26 130.26; 0.19291 73.92 128.04 128.04; 0.16212 75.706 127.15 127.15; 0.098656 79.36 125.32 125.32];
CODEXLAB = {'B-C','B-C','B-C','C-A','C-A','C-A','C-A','C-A','C-A'};
rows = {};
for ph = PH
    for k = 1:numel(RHO)
        Z2 = exp(1j*deg2rad(ph)); Rf = RHO(k);
        V1 = (Z2+Rf)/(2*Z2+Rf); V2 = Z2/(2*Z2+Rf);
        s = statwin(V1, V2, G.f0, G);
        q = V2/V1; cx = [NaN NaN NaN NaN]; cl = '';
        if ph == 85, cx = CODEX85(k,:); cl = CODEXLAB{k}; end
        rows(end+1,:) = {ph, Rf, abs(q), rad2deg(angle(q)), thpred(q), s.theta, wrap90(s.theta-thpred(q)), s.pct('B-C'), s.label, ...
            cx(1), cx(2), cx(4), cl}; %#ok<SAGROW>
    end
end
TB1 = cell2table(rows,'VariableNames',{'argZ2_deg','Rf_over_absZ2','absq','argq_deg','theta_pred_deg','theta_cls_median_deg','err_deg', ...
    'BC_percent','modal_label','codex_absq','codex_argq_deg','codex_theta_median_deg','codex_label'});
writetable(TB1,'faultpoint_codex36.csv');
k85 = TB1.argZ2_deg == 85;
fprintf('B1  36 casos: max |theta_cls - theta_pred| = %.3g grados\n', max(abs(TB1.err_deg)));
fprintf('    frente al log de Codex (85 grados): max dif |q| %.2g, arg q %.2g grados, theta %.2g grados; etiquetas iguales: %d de 9\n', ...
    max(abs(TB1.absq(k85)-TB1.codex_absq(k85))), max(abs(TB1.argq_deg(k85)-TB1.codex_argq_deg(k85))), ...
    max(abs(TB1.theta_cls_median_deg(k85)-TB1.codex_theta_median_deg(k85))), sum(strcmp(TB1.modal_label(k85),TB1.codex_label(k85))));
disp(TB1(k85,1:9));

% B2. comprobacion independiente en el dominio de fases (abc): red con dos
%     fuentes desfasadas, carga y falta b-c con Rf. Se comprueba que
%     V2/V1 = Z2th/(Z2th+Rf) en F para cualquier angulo prefalta y carga, y que
%     el solver de secuencias reproduce V1 y V2 en el punto de medida M.
a = exp(1j*2*pi/3); Aabc = [1; a^2; a];
rows = {}; emax = 0; emaxM = 0;
for trial = 1:300
    Zs = (0.05+rand)*exp(1j*deg2rad(20+70*rand));  Zc = (0.05+rand)*exp(1j*deg2rad(20+70*rand));
    Zm = (0.02+0.5*rand)*exp(1j*deg2rad(10+80*rand)); Zl = (1+5*rand)*exp(1j*deg2rad(40*rand));
    Rf = 10^(-2+3*rand); del = deg2rad(-60+120*rand); Ec = (0.9+0.2*rand)*exp(1j*del);
    % nodos: M (3 fases), F (3 fases). Carga Zl en F. Fuente de red 1/_0 tras Zs en F; fuente Ec tras Zc en M.
    I3 = eye(3); Gf = (1/Rf)*[0 0 0; 0 1 -1; 0 -1 1];
    Y = [ (1/Zc+1/Zm)*I3, -(1/Zm)*I3; -(1/Zm)*I3, (1/Zm+1/Zs+1/Zl)*I3 + Gf ];
    J = [ (Ec/Zc)*Aabc; (1/Zs)*Aabc ];
    Vn = Y\J; VM = Vn(1:3); VF = Vn(4:6);
    sq = @(v) [(v(1)+a*v(2)+a^2*v(3))/3, (v(1)+a^2*v(2)+a*v(3))/3];
    sF = sq(VF); sM = sq(VM);
    Zfg = 1/(1/Zs+1/Zl);                                    % red + carga vistas desde F
    Z2th = 1/(1/Zfg + 1/(Zm+Zc));
    e1 = abs(sF(2)/sF(1) - Z2th/(Z2th+Rf));
    P = struct('Eg', Zfg/Zs, 'Z_FG', Zfg, 'Z_MF', Zm, 'Rf', Rf, 'conv','gfm', 'Z2c', Zc, 'Ec', Ec, 'Zc', Zc);
    R = r2_sector_seqsolve(P);
    e2 = max(abs([R.V1M-sM(1), R.V2M-sM(2), R.V1F-sF(1), R.V2F-sF(2)]));
    emax = max(emax,e1); emaxM = max(emaxM,e2);
    if trial <= 20
        rows(end+1,:) = {trial, abs(Z2th), rad2deg(angle(Z2th)), Rf, rad2deg(del), abs(sF(2)/sF(1)), rad2deg(angle(sF(2)/sF(1))), ...
            abs(Z2th/(Z2th+Rf)), rad2deg(angle(Z2th/(Z2th+Rf))), e1, e2}; %#ok<SAGROW>
    end
end
TB2 = cell2table(rows,'VariableNames',{'trial','absZ2th','argZ2th_deg','Rf','prefault_angle_deg','absq_abc','argq_abc_deg','absq_formula','argq_formula_deg','err_q','err_solver'});
writetable(TB2,'faultpoint_abc_check.csv');
fprintf('B2  300 redes aleatorias en abc: max |q_abc - Z2th/(Z2th+Rf)| = %.2g; max error del solver de secuencias = %.2g\n', emax, emaxM);
B2err = [emax emaxM];

% B3. curvas de la region de cruce en (arg Z2, Rf/|Z2|)
phi = (30:0.25:100)';
rc = rho_iso(phi,60); r10 = rho_iso(phi,20); r20 = rho_iso(phi,40);       % Rf/|Z2| que da un giro de 30, 10 y 20 grados
qc = sind(phi-60)./sind(phi); qc(phi<=60) = NaN;
rdet = @(qm) -cosd(phi) + sqrt(cosd(phi).^2 - 1 + 1/qm^2);                  % |Z2/(Z2+Rf)| = qm
TB3 = table(phi, tand(min(phi,89.999)), rc, qc, r10, r20, rdet(QMIN_FIX), rdet(QMIN_AMB), phi/2, ...
    'VariableNames',{'argZ2_deg','X_over_R','rho_cross','absq_at_cross','rho_shift10','rho_shift20','rho_det_fixed','rho_det_ambient1pct','max_shift_deg'});
writetable(TB3,'crossing_curves.csv');
% poligono de la region de seleccion erronea detectable (para la figura)
kk = phi > 60.2 & phi <= 90;
rd = rdet(QMIN_FIX); poly = [phi(kk), min(rc(kk),rd(kk)); flipud(phi(kk)), flipud(rd(kk))];
writetable(array2table(poly,'VariableNames',{'argZ2_deg','rho'}),'fig_region_wrong.csv');
k9 = phi <= 90; polu = [phi(k9), rd(k9); flipud(phi(k9)), 1e4*ones(sum(k9),1)];
writetable(array2table(polu,'VariableNames',{'argZ2_deg','rho'}),'fig_region_undetected.csv');

% B4. el clasificador sobre la rejilla (arg Z2, Rf/|Z2|)
PHm = [30 45 55 60 62 65 70 75 80 83.13 85 88 90];
RHm = [0 0.05 0.1 0.2 0.5 1 1.5 1.8 2 2.1 2.2 2.5 3 3.5 4 5 7 10 20 50 100 150 180 190 195 205 250 500];
rows = {};
for ph = PHm
    for rh = RHm
        Z2 = exp(1j*deg2rad(ph)); V1 = (Z2+rh)/(2*Z2+rh); V2 = Z2/(2*Z2+rh); q = V2/V1;
        s = statwin(V1, V2, G.f0, G);
        lp = sector(thpred(q)); if (1-abs(q))/(1+abs(q)) > 0.99, lp = 'A-B-C'; end
        rows(end+1,:) = {ph, rh, abs(q), rad2deg(angle(q)), -0.5*rad2deg(angle(q)), thpred(q), s.theta, wrap90(s.theta-thpred(q)), ...
            s.round, lp, s.label, s.pct('B-C'), abs(V1)}; %#ok<SAGROW>
    end
end
TB4 = cell2table(rows,'VariableNames',{'argZ2_deg','Rf_over_absZ2','absq','argq_deg','shift_deg','theta_pred_deg','theta_cls_median_deg','err_deg', ...
    'roundness_cls','label_pred','label_cls','BC_percent','V1_pu'});
writetable(TB4,'crossing_map_classifier.csv');
det = ~strcmp(TB4.label_cls,'A-B-C');
fprintf('B4  %d puntos: etiqueta del clasificador distinta de la analitica en %d; max |err theta| (detectadas) = %.3g grados\n', ...
    height(TB4), sum(~strcmp(TB4.label_pred,TB4.label_cls)), max(abs(TB4.err_deg(det))));
disp(TB4(~strcmp(TB4.label_pred,TB4.label_cls), {'argZ2_deg','Rf_over_absZ2','absq','theta_pred_deg','theta_cls_median_deg','roundness_cls','label_pred','label_cls'}));
B4mis = sum(~strcmp(TB4.label_pred,TB4.label_cls)); B4err = max(abs(TB4.err_deg(det)));
km = ismember(TB4.argZ2_deg,[45 60 70 85]) & ismember(TB4.Rf_over_absZ2,[0.5 1 2 3 5 10 20 50 100]);
writetable(TB4(km & strcmp(TB4.label_cls,'B-C'), {'argZ2_deg','Rf_over_absZ2','shift_deg'}), 'fig_markers_cls_ok.csv');
writetable(TB4(km & strcmp(TB4.label_cls,'C-A'), {'argZ2_deg','Rf_over_absZ2','shift_deg'}), 'fig_markers_cls_wrong.csv');

%% ===================== C. MAGNITUDES DE INGENIERIA =========================
fprintf('=== C. Magnitudes de ingenieria ===\n');
% C1. red del paper: Z2th medida en el bus de falta (r2_sector_paper_network.m)
TZ = readtable('paper_network_z2.csv','Delimiter',',');
iz = find(strcmp(TZ.origen,'media RL (Pex)'));
Z2p = TZ.R_ohm(iz) + 1j*TZ.X_ohm(iz); php = rad2deg(angle(Z2p));
Zbase = 231/40;
Rf = logspace(-3, 3, 241)';
q = Z2p./(Z2p+Rf);
lab = arrayfun(@(x) string(sector(thpred(x))), q); lab((1-abs(q))./(1+abs(q)) > 0.99) = "A-B-C";
TC1 = table(Rf, Rf/abs(Z2p), Rf/Zbase, abs(q), rad2deg(angle(q)), -0.5*rad2deg(angle(q)), (1-abs(q))./(1+abs(q)), lab, ...
    'VariableNames',{'Rf_total_ohm','Rf_over_absZ2','Rf_pu_converter_base','absq','argq_deg','shift_deg','roundness','label_BC_fault'});
writetable(TC1,'eng_paper_network.csv');
rcp = sind(60)/sind(php-60); rdp = -cosd(php) + sqrt(cosd(php)^2 - 1 + 1/QMIN_FIX^2); rdp2 = -cosd(php) + sqrt(cosd(php)^2 - 1 + 1/QMIN_AMB^2);
fprintf('C1  red del paper: Z2th = %.3f ohm /_ %.2f grados (X/R = %.2f)\n', abs(Z2p), php, tand(php));
fprintf('    cruce: Rf > %.2f |Z2| = %.2f ohm (%.2f pu base convertidor), con |V2/V1| = %.3f\n', rcp, rcp*abs(Z2p), rcp*abs(Z2p)/Zbase, sind(php-60)/sind(php));
fprintf('    deja de detectarse: Rf > %.0f ohm (redondez 0.99) o > %.0f ohm (umbral 0.97)\n', rdp*abs(Z2p), rdp2*abs(Z2p));
PAPER = struct('absZ2',abs(Z2p),'phi',php,'rho_cross',rcp,'Rf_cross',rcp*abs(Z2p),'q_cross',sind(php-60)/sind(php),'Rf_det',rdp*abs(Z2p),'Rf_det2',rdp2*abs(Z2p));

% C2. redes realistas. Rf de cruce y comparacion con la resistencia de arco
%     (Warrington: R = 28707 L / I^1.4, L en m, I en A eficaces).
%     Parametros de red: CIGRE TB 575 (valores de la implementacion de
%     pandapower, cigre_networks.py) e IEC 60909-0 (R/X del equivalente de red).
%     Las separaciones entre fases y los niveles de cortocircuito de AT son
%     HIPOTESIS de trabajo, no datos normativos.
N = {};  % nombre, Un [kV], Z2 [ohm], separacion entre fases [m], fuente
zfeed = @(Un,Sk,rx) (Un^2/Sk)*(rx + 1j)/sqrt(1+rx^2);                % Un en kV, Sk en MVA -> ohm
ztr = @(Un,Sn,vk,vkr) (Un^2/Sn)*(vkr + 1j*sqrt(vk^2-vkr^2))/100;
Zlvbus = zfeed(0.4,100,1) + ztr(0.4,0.5,4.123106,1.0);
Zmvbus = zfeed(20,5000,0.1) + ztr(20,25,12.00107,0.16);
N(end+1,:) = {'BT 0.4 kV, bornes del trafo MT/BT (CIGRE BT)', 0.4, Zlvbus, 0.3, 'CIGRE TB 575 (pandapower)'};
N(end+1,:) = {'BT 0.4 kV, a 200 m de cable UG1', 0.4, Zlvbus + 0.2*(0.162+1j*0.0832), 0.3, 'CIGRE TB 575 (pandapower)'};
N(end+1,:) = {'Red del paper 0.283 kV, bus de falta (medido)', 0.2829, Z2p, 0.3, 'paper_network_z2.csv'};
N(end+1,:) = {'MT 20 kV, barras de subestacion (CIGRE MT)', 20, Zmvbus, 1.0, 'CIGRE TB 575 (pandapower)'};
N(end+1,:) = {'MT 20 kV, a 2 km de linea aerea', 20, Zmvbus + 2*(0.510+1j*0.366), 1.0, 'CIGRE TB 575 (pandapower)'};
N(end+1,:) = {'MT 20 kV, a 5 km de linea aerea', 20, Zmvbus + 5*(0.510+1j*0.366), 1.0, 'CIGRE TB 575 (pandapower)'};
N(end+1,:) = {'MT 20 kV, a 5 km de cable', 20, Zmvbus + 5*(0.501+1j*0.716), 1.0, 'CIGRE TB 575 (pandapower)'};
N(end+1,:) = {'AT 220 kV, Sk = 10 GVA (reactancia pura)', 220, zfeed(220,10000,0), 6.5, 'IEC 60909-0; Sk hipotesis'};
N(end+1,:) = {'AT 220 kV, Sk = 2 GVA + 50 km de linea', 220, zfeed(220,2000,0) + 50*(0.0653+1j*0.398), 6.5, 'IEC 60909-0; CIGRE TB 575; Sk hipotesis'};
N(end+1,:) = {'MAT 380 kV, Sk = 30 GVA (reactancia pura)', 380, zfeed(380,30000,0), 10, 'IEC 60909-0; Sk hipotesis'};
N(end+1,:) = {'MAT 380 kV, Sk = 5 GVA + 100 km de linea', 380, zfeed(380,5000,0) + 100*(0.0328+1j*0.312), 10, 'IEC 60909-0; CIGRE TB 575; Sk hipotesis'};
rows = {};
for i = 1:size(N,1)
    Un = 1e3*N{i,2}; Z2 = N{i,3}; L = N{i,4}; ph = rad2deg(angle(Z2));
    if ph > 60, rcx = sind(60)/sind(ph-60); qcx = sind(ph-60)/sind(ph); else, rcx = Inf; qcx = NaN; end
    rdx = -cosd(ph) + sqrt(cosd(ph)^2 - 1 + 1/QMIN_FIX^2);
    Ibolt = Un/abs(2*Z2);                                   % corriente de la falta L-L franca (Z1 = Z2)
    arc = zeros(1,2); iarc = zeros(1,2); sh = zeros(1,2); k = 0;
    for el = [1 3]                                          % arco de longitud L y alargado a 3 L
        k = k+1;
        % g(I) = I |2 Z2 + R(I)| - Un tiene dos raices o ninguna (caracteristica
        % descendente del arco); la fisica es la de mayor corriente. Sin raiz, el
        % arco de esa longitud no se sostiene a esa tension segun la formula.
        g = @(I) I.*abs(2*Z2 + 28707*el*L./I.^1.4) - Un;
        Im = fminbnd(g, 1e-5*Ibolt, Ibolt);
        if N{i,2} < 1 || g(Im) > 0
            arc(k) = NaN; iarc(k) = NaN; sh(k) = NaN;      % BT: fuera del rango de la formula de Warrington
        else
            Is = fzero(g, [Im, Ibolt]);
            arc(k) = 28707*el*L/Is^1.4; iarc(k) = Is;
            sh(k) = -0.5*rad2deg(angle(Z2/(Z2+arc(k))));
        end
    end
    rows(end+1,:) = {N{i,1}, N{i,2}, real(Z2), imag(Z2), abs(Z2), ph, imag(Z2)/real(Z2), 1e-3*Ibolt, rcx, rcx*abs(Z2), qcx, rdx*abs(Z2), ...
        L, arc(1), arc(1)/abs(Z2), sh(1), 1e-3*iarc(1), arc(2), arc(2)/abs(Z2), sh(2), N{i,5}}; %#ok<SAGROW>
end
TC2 = cell2table(rows,'VariableNames',{'network','Un_kV','R2_ohm','X2_ohm','absZ2_ohm','argZ2_deg','X_over_R','Ik_LL_bolted_kA','rho_cross','Rf_cross_ohm', ...
    'absq_at_cross','Rf_detect_limit_ohm','phase_spacing_m','Rarc_ohm','Rarc_over_absZ2','shift_arc_deg','Iarc_kA','Rarc3L_ohm','Rarc3L_over_absZ2','shift_arc3L_deg','source'});
writetable(TC2,'eng_networks.csv');
disp(TC2(:,{'network','absZ2_ohm','argZ2_deg','X_over_R','Rf_cross_ohm','absq_at_cross','Rf_detect_limit_ohm','Rarc_ohm','shift_arc_deg','Rarc3L_ohm','shift_arc3L_deg'}));
ka = ~isnan(TC2.Rarc_over_absZ2);
writetable(TC2(ka,{'argZ2_deg','Rarc_over_absZ2','Rarc3L_over_absZ2','Un_kV'}),'fig_markers_arc.csv');
writetable(table([php; php; php], [0.004; 0.006; 0.006]/abs(Z2p), 'VariableNames',{'argZ2_deg','Rf_over_absZ2'}),'fig_markers_dataset.csv');

%% ===================== D. LOS OCHO FACTORES ================================
fprintf('=== D. Ocho factores ===\n');
shiftF = @(ph,rh) -0.5*rad2deg(angle(exp(1j*deg2rad(ph))./(exp(1j*deg2rad(ph))+rh)));   % giro en el punto de falta (grados)
absqF  = @(ph,rh) abs(exp(1j*deg2rad(ph))./(exp(1j*deg2rad(ph))+rh));

% ---- F1. X/R de la fuente (punto de falta), para varios Rf/|Z2|
XR = [0.5 1 sqrt(3) 2 3 5 5.83 7 8.3 10 15 20 30 50]; RH = [0.1 0.5 1 2 5 10];
rows = {};
for xr = XR
    for rh = RH
        ph = atand(xr); [tc, lc] = clsF(ph, rh, G);
        rows(end+1,:) = {xr, ph, rh, absqF(ph,rh), shiftF(ph,rh), wrap90(tc-165), lc}; %#ok<SAGROW>
    end
end
T1 = cell2table(rows,'VariableNames',{'X_over_R','argZ2_deg','Rf_over_absZ2','absq','shift_pred_deg','shift_cls_deg','label_cls_BC_fault'});
writetable(T1,'f1_source_xr.csv');
e1 = errdet(T1.shift_pred_deg, T1.shift_cls_deg, T1.label_cls_BC_fault);
for rh = [0.1 1 10]
    k = T1.Rf_over_absZ2 == rh;
    SUM(end+1,:) = {1, 'Source X/R', sprintf('%.1f to %.0f', min(XR), max(XR)), sprintf('$R_f/|Z_2|=%g$', rh), max(abs(T1.shift_pred_deg(k))), ...
        any(abs(T1.shift_pred_deg(k)) > 30), 'f1_source_xr.csv', rh ~= 0.1}; %#ok<SAGROW>
end

% ---- F2. resistencia de falta (punto de falta), para varios X/R
RH2 = [0 logspace(-2, 2.7, 95)]; PH2 = [45 60 70 75 80 php 85 88 90];
rows = {};
for ph = PH2
    for rh = RH2
        rows(end+1,:) = {ph, tand(min(ph,89.999)), rh, absqF(ph,rh), shiftF(ph,rh), (1-absqF(ph,rh))/(1+absqF(ph,rh))}; %#ok<SAGROW>
    end
end
T2 = cell2table(rows,'VariableNames',{'argZ2_deg','X_over_R','Rf_over_absZ2','absq','shift_pred_deg','roundness'});
writetable(T2,'f2_fault_resistance.csv');
% verificacion con el clasificador en un subconjunto
rows = {};
for ph = [60 75 php 85 90]
    for rh = [0 0.1 0.3 1 2 3 5 10 30 100]
        [tc, lc] = clsF(ph, rh, G);
        rows(end+1,:) = {ph, rh, absqF(ph,rh), shiftF(ph,rh), wrap90(tc-165), lc}; %#ok<SAGROW>
    end
end
T2c = cell2table(rows,'VariableNames',{'argZ2_deg','Rf_over_absZ2','absq','shift_pred_deg','shift_cls_deg','label_cls_BC_fault'});
writetable(T2c,'f2_fault_resistance_classifier.csv');
e2 = errdet(T2c.shift_pred_deg, T2c.shift_cls_deg, T2c.label_cls_BC_fault);
k = T2.argZ2_deg == 45 & T2.Rf_over_absZ2 <= 100;
SUM(end+1,:) = {2, 'Fault resistance', '$R_f/|Z_2|$ from 0 to 100', '$X/R=1$', max(abs(T2.shift_pred_deg(k))), false, 'f2_fault_resistance.csv', true};
k = abs(T2.argZ2_deg - php) < 1e-9 & T2.Rf_over_absZ2 <= 100;
SUM(end+1,:) = {2, 'Fault resistance', '$R_f/|Z_2|$ from 0 to 100', sprintf('$X/R=%.1f$', tand(php)), max(abs(T2.shift_pred_deg(k))), true, 'f2_fault_resistance.csv', true};

% ---- F3. fortaleza de la fuente: |Z2| = 1/SCR (pu), con Rf fijo en pu. Interaccion Rf x SCR
SCR = [1.5 2 3 5 7 10 15 20 30 50]; RFPU = [0.01 0.05 0.1 0.2 0.5 1]; ph = atand(10);
rows = {};
for sc = SCR
    for rf = RFPU
        rh = rf*sc; [tc, lc] = clsF(ph, rh, G);
        rows(end+1,:) = {sc, rf, rh, absqF(ph,rh), shiftF(ph,rh), wrap90(tc-165), lc, sind(60)/sind(ph-60)/sc}; %#ok<SAGROW>
    end
end
T3 = cell2table(rows,'VariableNames',{'SCR','Rf_pu','Rf_over_absZ2','absq','shift_pred_deg','shift_cls_deg','label_cls_BC_fault','Rf_cross_pu'});
writetable(T3,'f3_source_strength.csv');   % es tambien la interaccion Rf x SCR
e3 = errdet(T3.shift_pred_deg, T3.shift_cls_deg, T3.label_cls_BC_fault);
for rf = [0.01 0.05 0.2]
    k = T3.Rf_pu == rf;
    SUM(end+1,:) = {3, 'Source strength', sprintf('SCR from %.1f to %.0f, $X/R=10$', min(SCR), max(SCR)), sprintf('$R_f=%g$~pu', rf), ...
        max(abs(T3.shift_pred_deg(k))), any(abs(T3.shift_pred_deg(k)) > 30), 'f3_source_strength.csv', rf ~= 0.05}; %#ok<SAGROW>
end

% ---- Caso base del modelo de bornes (pu sobre la potencia del convertidor)
%      red: SCR = 5, X/R = 10; linea M-G de 0.10 pu con X/R = 5; falta a mitad de linea
B.SCR = 5; B.xr_g = 10; B.zl = 0.10; B.xr_l = 5; B.m = 0.5; B.Imax = 1.1; B.gamma = -pi/2;
B.Zc = 0.15*exp(1j*atan(10));         % GFM pasivo: fuente tras 0.15 pu, X/R = 10
mk = @(SCR,xrg,zl,xrl,m) struct('Zg',(1/SCR)*exp(1j*atan(xrg)), 'Zl', zl*exp(1j*atan(xrl)), 'm', m);

% ---- F4. angulo de potencia prefalta
%      (a) en el punto de falta: no interviene (B2). (b) en bornes de un GFM
%      pasivo, con la f.e.m. interna adelantada delta. (c) en bornes de un GFL,
%      con el angulo gamma de la corriente inyectada.
rows = {};
net = mk(B.SCR,B.xr_g,B.zl,B.xr_l,B.m);
for rf = [0 0.05 0.2 0.5]
    for del = -40:5:40
        P = struct('Eg',1,'Z_FG',(1-net.m)*net.Zl+net.Zg,'Z_MF',net.m*net.Zl,'Rf',rf,'conv','gfm','Z2c',B.Zc,'Ec',exp(1j*deg2rad(del)),'Zc',B.Zc);
        R = r2_sector_seqsolve(P);
        [tc, lc] = clsQ(R.V1M, R.V2M, G);
        rows(end+1,:) = {'gfm', rf, del, -0.5*rad2deg(angle(R.qF)), -0.5*rad2deg(angle(R.qM)), wrap90(tc-165), abs(R.qM), lc}; %#ok<SAGROW>
    end
    for gam = -180:15:165
        P = struct('Eg',1,'Z_FG',(1-net.m)*net.Zl+net.Zg,'Z_MF',net.m*net.Zl,'Rf',rf,'conv','gfl','Z2c',Inf,'Imax',B.Imax,'gamma',deg2rad(gam));
        R = r2_sector_seqsolve(P);
        [tc, lc] = clsQ(R.V1M, R.V2M, G);
        rows(end+1,:) = {'gfl', rf, gam, -0.5*rad2deg(angle(R.qF)), -0.5*rad2deg(angle(R.qM)), wrap90(tc-165), abs(R.qM), lc}; %#ok<SAGROW>
    end
end
T4 = cell2table(rows,'VariableNames',{'converter','Rf_pu','angle_deg','shift_F_deg','shift_M_pred_deg','shift_M_cls_deg','absq_M','label_cls_BC_fault'});
writetable(T4,'f4_prefault_angle.csv');
e4 = errdet(T4.shift_M_pred_deg, T4.shift_M_cls_deg, T4.label_cls_BC_fault);
k = strcmp(T4.converter,'gfm');
dF = 0; for rf = [0 0.05 0.2 0.5], kk2 = k & T4.Rf_pu == rf; dF = max(dF, max(T4.shift_F_deg(kk2)) - min(T4.shift_F_deg(kk2))); end
SUM(end+1,:) = {4, 'Pre-fault power angle', '$\delta$ from $-40^\circ$ to $40^\circ$', 'fault point', dF, false, 'f4_prefault_angle.csv', true};
k4 = k & T4.Rf_pu == 0.05; k40 = k4 & T4.angle_deg == 0;
SUM(end+1,:) = {4, 'Pre-fault power angle', '$\delta$ from $-40^\circ$ to $40^\circ$', 'GFM terminals', ...
    max(abs(T4.shift_M_pred_deg(k4))), any(abs(T4.shift_M_pred_deg(k4)) > 30), 'f4_prefault_angle.csv', true};
k4 = strcmp(T4.converter,'gfl') & T4.Rf_pu == 0.05;
SUM(end+1,:) = {4, 'Pre-fault power angle', 'injected-current angle over $360^\circ$', 'GFL terminals', ...
    max(abs(T4.shift_M_pred_deg(k4))), any(abs(T4.shift_M_pred_deg(k4)) > 30), 'f4_prefault_angle.csv', false};

% ---- F5. frecuencia 47.5 a 51.5 Hz con el clasificador fijado a 50 Hz
FQ = 47.5:0.5:51.5;
CASE5 = {'Pex RL medido', 0.754*exp(-1j*deg2rad(1.43)); 'phi 85, rho 1.9', exp(1j*deg2rad(85))/(exp(1j*deg2rad(85))+1.9); ...
         'firma debil 5 %', 0.05*exp(1j*deg2rad(20)); 'firma debil 2 %', 0.02*exp(1j*deg2rad(-30))};
rows = {};
for c = 1:size(CASE5,1)
    q = CASE5{c,2}; s50 = statwin(1, q, 50, G);
    for fq = FQ
        s = statwin(1, q, fq, G);
        rows(end+1,:) = {CASE5{c,1}, abs(q), fq, thpred(q), s.theta, wrap90(s.theta - s50.theta), wrap90(s.thmin - s50.theta), wrap90(s.thmax - s50.theta), s.round, s.label}; %#ok<SAGROW>
    end
end
T5 = cell2table(rows,'VariableNames',{'case','absq','f_Hz','theta_pred_deg','theta_cls_median_deg','dtheta_vs_50Hz_deg','dtheta_win_min_deg','dtheta_win_max_deg','roundness_cls','label_cls'});
writetable(T5,'f5_frequency.csv');
% efecto de la frecuencia a traves de X = wL (arg Z2 y |Z2| cambian con f)
rows = {};
for ph0 = [70 php 85 88]
    for rh0 = [0.5 1 1.9 3]
        for fq = FQ
            Z = cosd(ph0) + 1j*sind(ph0)*fq/50;
            rows(end+1,:) = {ph0, rh0, fq, -0.5*rad2deg(angle(Z/(Z+rh0))) - shiftF(ph0,rh0)}; %#ok<SAGROW>
        end
    end
end
T5b = cell2table(rows,'VariableNames',{'argZ2_at_50Hz_deg','Rf_over_absZ2_at_50Hz','f_Hz','dshift_deg'});
writetable(T5b,'f5_frequency_reactance.csv');
k5 = T5.absq >= 0.05;
SUM(end+1,:) = {5, 'System frequency', '47.5 to 51.5 Hz, classifier fixed at 50 Hz', '$|V_2/V_1|\geq 0.05$', max(abs([T5.dtheta_win_min_deg(k5); T5.dtheta_win_max_deg(k5)])), false, 'f5_frequency.csv', true};
SUM(end+1,:) = {5, 'System frequency', '47.5 to 51.5 Hz', 'through $X=\omega L$', max(abs(T5b.dshift_deg)), false, 'f5_frequency_reactance.csv', false};

% ---- F6. control del convertidor (medida en bornes M, falta a mitad de linea)
%      Z2c de la inyeccion VDE-AR-N 4120 segun Haddadi et al. 2020 (tabla 1):
%      k = 2: |V2| 0.29, |I2| 0.53, I2 adelanta 85.0 a V2  -> Z2c = 0.547 /_ 95.0
%      k = 6: |V2| 0.24, |I2| 0.67, I2 adelanta 100.6 a V2 -> Z2c = 0.358 /_ 79.4
%      GFM con control de tension (dataset del paper, DG1/DG2): Z2c = 0.0147 pu /_ 144
CTRL = { ...
 'GFL, I2 = 0, reactive priority',   'gfl', Inf, -90; ...
 'GFL, I2 = 0, active current',      'gfl', Inf, 0; ...
 'GFL, I2 = 0, worst current angle', 'gfl', Inf, NaN; ...
 'GFL, k = 2 (VDE-AR-N 4120)',       'gfl', 0.29/0.53*exp(1j*deg2rad(180-85.0)), -90; ...
 'GFL, k = 6 (VDE-AR-N 4120)',       'gfl', 0.24/0.67*exp(1j*deg2rad(180-100.6)), -90; ...
 'GFM, source behind 0.15 pu',       'gfm', B.Zc, 0; ...
 'GFM, voltage-controlled (measured Z2)', 'gfm', 0.085/Zbase*exp(1j*deg2rad(144)), 0};
rows = {};
for xl = [5 0.5]                            % linea inductiva (AT) y resistiva (BT)
    net = mk(B.SCR,B.xr_g,B.zl,xl,B.m);
    for rf = [0 0.05 0.2 0.5]
        for c = 1:size(CTRL,1)
            P = struct('Eg',1,'Z_FG',(1-net.m)*net.Zl+net.Zg,'Z_MF',net.m*net.Zl,'Rf',rf,'conv',CTRL{c,2},'Z2c',CTRL{c,3}, ...
                'Imax',B.Imax,'gamma',deg2rad(CTRL{c,4}),'Ec',1,'Zc',B.Zc);
            if isnan(CTRL{c,4})                  % peor angulo de corriente
                best = 0;
                for gam = -180:5:175
                    P.gamma = deg2rad(gam); Rt = r2_sector_seqsolve(P);
                    if abs(angle(Rt.qM)) >= best, best = abs(angle(Rt.qM)); R = Rt; end
                end
            else
                R = r2_sector_seqsolve(P);
            end
            [tc, lc] = clsQ(R.V1M, R.V2M, G);
            rows(end+1,:) = {CTRL{c,1}, xl, rf, abs(R.Z2th), rad2deg(angle(R.Z2th)), -0.5*rad2deg(angle(R.qF)), ...
                -0.5*rad2deg(angle(R.k2)), 0.5*rad2deg(angle(R.k1)), -0.5*rad2deg(angle(R.qM)), wrap90(tc-165), abs(R.qM), lc}; %#ok<SAGROW>
        end
    end
end
T6 = cell2table(rows,'VariableNames',{'control','line_X_over_R','Rf_pu','absZ2th_pu','argZ2th_deg','shift_F_deg','shift_k2_deg','shift_k1_deg', ...
    'shift_M_pred_deg','shift_M_cls_deg','absq_M','label_cls_BC_fault'});
writetable(T6,'f6_converter_control.csv');
e6 = errdet(T6.shift_M_pred_deg, T6.shift_M_cls_deg, T6.label_cls_BC_fault);
disp(T6(T6.Rf_pu == 0.05, :));
kv = contains(T6.control,'voltage-controlled'); kb = T6.Rf_pu == 0.05;
k6 = kb & contains(T6.control,'I2 = 0');
SUM(end+1,:) = {6, 'Converter control', 'converter terminals', 'GFL with $I_2=0$', max(abs(T6.shift_M_pred_deg(k6))), any(abs(T6.shift_M_pred_deg(k6)) > 30), 'f6_converter_control.csv', true};
k6 = kb & contains(T6.control,'VDE');
SUM(end+1,:) = {6, 'Converter control', 'converter terminals', 'GFL with $I_2$ injection, $k=2$ and 6', max(abs(T6.shift_M_pred_deg(k6))), any(abs(T6.shift_M_pred_deg(k6)) > 30), 'f6_converter_control.csv', true};
k6 = kb & contains(T6.control,'behind');
SUM(end+1,:) = {6, 'Converter control', 'converter terminals', 'GFM behind 0.15~pu', max(abs(T6.shift_M_pred_deg(k6))), any(abs(T6.shift_M_pred_deg(k6)) > 30), 'f6_converter_control.csv', true};
k6 = kb & kv;
SUM(end+1,:) = {6, 'Converter control', 'converter terminals', 'GFM with voltage loop, measured $Z_2$', max(abs(T6.shift_M_pred_deg(k6))), any(abs(T6.shift_M_pred_deg(k6)) > 30), 'f6_converter_control.csv', true};

% ---- F7. aportacion remota: fortaleza del extremo remoto y posicion de la falta
rows = {};
for cv = {'gfl','gfm'}
    for sc = [1.5 2 3 5 10 20 50]
        for mm = [0.1 0.5 0.9]
            for rf = [0.05 0.2]
                net = mk(sc,B.xr_g,B.zl,B.xr_l,mm);
                P = struct('Eg',1,'Z_FG',(1-net.m)*net.Zl+net.Zg,'Z_MF',net.m*net.Zl,'Rf',rf,'conv',cv{1},'Z2c',Inf, ...
                    'Imax',B.Imax,'gamma',B.gamma,'Ec',1,'Zc',B.Zc);
                if strcmp(cv{1},'gfm'), P.Z2c = B.Zc; end
                R = r2_sector_seqsolve(P);
                [tc, lc] = clsQ(R.V1M, R.V2M, G);
                rows(end+1,:) = {cv{1}, sc, mm, rf, abs(R.Z2th), rad2deg(angle(R.Z2th)), rf/abs(R.Z2th), -0.5*rad2deg(angle(R.qF)), ...
                    -0.5*rad2deg(angle(R.qM)), wrap90(tc-165), abs(R.qM), lc}; %#ok<SAGROW>
            end
        end
    end
end
T7 = cell2table(rows,'VariableNames',{'converter','remote_SCR','fault_position_m','Rf_pu','absZ2th_pu','argZ2th_deg','Rf_over_absZ2th', ...
    'shift_F_deg','shift_M_pred_deg','shift_M_cls_deg','absq_M','label_cls_BC_fault'});
writetable(T7,'f7_remote_infeed.csv');
e7 = errdet(T7.shift_M_pred_deg, T7.shift_M_cls_deg, T7.label_cls_BC_fault);
for rf = [0.05 0.2]
    k = T7.Rf_pu == rf;
    SUM(end+1,:) = {7, 'Remote infeed', 'remote SCR from 1.5 to 50, fault at 10 to 90\,\% of the line', sprintf('$R_f=%g$~pu', rf), ...
        max(abs(T7.shift_M_pred_deg(k))), any(abs(T7.shift_M_pred_deg(k)) > 30), 'f7_remote_infeed.csv', rf == 0.05}; %#ok<SAGROW>
end

% ---- F8. errores de fase desiguales entre canales (0 a 2 grados, uno y dos canales)
%      Analitico: cada fasor de fase se retrasa delta_k; el angulo se evalua en
%      el plano de Kirchhoff (se descarta la V0 aparente). Clasificador: misma
%      senal con la componente comun retirada, para aislar el efecto sobre el
%      angulo; aparte, fraccion de ventanas que NO llega a la rama angular si
%      no se retira (puerta de secuencia cero, tarea T51/T52).
DL = [0.25 0.5 1 1.5 2]; QA = [0.754 0.42 0.2 0.1 0.05 0.02];
PAT = {'A',[1 0 0]; 'B',[0 1 0]; 'C',[0 0 1]; 'A+ B-',[1 -1 0]; 'B+ C-',[0 1 -1]; 'A+ C-',[1 0 -1]};
a = exp(1j*2*pi/3);
rows = {}; rowsc = {};
for qa = QA
    for dl = DL
        worst = 0; wpat = ''; warg = 0;
        for p = 1:size(PAT,1)
            for ag = 0:5:355
                q = qa*exp(1j*deg2rad(ag));
                Vabc = [1+q; a^2+a*q; a+a^2*q] .* exp(-1j*deg2rad(dl*PAT{p,2}(:)));
                V1e = (Vabc(1)+a*Vabc(2)+a^2*Vabc(3))/3; V2e = (Vabc(1)+a^2*Vabc(2)+a*Vabc(3))/3;
                d = -0.5*rad2deg(angle((V2e/V1e)/q));
                if abs(d) > abs(worst), worst = d; wpat = PAT{p,1}; warg = ag; end
            end
        end
        rows(end+1,:) = {qa, dl, abs(worst), wpat, warg}; %#ok<SAGROW>
        % clasificador en el peor caso encontrado
        p = find(strcmp(PAT(:,1), wpat)); q = qa*exp(1j*deg2rad(warg));
        tt = (0:159)'/G.fs; e = exp(1j*2*pi*G.f0*tt);
        Vabc = [1+q; a^2+a*q; a+a^2*q] .* exp(-1j*deg2rad(dl*PAT{p,2}(:)));
        V = real(e*Vabc.');
        S0 = r2_sector_scan(V, tt, G.W, G.STRIDE, G.fs, G.f0);
        Vk = V - mean(V,2);
        S1 = r2_sector_scan(Vk, tt, G.W, G.STRIDE, G.fs, G.f0);
        gnd = ismember(S0.label, {'A-G','B-G','C-G','A-B-G','B-C-G','C-A-G'});
        rowsc(end+1,:) = {qa, dl, wpat, warg, thpred(q), cmed(S1.theta), wrap90(cmed(S1.theta)-thpred(q)), worst, 100*mean(gnd)}; %#ok<SAGROW>
    end
end
T8 = cell2table(rows,'VariableNames',{'absq','phase_error_deg','max_abs_dtheta_deg','worst_pattern','worst_argq_deg'});
writetable(T8,'f8_channel_phase.csv');
T8c = cell2table(rowsc,'VariableNames',{'absq','phase_error_deg','pattern','argq_deg','theta_nominal_deg','theta_cls_deg','dtheta_cls_deg','dtheta_pred_deg','pct_windows_ground_label_raw'});
writetable(T8c,'f8_channel_phase_classifier.csv');
e8 = max(abs(T8c.dtheta_cls_deg - T8c.dtheta_pred_deg));
disp(unstack(T8(:,1:3),'max_abs_dtheta_deg','phase_error_deg'));
for qa = [0.42 0.05 0.02]
    k = T8.absq == qa;
    SUM(end+1,:) = {8, 'Channel phase errors', '0 to $2^\circ$ in one or two channels', sprintf('$|V_2/V_1|=%g$', qa), ...
        max(T8.max_abs_dtheta_deg(k)), any(T8.max_abs_dtheta_deg(k) > 30), 'f8_channel_phase.csv', qa ~= 0.02}; %#ok<SAGROW>
end
fprintf('D   max |giro analitico - giro del clasificador| por factor: F1 %.2g, F2 %.2g, F3 %.2g, F4 %.2g, F6 %.2g, F7 %.2g, F8 %.2g grados\n', e1,e2,e3,e4,e6,e7,e8);
ERRF = [e1 e2 e3 e4 e6 e7 e8];

%% ===================== E. FUERA DE LA REGION DE VALIDEZ: EVENTOS ===========
% Senal sintetica en voltios (pico 231 V, como la simulacion del paper):
% 0.30 s de prefalta equilibrada, falta B-C de 0.16 s en el punto de falta
% (Z1 = Z2, cambio instantaneo de fasores, sin transitorio) y 0.06 s de
% posfalta. Decision: cadena de referencia de T03 (r2_decision_chain, valores
% por defecto: puerta adaptativa causal, M = 10) si esta disponible; ademas se
% guardan las metricas ventana a ventana propias (r2_sector_event) como contraste.
fprintf('=== E. Eventos (prefalta 0.30 s, falta 0.16 s, M = 10) ===\n');
USE_T03 = exist(fullfile(RAIZ,'_revision_TPWRD','ronda2','informes','T03_API.md'),'file') == 2 && exist('r2_decision_chain','file') == 2;
fprintf('    cadena de decision de T03 disponible: %d\n', USE_T03);
VPK = 231; ph0 = exp(1j*0.3);
t = (0:round(0.52*G.fs)-1)'/G.fs; t_on = 0.30; t_off = 0.46;
inF = t > t_on & t <= t_off;           % como en la simulacion: la muestra de t_on es la ultima de prefalta
PHe = [45 60 70 php 85 90]; RHe = [0 0.5 1 2 2.2 2.5 3 5 10 30 100 150 180 195 220 300];
PAIRQ = {'B-C',0,'C-A'; 'A-B',-120,'B-C'; 'C-A',120,'A-B'};
rows = {};
for pr = 1:3
    for ph = PHe
        if pr > 1 && ~ismember(ph, [php 85]), continue; end
        for rh = RHe
            if pr > 1 && ~ismember(rh, [0.5 2 3 10]), continue; end
            Z2 = exp(1j*deg2rad(ph)); rot = exp(1j*deg2rad(PAIRQ{pr,2}));
            V1f = (Z2+rh)/(2*Z2+rh); V2f = rot*Z2/(2*Z2+rh);
            V1 = ones(size(t)); V1(inF) = V1f; V2 = zeros(size(t)); V2(inF) = V2f;
            V = VPK*r2_sector_wave(t, ph0*V1, ph0*V2, G.f0);
            S = r2_sector_scan(V, t, G.W, G.STRIDE, G.fs, G.f0);
            E = r2_sector_event(S, PAIRQ{pr,1}, t_on, t_off, G.M);
            q = V2f/V1f;
            r = {PAIRQ{pr,1}, ph, rh, ifelse(abs(ph-php)<1e-9, rh*PAPER.absZ2, NaN), abs(q), -0.5*rad2deg(angle(q/rot)), thpred(q), sector(thpred(q))};
            if USE_T03
                o = r2_decision_chain(V, t, 'truth', PAIRQ{pr,1}, 't_on', t_on, 't_off', t_off);
                e = o.event;
                r = [r, {e.outcome, e.declared, e.decl_class, e.decl_ms, e.first_label_ms, e.correct_decl_ms, e.n_decl_correct, e.n_decl_wrong, ...
                    e.wrong_classes, e.max_run_correct, e.retention_pct, e.full_healthy_pct, e.full_modal_label, e.full_angle_med, e.full_roundness_med, ...
                    e.full_thr_med, e.n_decl_pre}]; %#ok<AGROW>
            else
                r = [r, {'', false, '', NaN, NaN, NaN, NaN, NaN, '', NaN, NaN, NaN, '', NaN, NaN, NaN, NaN}]; %#ok<AGROW>
            end
            r = [r, {E.declared, E.declared_class, E.t_declared_ms, E.t_first_correct_ms, E.n_wrong_secured, E.max_correct_run, ...
                E.pct_full_windows_correct, E.theta_median_deg, E.round_median}]; %#ok<AGROW>
            rows(end+1,:) = r; %#ok<SAGROW>
        end
    end
end
TE = cell2table(rows,'VariableNames',{'true_pair','argZ2_deg','Rf_over_absZ2','Rf_total_ohm_paper_network','absq','shift_pred_deg','theta_pred_deg','sector_pred', ...
    'outcome','declared','decl_class','decl_ms','first_correct_label_ms','correct_decl_ms','n_decl_correct','n_decl_wrong','wrong_classes','max_run_correct', ...
    'retention_pct','full_healthy_pct','full_modal_label','full_angle_med','full_roundness_med','full_thr_med','n_decl_prefault', ...
    'w_declared','w_decl_class','w_decl_ms','w_first_correct_ms','w_n_wrong_secured','w_max_correct_run','w_pct_full_correct','w_theta_median_deg','w_roundness_median'});
writetable(TE,'events_out_of_validity.csv');
kE = strcmp(TE.true_pair,'B-C') & ismember(TE.argZ2_deg,[php 85]);
disp(TE(kE, {'argZ2_deg','Rf_over_absZ2','Rf_total_ohm_paper_network','absq','full_angle_med','outcome','decl_class','decl_ms','first_correct_label_ms','n_decl_wrong','max_run_correct','retention_pct','w_decl_class','w_decl_ms'}));
disp(TE(~strcmp(TE.true_pair,'B-C'), {'true_pair','argZ2_deg','Rf_over_absZ2','absq','full_angle_med','outcome','decl_class','decl_ms'}));
if USE_T03
    fprintf('E   cadena T03 frente a metricas ventana a ventana: clase declarada distinta en %d de %d eventos\n', ...
        sum(~strcmp(TE.decl_class, TE.w_decl_class)), height(TE));
end
% resumen de la franja de seleccion erronea por angulo de Z2 (falta B-C)
rows = {};
for ph = PHe
    k = strcmp(TE.true_pair,'B-C') & TE.argZ2_deg == ph;
    if USE_T03, oc = TE.outcome(k); else, oc = repmat({''},sum(k),1); end
    rh = TE.Rf_over_absZ2(k);
    rows(end+1,:) = {ph, max([rh(strcmp(oc,'correcta')); NaN]), min([rh(strcmp(oc,'erronea')); NaN]), max([rh(strcmp(oc,'erronea')); NaN]), ...
        min([rh(strcmp(oc,'sin_declaracion')); NaN]), sum(strcmp(oc,'correcta')), sum(strcmp(oc,'erronea')), sum(strcmp(oc,'sin_declaracion'))}; %#ok<SAGROW>
end
TEs = cell2table(rows,'VariableNames',{'argZ2_deg','max_rho_correct','min_rho_wrong','max_rho_wrong','min_rho_no_declaration','n_correct','n_wrong','n_no_declaration'});
writetable(TEs,'events_summary.csv'); disp(TEs);

%% ===================== F. TABLA RESUMEN ====================================
TS = cell2table(SUM,'VariableNames',{'factor_id','factor','tested_range','variant','max_abs_shift_deg','crosses_30deg','csv','in_manuscript_table'});
writetable(TS,'summary_8factors.csv'); disp(TS);
cap = ['\\rtwo{Largest displacement of the line-to-line sector angle from its centre for the eight factors of the portability analysis. ' ...
       '$R_f$ is the total resistance between the two faulted phases and $Z_2$ the negative-sequence Th\\''evenin impedance seen from the fault point. ' ...
       'Factors 4, 6 and 7 use a two-terminal model (grid with $\\mathrm{SCR}=5$ and $X/R=10$, line of 0.10~pu, fault at mid-line with $R_f=0.05$~pu) ' ...
       'and the voltages at the converter terminals.}'];
shs = @(x) regexprep(sprintf('$%.1f^\\circ$', x), '^\$(89\.9|9\d\.\d)\^\\circ\$$', 'any');
yns = {'no','yes'};
for vv = 1:2
    if vv == 1, fn = 'T21_tabla.tex'; else, fn = 'T21_tabla_carta.tex'; end
    fid = fopen(fullfile(TEXDIR,fn),'w','n','UTF-8');
    fprintf(fid,'%%%% %s. Generado por experiments/r2_sector_analytic.m desde resultados/T21/summary_8factors.csv. No editar a mano.\n', fn);
    if vv == 1
        fprintf(fid,'%%%% Destino en main_rev2.tex: seccion del ajuste GAC, tras el parrafo de la linea 573 (derivacion de los centros de sector), junto a la figura de T21.\n');
        fprintf(fid,'%%%% Version compacta (una fila por factor). Si no cabe en el presupuesto de paginas (D8), va a la carta de respuesta (R1.2).\n');
        fprintf(fid,'%%%% Requiere booktabs y la macro \\rtwo del manuscrito.\n');
    else
        fprintf(fid,'%%%% Version completa (una fila por variante) para la carta de respuesta (R1.2). Requiere booktabs y \\rtwo.\n');
    end
    fprintf(fid,'\\begin{table}[tb]\n\\centering\n');
    fprintf(fid,['\\caption{' cap '}\n']);
    if vv == 1, fprintf(fid,'\\label{tab:sector_portability}\n'); else, fprintf(fid,'\\label{tab:sector_portability_full}\n'); end
    fprintf(fid,'\\footnotesize\n\\setlength{\\tabcolsep}{3pt}\n');
    if vv == 1
        fprintf(fid,'\\begin{tabular}{@{}>{\\raggedright\\arraybackslash}p{0.17\\columnwidth}>{\\raggedright\\arraybackslash}p{0.43\\columnwidth}>{\\raggedright\\arraybackslash}p{0.18\\columnwidth}>{\\raggedright\\arraybackslash}p{0.13\\columnwidth}@{}}\n\\toprule\n');
        fprintf(fid,'\\textbf{Factor} & \\textbf{Tested range} & \\textbf{Max. shift} & \\textbf{Crosses $\\pm30^\\circ$} \\\\\n\\midrule\n');
        for id = unique(TS.factor_id)'
            k = find(TS.factor_id == id & TS.in_manuscript_table);
            vr = strjoin(TS.variant(k)', '; ');
            fprintf(fid,'%d. %s & %s (%s) & %s & %s \\\\\n', id, TS.factor{k(1)}, TS.tested_range{k(1)}, vr, ...
                strjoin(arrayfun(shs, TS.max_abs_shift_deg(k)', 'UniformOutput', false), '; '), strjoin(yns(TS.crosses_30deg(k)'+1), '; '));
            if id < max(TS.factor_id), fprintf(fid,'\\addlinespace[2pt]\n'); end
        end
    else
        fprintf(fid,'\\begin{tabular}{@{}>{\\raggedright\\arraybackslash}p{0.19\\columnwidth}>{\\raggedright\\arraybackslash}p{0.47\\columnwidth}rc@{}}\n\\toprule\n');
        fprintf(fid,'\\textbf{Factor} & \\textbf{Tested range} & \\textbf{Shift} & \\textbf{Crosses} \\\\\n\\midrule\n');
        last = 0;
        for i = 1:height(TS)
            nm = ''; if TS.factor_id(i) ~= last, nm = sprintf('%d. %s', TS.factor_id(i), TS.factor{i}); if last > 0, fprintf(fid,'\\addlinespace[1pt]\n'); end, end
            last = TS.factor_id(i);
            fprintf(fid,'%s & %s, %s & %s & %s \\\\\n', nm, TS.tested_range{i}, TS.variant{i}, shs(TS.max_abs_shift_deg(i)), yns{TS.crosses_30deg(i)+1});
        end
    end
    fprintf(fid,'\\bottomrule\n\\end{tabular}\n\\end{table}\n');
    fclose(fid);
end

%% ===================== G. CASOS PROPUESTOS PARA LA CAMPANA EMT =============
% Parametros en el formato de r2_emt_run_case.m (T22). OJO: alli Rf es la
% resistencia POR FASE; la total entre fases es 2 Rf. La prediccion usa la ley
% del punto de falta con la Z2th calibrada: Z2th = Zg || (resto de red), donde
% el resto (dos GFM tras su inductancia de red) se dedujo de los registros RL.
% Con una DG desconectada se supone que el resto se duplica; con las dos, que
% queda solo Zg (se desprecia Pex, |Z2| ~ 29 ohm). Son hipotesis a confirmar.
iz2 = find(strcmp(TZ.origen,'resto de red deducido (Z2m sin Zg)'));
Zrest = TZ.R_ohm(iz2) + 1j*TZ.X_ohm(iz2);
wn = 2*pi*50; Rg0 = 0.8929; Lg0 = 16.58e-3;
z2net = @(Rg,Lg,fr,ndg) 1./(1./(Rg + 1j*2*pi*fr*Lg) + (ndg/2)./(real(Zrest) + 1j*imag(Zrest)*fr/50));
C = {};   % grupo, factor, tipo, Rf_total, Lf, red_R, red_L, red_f, Id_ref_pu, dg1_on, dg2_on, nota
for rt = [0.006 0.5 1.0 1.9 2.3 3 5 10 30 100 190 250]
    C(end+1,:) = {'G1','Rf','B-C',rt,0,Rg0,Lg0,50,-0.05,1,1,'red base; barrido de Rf que atraviesa la frontera (2.12 ohm) y el limite de deteccion'}; %#ok<SAGROW>
end
for tp = {'A-B','C-A'}
    for rt = [1.0 3.0]
        C(end+1,:) = {'G1','Rf (otros pares)',tp{1},rt,0,Rg0,Lg0,50,-0.05,1,1,'comprobar el sentido ciclico del error de par'}; %#ok<SAGROW>
    end
end
C(end+1,:) = {'G1','falta inductiva','B-C',0.006,0.5e-3,Rg0,Lg0,50,-0.05,1,1,'caso RL del dataset (giro de signo contrario)'};
for xr = [1 sqrt(3) 3 5.83 10 20]                 % X/R de la fuente con las DG desconectadas (Z2th = Zg), |Zg| = 1 ohm
    Rg = 1/sqrt(1+xr^2); Lg = xr*Rg/wn;
    for mult = [1 3 10]
        C(end+1,:) = {'G2','X/R de la fuente','B-C',mult*1,0,Rg,Lg,50,-0.05,0,0,sprintf('X/R = %.2f, |Zg| = 1 ohm, DG desconectadas',xr)}; %#ok<SAGROW>
    end
end
for zg = [0.25 1 5.2847 20]                        % fortaleza de la fuente, X/R = 5.83, Rf total fijo 2 ohm
    Rg = zg/sqrt(1+5.8335^2); Lg = 5.8335*Rg/wn;
    C(end+1,:) = {'G3','fortaleza de la fuente','B-C',2,0,Rg,Lg,50,-0.05,0,0,sprintf('|Zg| = %.2f ohm, X/R = 5.83, DG desconectadas',zg)}; %#ok<SAGROW>
    C(end+1,:) = {'G3','fortaleza de la fuente','B-C',2,0,Rg,Lg,50,-0.05,1,1,sprintf('|Zg| = %.2f ohm, X/R = 5.83, DG conectadas',zg)}; %#ok<SAGROW>
end
for dg = [1 1; 1 0; 0 1; 0 0]'                     % aportacion remota y control (GFM presentes o no)
    for rt = [1.0 3.0 8.0]
        C(end+1,:) = {'G4','aportacion remota / control','B-C',rt,0,Rg0,Lg0,50,-0.05,dg(1),dg(2),'medir tambien en bornes de DG1 y DG2 (V2, tercer armonico y angulo)'}; %#ok<SAGROW>
    end
end
for id = [-0.15 0 0.15]                            % punto de operacion prefalta
    for rt = [1.0 3.0]
        C(end+1,:) = {'G5','angulo de potencia prefalta','B-C',rt,0,Rg0,Lg0,50,id,1,1,'consigna de corriente activa de Pex'}; %#ok<SAGROW>
    end
end
for fr = [47.5 49 51 51.5]                         % frecuencia de la red
    for rt = [1.0 1.9]
        C(end+1,:) = {'G6','frecuencia','B-C',rt,0,Rg0,Lg0,fr,-0.05,1,1,'clasificador fijado a 50 Hz'}; %#ok<SAGROW>
    end
end
rows = {};
for i = 1:size(C,1)
    ndg = C{i,10} + C{i,11};
    Z2 = z2net(C{i,6}, C{i,7}, C{i,8}, ndg);
    Zf = C{i,4} + 1j*2*pi*C{i,8}*2*C{i,5};          % total entre fases (Lf por fase, dos en serie)
    q = Z2/(Z2+Zf); sh = -0.5*rad2deg(angle(q));
    ctr = containers.Map({'B-C','A-B','C-A'},{165,45,105});
    lab = sector(ctr(C{i,3}) + sh); if (1-abs(q))/(1+abs(q)) > 0.99, lab = 'no detectada'; end
    rows(end+1,:) = {sprintf('%s_%03d',C{i,1},i), C{i,2}, C{i,3}, C{i,4}/2, C{i,4}, C{i,5}, C{i,6}, C{i,7}, C{i,8}, C{i,9}, C{i,10}, C{i,11}, ...
        abs(Z2), rad2deg(angle(Z2)), C{i,4}/abs(Z2), abs(q), sh, mod(ctr(C{i,3})+sh,180), lab, strcmp(lab,C{i,3}), C{i,12}}; %#ok<SAGROW>
end
TG = cell2table(rows,'VariableNames',{'id','factor','tipo','Rf_ohm_por_fase','Rf_total_ohm','Lf_H','red_R_ohm','red_L_H','red_f_Hz','Id_ref_pu','dg1_on','dg2_on', ...
    'pred_absZ2th_ohm','pred_argZ2th_deg','pred_Rf_over_absZ2','pred_absq_Pex','pred_shift_deg','pred_theta_Pex_deg','pred_label_Pex','pred_correct','nota'});
writetable(TG,'emt_cases_proposed.csv');
fprintf('G   %d casos EMT propuestos: prediccion correcta %d, par adyacente %d, no detectada %d\n', height(TG), sum(TG.pred_correct), ...
    sum(~TG.pred_correct & ~strcmp(TG.pred_label_Pex,'no detectada')), sum(strcmp(TG.pred_label_Pex,'no detectada')));

% numeros clave para el informe
K = struct('n_emt',height(TG),'A1max',A1max,'B2err',B2err,'B4mis',B4mis,'B4err',B4err,'ERRF',ERRF,'PAPER',PAPER,'QMIN_FIX',QMIN_FIX,'QMIN_AMB',QMIN_AMB,'USE_T03',USE_T03);
fid = fopen('key_numbers.json','w'); fprintf(fid,'%s',jsonencode(K,'PrettyPrint',true)); fclose(fid);
fprintf('Hecho en %.0f s\n', toc);

%% ===================== funciones locales ====================================
function th = thpred(q)
% angulo del eje mayor en el marco del clasificador (grados, mod 180)
th = mod(165 - 0.5*rad2deg(angle(q)), 180);
end
function d = wrap90(d), d = mod(d+90,180)-90; end
function v = ifelse(c, a, b), if c, v = a; else, v = b; end, end
function r = rho_iso(phi, ps)
% Rf/|Z2| que produce un giro psi = ps (grados) de arg(V2/V1) en el punto de falta
r = NaN(size(phi)); k = phi > ps; r(k) = sind(ps)./sind(phi(k)-ps);
end
function e = errdet(pred, cls, lab)
% error maximo entre giro analitico y giro del clasificador, solo en ventanas
% que llegan a la rama angular (se excluyen las etiquetas sanas por redondez)
k = ~ismember(lab, {'A-B-C','A-B-C-G'});
e = max(abs(mod(pred(k)-cls(k)+90,180)-90));
end
function lab = sector(th)
% misma regla que classify_by_angle de ga_fault_classifier_pure.m
th = mod(th,180);
if th >= 15 && th <= 75, lab = 'A-B'; elseif th >= 75 && th <= 135, lab = 'C-A'; else, lab = 'B-C'; end
end
function m = cmed(th)
% mediana de orientaciones (periodo 180 grados)
mu = 0.5*rad2deg(angle(mean(exp(1j*2*deg2rad(th)))));
m = mod(mu + median(mod(th-mu+90,180)-90), 180);
end
function s = statwin(V1, V2, fsig, G)
% tramo estacionario de 40 ms a fs, clasificador por ventanas (f fijado a G.f0)
% (origen desplazado 0.37 muestras: con la falta franca ideal las tres tensiones
% se anulan a la vez y una muestra exactamente en ese instante da 0/0 al
% normalizar; es un artefacto de la senal ideal, ver informe)
tt = ((0:159)' + 0.37)/G.fs;
V = r2_sector_wave(tt, V1, V2, fsig);
S = r2_sector_scan(V, tt, G.W, G.STRIDE, G.fs, G.f0);
s.theta = cmed(S.theta);
d = mod(S.theta - s.theta + 90, 180) - 90;
s.thmin = mod(s.theta + min(d), 180); s.thmax = mod(s.theta + max(d), 180);
s.round = median(S.round);
[u,~,j] = unique(S.label); c = accumarray(j(:),1); [~,im] = max(c); s.label = u{im};
s.pct = @(x) 100*mean(strcmp(S.label, x));
s.pline = 100*mean(S.isline);
end
function [tc, lc] = clsF(ph, rh, G)
% falta B-C en el punto de falta con Z1 = Z2 = 1/_ph y Rf = rh: angulo y etiqueta del clasificador
Z2 = exp(1j*deg2rad(ph)); s = statwin((Z2+rh)/(2*Z2+rh), Z2/(2*Z2+rh), G.f0, G);
tc = s.theta; lc = s.label;
end
function [tc, lc] = clsQ(V1, V2, G)
s = statwin(V1, V2, G.f0, G); tc = s.theta; lc = s.label;
end
