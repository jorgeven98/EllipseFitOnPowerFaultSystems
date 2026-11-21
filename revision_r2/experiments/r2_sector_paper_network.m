%% R2_SECTOR_PAPER_NETWORK  (T21, R1.2 / X4)
% Contraste de la ley angular con el dataset simulado del paper y calibracion
% de la impedancia de secuencia inversa vista desde el bus de falta (bus de Pex).
%
% Para cada registro L-L (9 archivos) y cada convertidor (Pex, DG1, DG2):
%   - fasores de secuencia V1, V2 (e I1, I2) por minimos cuadrados sobre dos
%     ciclos de falta estacionaria, con v(t) = Re{V e^{jwt}};
%   - angulo predicho  theta = 165 - (1/2) arg(V2/V1)  (mod 180);
%   - angulo medido por ga_fault_classifier_pure (mediana de las ventanas de
%     un cuarto de ciclo dentro del tramo de falta).
% Con los registros RL (falta con rama serie R+L conocida) se invierte
%   q = Z2th/(Z2th+Zf)  ->  Z2th = q Zf/(1-q)
% tomando Zf como impedancia TOTAL entre las dos fases en falta:
%   Zf = 2 (Ron + R_FaultLoad + j w L_FaultLoad).
%
% Ejecutar desde ronda2\resultados\T21 con R2024b:
%   matlab -batch "addpath('<RAIZ>\experiments'); r2_sector_paper_network"
% Salidas (carpeta actual): paper_network_seq.csv, paper_network_z2.csv

RAIZ = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(RAIZ,'src','fault_classifiers'));
addpath(fullfile(RAIZ,'src','utils'));
warning('off','all');
fs = 4000; f = 50; W = 20; STRIDE = 5; w = 2*pi*f;
a = exp(1j*2*pi/3);

% archivo, par en falta, Ron por fase, R Fault Load por fase, L Fault Load por fase
% (valores de INFORME_proc_sim.md, seccion 3.2: DOC del bloque, DED por archivo)
F = { ...
 'DataGridRFAB.mat','A-B',0.002,0,0; 'DataGridRFBC.mat','B-C',0.002,0,0; 'DataGridRFCA.mat','C-A',0.002,0,0; ...
 'DataGridRL_AB.mat','A-B',0.002,0.001,0.5e-3; 'DataGridRL_BC.mat','B-C',0.002,0.001,0.5e-3; 'DataGridRL_AC.mat','C-A',0.002,0.001,0.5e-3; ...
 'DataGrid_AB.mat','A-B',0.002,0.001,0; 'DataGrid_BC.mat','B-C',0.002,0.001,0; 'DataGrid_AC.mat','C-A',0.002,0.001,0};
CONV = {'Pex','DG1','DG2'};
TPH = [3.06 3.10];      % dos ciclos de falta estacionaria (valido para las tres familias)
TCL = [3.045 3.11];     % ventanas del clasificador completamente dentro de la falta
qbolt = containers.Map({'B-C','A-B','C-A'}, {0, -120, 120});   % arg(V2/V1) de la falta franca
cent  = containers.Map({'B-C','A-B','C-A'}, {165, 45, 105});

rows = {};
for i = 1:size(F,1)
    S = load(fullfile(RAIZ,'data','simulation_curated',F{i,1}));
    Zf = 2*(F{i,3} + F{i,4} + 1j*w*F{i,5});          % total entre fases
    for c = 1:numel(CONV)
        cv = CONV{c};
        [V,t] = extract_converter_data(S, cv);
        I = S.GridConvs.(cv).Io.Data;
        npk = max(abs(V(t>=2.4 & t<2.8,:)),[],'all');
        m = t>=TPH(1) & t<TPH(2);
        Vp = phasors(V(m,:), t(m), w);  Ip = phasors(I(m,:), t(m), w);
        V1 = (Vp(1)+a*Vp(2)+a^2*Vp(3))/3;  V2 = (Vp(1)+a^2*Vp(2)+a*Vp(3))/3;  V0 = sum(Vp)/3;
        I1 = (Ip(1)+a*Ip(2)+a^2*Ip(3))/3;  I2 = (Ip(1)+a^2*Ip(2)+a*Ip(3))/3;
        q = V2/V1;
        th_pred = mod(165 - 0.5*rad2deg(angle(q)), 180);
        psi = wrap180(rad2deg(angle(q)) - qbolt(F{i,2}));       % giro respecto a la falta franca
        % clasificador real, ventana a ventana
        i0 = find(t>=TCL(1),1); i1 = find(t<TCL(2),1,'last');
        th = []; rn = []; lab = {};
        for s = i0:STRIDE:(i1-W+1)
            [ft,P] = ga_fault_classifier_pure(V(s:s+W-1,:)/npk, fs, f);
            lab{end+1} = char(ft); th(end+1) = P.angle_deg; rn(end+1) = P.roundness; %#ok<SAGROW>
        end
        ok = strcmp(lab, F{i,2});
        nong = ~ismember(lab, {'A-G','B-G','C-G','A-B-G','B-C-G','C-A-G'});
        th_med = cmedian180(th(nong));
        z2c = -V2/I2;     % impedancia inversa aparente del convertidor con Io saliente (ver signo de P prefalta)
        % potencia activa prefalta con el convenio de Io (positiva = Io sale del convertidor hacia la red)
        mp = t>=2.60 & t<2.64; Vq = phasors(V(mp,:), t(mp), w); Iq = phasors(I(mp,:), t(mp), w);
        Ppre = 0.5*real(sum(Vq.*conj(Iq)));
        % armonicos del vector espacial: s(t) = sum_k c_k e^{j k w t}, k = 1,-1,3,-3,2,-2,0
        sv = (2/3)*(V(m,1) + a*V(m,2) + a^2*V(m,3)); tm = t(m);
        KH = [1 -1 3 -3 2 -2 0]; ck = exp(1j*w*tm*KH)\sv;
        V1s = ck(1); V2s = conj(ck(2)); V3s = ck(3);
        qef = (V2s + V3s*exp(-1j*2*angle(V1s)))/V1s;      % firma efectiva de la orientacion (primer orden)
        th_eff = mod(165 - 0.5*rad2deg(angle(qef)), 180);
        % distancia minima de las ventanas a la frontera de sector mas cercana
        dctr = mod(th(nong) - cent(F{i,2}) + 90, 180) - 90;
        rows(end+1,:) = {F{i,1}, F{i,2}, cv, real(Zf), imag(Zf), abs(V1)/npk, abs(q), rad2deg(angle(q)), psi, ...
            abs(V0)/abs(V1), th_pred, th_med, wrap180(th_med-th_pred), th_pred-cent(F{i,2}), median(rn), ...
            (1-abs(q))/(1+abs(q)), 100*mean(ok), numel(lab), abs(I1), abs(I2), abs(z2c), rad2deg(angle(z2c)), 1e-3*Ppre, ...
            abs(V3s)/abs(V1s), abs(ck(5))/abs(V1s), abs(ck(7))/abs(V1s), abs(qef), th_eff, wrap180(th_med-th_eff), ...
            min(dctr), max(dctr), 30-max(abs(dctr))}; %#ok<SAGROW>
    end
end
T = cell2table(rows,'VariableNames',{'file','pair','conv','Zf_total_R_ohm','Zf_total_X_ohm','V1_pu','absq','argq_deg', ...
    'psi_deg','V0_over_V1','theta_pred_deg','theta_cls_median_deg','theta_err_deg','theta_shift_from_centre_deg', ...
    'roundness_cls_median','roundness_from_q','pct_windows_correct','n_windows','I1_A','I2_A','Z2conv_abs_ohm','Z2conv_arg_deg','P_prefault_kW', ...
    'V3pos_over_V1','h2_over_V1','dc_over_V1','absq_eff','theta_pred_eff_deg','theta_err_eff_deg', ...
    'win_dev_from_centre_min_deg','win_dev_from_centre_max_deg','min_margin_to_boundary_deg'});
writetable(T,'paper_network_seq.csv');
disp(T(:,{'file','conv','absq','psi_deg','theta_pred_deg','theta_cls_median_deg','theta_err_deg','V3pos_over_V1','absq_eff','theta_pred_eff_deg','theta_err_eff_deg','min_margin_to_boundary_deg','pct_windows_correct','Z2conv_abs_ohm','Z2conv_arg_deg','P_prefault_kW'}));
% resumen por convertidor
G = groupsummary(T, 'conv', {'min','max'}, {'absq','psi_deg','theta_shift_from_centre_deg','theta_err_deg','theta_err_eff_deg','V3pos_over_V1','absq_eff','min_margin_to_boundary_deg','pct_windows_correct','Z2conv_abs_ohm','Z2conv_arg_deg','P_prefault_kW'});
writetable(G,'paper_network_summary.csv'); disp(rows2vars(G));

%% Calibracion de Z2th en el bus de falta con los registros RL (medida en Pex = bus de falta)
Rg = 0.8929; Lg = 16.58e-3;                 % Three-Phase Source (SpecifyImpedance = off: se usan R y L)
Zg = Rg + 1j*w*Lg;
Lvsi = 0.085179365*5*(231/40)/314.15; Rvsi = 0.005*231/40;   % VSIdesign.m
Zdg = Rvsi + 1j*w*Lvsi;                     % impedancia serie del filtro de cada GFM
rz = {};
for i = 4:6
    k = strcmp(T.file,F{i,1}) & strcmp(T.conv,'Pex');
    qn = T.absq(k)*exp(1j*deg2rad(T.psi_deg(k)));      % q referido a la falta franca del par
    Zf = T.Zf_total_R_ohm(k) + 1j*T.Zf_total_X_ohm(k);
    Z2 = qn*Zf/(1-qn);
    rz(end+1,:) = {F{i,1}, abs(qn), T.psi_deg(k), real(Z2), imag(Z2), abs(Z2), rad2deg(angle(Z2)), imag(Z2)/real(Z2)}; %#ok<SAGROW>
end
Z2m = mean(cellfun(@(r,x) r+1j*x, rz(:,4), rz(:,5)));
rz(end+1,:) = {'media RL (Pex)', NaN, NaN, real(Z2m), imag(Z2m), abs(Z2m), rad2deg(angle(Z2m)), imag(Z2m)/real(Z2m)};
rz(end+1,:) = {'Zg (R=0.8929, L=16.58 mH)', NaN, NaN, real(Zg), imag(Zg), abs(Zg), rad2deg(angle(Zg)), imag(Zg)/real(Zg)};
Zp = 1/(1/Zg + 2/Zdg);
rz(end+1,:) = {'Zg || Zfiltro_DG1 || Zfiltro_DG2', NaN, NaN, real(Zp), imag(Zp), abs(Zp), rad2deg(angle(Zp)), imag(Zp)/real(Zp)};
% impedancia equivalente del resto de la red (todo menos Zg) deducida de la medida
Zrest = 1/(1/Z2m - 1/Zg);
rz(end+1,:) = {'resto de red deducido (Z2m sin Zg)', NaN, NaN, real(Zrest), imag(Zrest), abs(Zrest), rad2deg(angle(Zrest)), imag(Zrest)/real(Zrest)};
TZ = cell2table(rz,'VariableNames',{'origen','absq','psi_deg','R_ohm','X_ohm','absZ_ohm','argZ_deg','X_over_R'});
writetable(TZ,'paper_network_z2.csv'); disp(TZ);
phi = angle(Z2m);
if rad2deg(phi) > 60
    rc = sind(60)/sin(phi-pi/3);
    fprintf('Red del paper: |Z2th| = %.3f ohm, arg = %.2f deg, X/R = %.2f\n', abs(Z2m), rad2deg(phi), tan(phi));
    fprintf('  cruce si Rf_total > %.2f |Z2| = %.2f ohm; en el cruce |V2/V1| = %.3f\n', rc, rc*abs(Z2m), sin(phi-pi/3)/sin(phi));
end

function P = phasors(X, t, w)
% ajuste LS de x(t) = Re{P e^{jwt}} + dc, por columna
A = [cos(w*t), -sin(w*t), ones(size(t))];
c = A\X;  P = (c(1,:) + 1j*c(2,:)).';
end
function d = wrap180(d), d = mod(d+180,360)-180; end
function m = cmedian180(th)
% mediana de orientaciones (periodo 180): se centra en la media circular
if isempty(th), m = NaN; return; end
mu = 0.5*rad2deg(angle(mean(exp(1j*2*deg2rad(th)))));
d = mod(th-mu+90,180)-90;  m = mod(mu+median(d),180);
end
