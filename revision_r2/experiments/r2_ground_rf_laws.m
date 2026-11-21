function r2_ground_rf_laws()
% R2_GROUND_RF_LAWS  (T23b) Mecanismo y region de validez de la rama de tierra.
%  A. Ondas sinteticas desde fasores: la etiqueta del clasificador publicado
%     coincide con la que predice la formula de la normal (r2_ground_rf_legacy)
%     y la regla nueva (r2_ground_rf_rule) con el clasificador r2_ga_classifier_ground.
%  B. Mapa de classify_by_bivector en el plano de la inclinacion compleja
%     tau = rho e^{j delta} y trayectorias de la red del paper.
%  C. A-B-G en la red del paper: mapa Rf x Rtierra.
%  D. Ley general de la L-G con el clasificador publicado: Rf de cruce en
%     funcion de arg Z0 y |Z2/Z0|; formula de inclinacion pequena.
%  E. Ley en bornes de un convertidor que regula su tension (DG1): contraste
%     con los casos EMT.
%  F. Barrido de la regla nueva sobre redes (arg Z2, arg Z0 - arg Z2, |Z0/Z2|)
%     y resistencias (Rf, Rtierra).
%  G. Redes de T21 con hipotesis de puesta a tierra: valores en ohmios.
% Salidas en ronda2\resultados\T23b\: synth_check.csv, legacy_tau_map.csv,
% tray_paper.csv, abg_rf_rg_map.csv, abg_rf_rg_frontera.csv,
% legacy_lg_validity.csv, dg1_ley.csv, newrule_scan.csv,
% newrule_scan_resumen.csv, redes_tierra.csv, laws_log.txt
raiz  = fileparts(fileparts(mfilename('fullpath')));
raiz2 = fullfile(raiz,'_revision_TPWRD','ronda2');
res   = fullfile(raiz2,'resultados','T23b');
addpath(fullfile(raiz,'experiments'), fullfile(raiz,'src','fault_classifiers'), fullfile(raiz,'src','utils'));
logf = fullfile(res,'laws_log.txt'); if exist(logf,'file'), delete(logf); end
diary(logf); cu = onCleanup(@() diary('off')); warning('off','all');
fprintf('R2_GROUND_RF_LAWS  %s  MATLAB %s\n', char(datetime('now')), version);
Zt = readtable(fullfile(res,'emt_z0.csv'));
Z2 = Zt.Z2_R + 1j*Zt.Z2_X; Z0 = Zt.Z0_R + 1j*Zt.Z0_X; Z1 = Z2; E = 231; fs = 4000; f = 50;
fprintf('Red del paper: Z2 = %.3f%+.3fj, Z0 = %.3f%+.3fj (|Z0/Z2| = %.2f, arg Z0 - arg Z2 = %.1f grados)\n', real(Z2), imag(Z2), real(Z0), imag(Z0), abs(Z0/Z2), rad2deg(angle(Z0/Z2)));
TIPOS = {'A-G','B-G','C-G','B-C-G','C-A-G','A-B-G'};
P0 = struct('u2_min', 0, 'v2v0_min', 0);

%% A. ondas sinteticas
RF = [0.003 0.01 0.03 0.05 0.1 0.3 0.5 1 1.2 1.5 2 3 5 10 30 100]; RG = [0.01 1 10 100];
rows = {};
for it = 1:6
    for Rg = RG
        for Rf = RF
            if it <= 3 && Rg ~= 0.01, continue, end
            [V, ~, u0, u2] = r2_ground_rf_synth(TIPOS{it}, Z0, Z1, Z2, Rf, Rg, E, 0.1, fs, f);
            L = cell(0,3);
            for s = 1:5:(size(V,1) - 19)
                Vw = V(s:s+19,:);
                l0 = char(ga_fault_classifier_pure(Vw, fs, f));
                l1 = r2_ga_classifier_ground(Vw, fs, f, P0);
                Pg = P0; Pg.metodo = 'geom'; l2 = r2_ga_classifier_ground(Vw, fs, f, Pg);
                L(end+1,:) = {l0, l1, l2}; %#ok<AGROW>
            end
            md = @(c) char(mode(categorical(c)));
            lf = r2_ground_rf_legacy(u0, u2); lr = r2_ground_rf_rule(u0, u2);
            rows(end+1,:) = {TIPOS{it}, Rf, Rg, abs(u0), abs(u2), md(L(:,1)), 100*mean(strcmp(L(:,1), md(L(:,1)))), lf, ...
                md(L(:,2)), 100*mean(strcmp(L(:,2), TIPOS{it})), md(L(:,3)), 100*mean(strcmp(L(:,3), TIPOS{it})), lr}; %#ok<AGROW>
        end
    end
end
A = cell2table(rows, 'VariableNames', {'tipo','Rf_ohm','Rg_ohm','absu0','absu2','legacy_cls','legacy_cls_pct','legacy_formula', ...
    'nueva_conj','nueva_conj_pct_ok','nueva_geom','nueva_geom_pct_ok','nueva_ideal'});
writetable(A, fullfile(res,'synth_check.csv'));
fprintf('\nA. %d casos sinteticos (6 tipos, red del paper)\n', height(A));
fprintf('   etiqueta mayoritaria del clasificador publicado = formula de la normal: %d de %d (minimo de ventanas con la etiqueta mayoritaria %.0f %%)\n', ...
    nnz(strcmp(A.legacy_cls, A.legacy_formula)), height(A), min(A.legacy_cls_pct));
fprintf('   clasificador publicado correcto: %d; regla nueva ideal correcta: %d; r2_ga_classifier_ground conj: %d; geom: %d\n', ...
    nnz(strcmp(A.legacy_cls, A.tipo)), nnz(strcmp(A.nueva_ideal, A.tipo)), nnz(strcmp(A.nueva_conj, A.tipo)), nnz(strcmp(A.nueva_geom, A.tipo)));
for tp = {'A-G','B-G','C-G'}
    k = strcmp(A.tipo, tp{1}); B = A(k,:);
    j = find(~strcmp(B.legacy_cls, B.tipo), 1);
    fprintf('   %s: publicado correcto hasta Rf = %g ohm; con %g ohm declara %s\n', tp{1}, B.Rf_ohm(j-1), B.Rf_ohm(j), B.legacy_cls{j});
end
for tp = {'B-C-G','C-A-G','A-B-G'}
    for Rg = RG
        k = strcmp(A.tipo, tp{1}) & A.Rg_ohm == Rg; B = A(k,:); ok = strcmp(B.legacy_cls, B.tipo);
        j = find(~ok, 1); s = 'correcto en todo el barrido';
        if ~isempty(j) && j > 1, s = sprintf('correcto hasta Rf = %g ohm; con %g ohm declara %s', B.Rf_ohm(j-1), B.Rf_ohm(j), B.legacy_cls{j});
        elseif ~isempty(j), s = sprintf('incorrecto desde Rf = %g ohm (declara %s)', B.Rf_ohm(1), B.legacy_cls{1}); end
        fprintf('   %s, Rtierra = %g ohm: publicado %s\n', tp{1}, Rg, s);
    end
end

%% B. mapa en el plano tau y trayectorias
dl = -180:2:180; rh = logspace(-2, 1.5, 71); rows = {};
for d = dl, for r = rh
    lab = legacy_tau_(r*exp(1j*deg2rad(d)));
    rows(end+1,:) = {d, r, lab}; %#ok<AGROW>
end, end
writetable(cell2table(rows, 'VariableNames', {'delta_deg','rho','etiqueta'}), fullfile(res,'legacy_tau_map.csv'));
% sectores con inclinacion pequena
fprintf('\nB. Mapa de classify_by_bivector: con rho = 0.05 la etiqueta por delta es\n   ');
prev = ''; for d = -180:1:180
    lab = legacy_tau_(0.05*exp(1j*deg2rad(d)));
    if ~strcmp(lab, prev), fprintf('[%d] %s  ', d, lab); prev = lab; end
end, fprintf('\n');
for r = [0.3 0.65 1 1.5 2]
    dd = 0; while dd < 60 && strcmp(legacy_tau_(r*exp(-1j*deg2rad(dd))), 'A-G'), dd = dd + 0.1; end
    fprintf('   rho = %.2f: la etiqueta A-G se mantiene hasta |delta| = %.1f grados\n', r, dd);
end
RFt = logspace(-3, 3, 121); rows = {};
for Rf = RFt
    [u0, u2, ax] = r2_ground_rf_seqlaw('A-G', Z0, Z2, Rf, 0.01);
    rows(end+1,:) = {'A-G', Rf, 0.01, abs(ax.tau), rad2deg(angle(ax.tau)), abs(u0), abs(u2), r2_ground_rf_legacy(u0,u2), r2_ground_rf_rule(u0,u2)}; %#ok<AGROW>
    [u0, u2, ax] = r2_ground_rf_seqlaw('A-B-G', Z0, Z2, Rf, 0.01);
    rows(end+1,:) = {'A-B-G', Rf, 0.01, abs(ax.tau), rad2deg(angle(ax.tau)), abs(u0), abs(u2), r2_ground_rf_legacy(u0,u2), r2_ground_rf_rule(u0,u2)}; %#ok<AGROW>
    [u0, u2, ax] = r2_ground_rf_seqlaw('A-B-G', Z0, Z2, 0.003, Rf);
    rows(end+1,:) = {'A-B-G', 0.003, Rf, abs(ax.tau), rad2deg(angle(ax.tau)), abs(u0), abs(u2), r2_ground_rf_legacy(u0,u2), r2_ground_rf_rule(u0,u2)}; %#ok<AGROW>
end
TR = cell2table(rows, 'VariableNames', {'tipo','Rf_ohm','Rg_ohm','rho','delta_deg','absu0','absu2','legacy','nueva'});
writetable(TR, fullfile(res,'tray_paper.csv'));
k = strcmp(TR.tipo,'A-G'); B = TR(k,:); j = find(~strcmp(B.legacy,'A-G'), 1);
Rc_AG = fzero(@(x) double(strcmp(lg_(Z0, Z2, x), 'A-G')) - 0.5, [B.Rf_ohm(j-1) B.Rf_ohm(j)], optimset('TolX', 1e-4));
[u0c, ~, axc] = r2_ground_rf_seqlaw('A-G', Z0, Z2, Rc_AG, 0);
fprintf('   A-G en el bus de la falta de la red del paper: el clasificador publicado cambia de A-G a %s en Rf = %.3f ohm (3Rf/|Z0| = %.3f; |V0/V1| = %.3f; rho = %.2f, delta = %.1f grados)\n', ...
    B.legacy{j}, Rc_AG, 3*Rc_AG/abs(Z0), abs(u0c), abs(axc.tau), rad2deg(angle(axc.tau)));
xdet = fzero(@(x) abs(r2_ground_rf_seqlaw('A-G', Z0, Z2, x, 0)) - 0.01, [1 1e4]);
fprintf('   |V0/V1| baja al 1 %% en Rf = %.0f ohm (3Rf/|Z0| = %.0f)\n', xdet, 3*xdet/abs(Z0));

%% C. A-B-G en la red del paper: Rf x Rtierra
RFm = logspace(-3, 3, 61); RGm = [0 logspace(-3, 3, 61)]; rows = {}; fr = {};
for Rg = RGm
    okv = false(size(RFm));
    for i = 1:numel(RFm)
        [u0, u2] = r2_ground_rf_seqlaw('A-B-G', Z0, Z2, RFm(i), Rg);
        ll = r2_ground_rf_legacy(u0, u2); ln = r2_ground_rf_rule(u0, u2); okv(i) = strcmp(ll, 'A-B-G');
        rows(end+1,:) = {RFm(i), Rg, abs(u0), abs(u2), ll, ln}; %#ok<AGROW>
    end
    j = find(~okv, 1); if isempty(j), rfmax = Inf; elseif j == 1, rfmax = 0; else, rfmax = sqrt(RFm(j-1)*RFm(j)); end
    j2 = find(okv, 1, 'last');
    fr(end+1,:) = {Rg, rfmax, rfmax/max(Rg, eps), RFm(j2), nnz(okv)}; %#ok<AGROW>
end
writetable(cell2table(rows, 'VariableNames', {'Rf_ohm','Rg_ohm','absu0','absu2','legacy','nueva'}), fullfile(res,'abg_rf_rg_map.csv'));
FR = cell2table(fr, 'VariableNames', {'Rg_ohm','Rf_max_correcto_ohm','Rf_max_sobre_Rg','ultimo_Rf_correcto_ohm','n_Rf_correctos'});
writetable(FR, fullfile(res,'abg_rf_rg_frontera.csv'));
M = cell2table(rows, 'VariableNames', {'Rf_ohm','Rg_ohm','absu0','absu2','legacy','nueva'});
fprintf('\nC. A-B-G, red del paper, %d x %d puntos (Rf por fase x Rtierra)\n', numel(RFm), numel(RGm));
fprintf('   clasificador publicado correcto en %d de %d; regla nueva ideal en %d de %d\n', nnz(strcmp(M.legacy,'A-B-G')), height(M), nnz(strcmp(M.nueva,'A-B-G')), height(M));
for Rg = [0 0.001 0.01 0.1 1 10 100]
    [~, i] = min(abs(FR.Rg_ohm - Rg));
    fprintf('   Rtierra = %-6g ohm: publicado correcto hasta Rf = %.4g ohm por fase (Rf/Rtierra = %.2f)\n', FR.Rg_ohm(i), FR.Rf_max_correcto_ohm(i), FR.Rf_max_sobre_Rg(i));
end
u = unique(M.legacy(~strcmp(M.legacy,'A-B-G'))); fprintf('   etiquetas erroneas del publicado: %s\n', strjoin(u', ', '));
bad = M(~strcmp(M.nueva,'A-B-G'),:);
if ~isempty(bad), fprintf('   regla nueva falla en %d puntos: Rf de %.3g a %.3g, Rtierra de %.3g a %.3g, |u2| de %.3f a %.3f (%s)\n', height(bad), min(bad.Rf_ohm), max(bad.Rf_ohm), min(bad.Rg_ohm), max(bad.Rg_ohm), min(bad.absu2), max(bad.absu2), strjoin(unique(bad.nueva)', ', ')); end

%% D. ley general de la L-G (clasificador publicado), en el punto de falta
PH0 = -90:5:90; MM = [0 0.05 0.1 0.235 1/3 0.5 1 2]; ph2 = 85; rows = {};
X = logspace(-3, 4, 1401);    % x = 3 Rf / |Z0|
for m = MM
    for p0 = PH0
        z0 = exp(1j*deg2rad(p0)); z2 = m*exp(1j*deg2rad(ph2));
        okv = false(size(X)); lab = cell(size(X)); a0 = zeros(size(X));
        for i = 1:numel(X)
            [u0, u2] = r2_ground_rf_seqlaw('A-G', z0, z2, X(i)/3, 0);
            lab{i} = r2_ground_rf_legacy(u0, u2); okv(i) = strcmp(lab{i}, 'A-G'); a0(i) = abs(u0);
        end
        j = find(~okv, 1); idet = find(a0 < 0.01, 1);
        if isempty(j) || (~isempty(idet) && j >= idet), xc = Inf; la = ''; u0c = NaN;
        else, xc = X(j); la = lab{j}; u0c = a0(j); end
        xs = Inf; if abs(p0) > 30, xs = (1 + 2*m*cosd(p0 - ph2))/(2*sind(abs(p0) - 30)); end
        rows(end+1,:) = {m, p0, ph2, xc, la, u0c, xs, X(idet)}; %#ok<AGROW>
    end
end
D = cell2table(rows, 'VariableNames', {'absZ2_sobre_absZ0','argZ0_deg','argZ2_deg','x_cruce_3Rf_sobre_Z0','etiqueta_tras_cruce','absu0_en_cruce','x_formula_inclinacion_pequena','x_deteccion_u0_1pct'});
writetable(D, fullfile(res,'legacy_lg_validity.csv'));
fprintf('\nD. L-G con el clasificador publicado, punto de falta (x = 3Rf/|Z0|, arg Z2 = %d grados)\n', ph2);
fprintf('   arg Z0 :'); fprintf(' %6d', PH0(1:2:end)); fprintf('\n');
for m = [0 0.235 1/3 1]
    k = abs(D.absZ2_sobre_absZ0 - m) < 1e-9; B = D(k,:);
    fprintf('   |Z2/Z0| = %.3f:', m); fprintf(' %6.2f', B.x_cruce_3Rf_sobre_Z0(1:2:end)); fprintf('\n');
end
k = isfinite(D.x_cruce_3Rf_sobre_Z0);
fprintf('   sin cruce hasta el suelo de deteccion: |arg Z0| de %d a %d grados (todas las |Z2/Z0|); etiquetas tras el cruce con Z0 inductiva: %s; capacitiva: %s\n', ...
    min(D.argZ0_deg(~k)), max(D.argZ0_deg(~k)), strjoin(unique(D.etiqueta_tras_cruce(k & D.argZ0_deg > 0))', ', '), strjoin(unique(D.etiqueta_tras_cruce(k & D.argZ0_deg < 0))', ', '));
kk = k & D.absZ2_sobre_absZ0 <= 0.1;
fprintf('   formula de inclinacion pequena frente al clasificador, |Z2/Z0| <= 0.1: cociente de %.2f a %.2f (mediana %.2f)\n', ...
    min(D.x_cruce_3Rf_sobre_Z0(kk)./D.x_formula_inclinacion_pequena(kk)), max(D.x_cruce_3Rf_sobre_Z0(kk)./D.x_formula_inclinacion_pequena(kk)), median(D.x_cruce_3Rf_sobre_Z0(kk)./D.x_formula_inclinacion_pequena(kk)));

%% E. bornes de un convertidor que regula su tension (DG1): giro de arg(V0/V1)
S = readtable(fullfile(res,'emt_seq.csv'), 'TextType','char');
k = strcmp(S.conv,'DG1') & strcmp(S.tipo,'A-G') & ismember(S.conjunto, {'S2a','S2d'}) & S.absu0 > 0.02; B = sortrows(S(k,:), 'Rf_ohm');
ZS = Z0 + Z1 + Z2; rot_m = B.argu0_deg - B.argu0_deg(1); rot_m = mod(rot_m + 180, 360) - 180;
rot_p = rad2deg(angle(ZS + 3*(B.Rf_ohm(1) + 0.01)) - angle(ZS + 3*(B.Rf_ohm + 0.01)));
writetable(table(B.Rf_ohm, B.argu0_deg, rot_m, rot_p, rot_m - rot_p, B.etq_mayor, 'VariableNames', {'Rf_ohm','argu0_DG1_deg','giro_medido_deg','giro_ley_deg','dif_deg','etiqueta'}), fullfile(res,'dg1_ley.csv'));
fprintf('\nE. A-G medida en DG1 (tension regulada): giro de arg(V0/V1) respecto de la falta franca frente a arg(Z0+Z1+Z2) - arg(Z0+Z1+Z2+3Rf)\n');
fprintf('   %d casos, Rf de %g a %g ohm: diferencia maxima %.1f grados (mediana %.1f)\n', height(B), min(B.Rf_ohm), max(B.Rf_ohm), max(abs(rot_m - rot_p)), median(abs(rot_m - rot_p)));
xg = 0.5/sind(rad2deg(angle(ZS)) - 30);
fprintf('   cruce de 30 grados previsto en 3Rf/|Z0+Z1+Z2| = %.3f, Rf = %.2f ohm (EMT: correcta con 1 ohm, C-A-G con 1.5 ohm)\n', xg, xg*abs(ZS)/3);

%% F. barrido de la regla nueva sobre redes
PH2 = [60 70 80 85 90]; DEL = [-90 -75 -60 -45 -30 -15 0 15]; M0 = [0.3 0.5 1 2 3 4.45 6.7 10 20];
RFs = logspace(-3, 3, 25); RGs = [0 logspace(-3, 3, 13)]; rows = {};
for p2 = PH2, for dlt = DEL, for m0 = M0
    p0 = p2 + dlt; if p0 > 90 || p0 < -90, continue, end
    z2 = exp(1j*deg2rad(p2)); z0 = m0*exp(1j*deg2rad(p0));
    n = zeros(2,3);   % familia x (total, nueva ok, publicado ok)
    na = zeros(2,2);  % con la condicion |V2/V0| >= 0.15: activos, activos correctos
    for it = 1:6
        fam = 1 + (it > 3);
        for Rf = RFs
            if fam == 1, RGl = 0; else, RGl = RGs; end
            for Rg = RGl
                [u0, u2] = r2_ground_rf_seqlaw(TIPOS{it}, z0, z2, Rf, Rg);
                okn = strcmp(r2_ground_rf_rule(u0, u2), TIPOS{it}); okl = strcmp(r2_ground_rf_legacy(u0, u2), TIPOS{it});
                n(fam,:) = n(fam,:) + [1, okn, okl];
                if abs(u2) >= 0.15*abs(u0), na(fam,:) = na(fam,:) + [1, okn]; end
            end
        end
    end
    rows(end+1,:) = {p2, dlt, m0, n(1,1), 100*n(1,2)/n(1,1), 100*n(1,3)/n(1,1), n(2,1), 100*n(2,2)/n(2,1), 100*n(2,3)/n(2,1), ...
        100*na(1,1)/n(1,1), 100*na(1,2)/max(na(1,1),1), 100*na(2,1)/n(2,1), 100*na(2,2)/max(na(2,1),1)}; %#ok<AGROW>
end, end, end
Fs = cell2table(rows, 'VariableNames', {'argZ2_deg','argZ0_menos_argZ2_deg','absZ0_sobre_absZ2','n_LG','LG_nueva_pct','LG_publicado_pct','n_LLG','LLG_nueva_pct','LLG_publicado_pct', ...
    'LG_activa_pct','LG_activa_correcta_pct','LLG_activa_pct','LLG_activa_correcta_pct'});
writetable(Fs, fullfile(res,'newrule_scan.csv'));
G = groupsummary(Fs, {'argZ0_menos_argZ2_deg'}, {'min','mean'}, {'LG_nueva_pct','LLG_nueva_pct','LG_publicado_pct','LLG_publicado_pct'});
writetable(G, fullfile(res,'newrule_scan_resumen.csv'));
fprintf('\nF. Regla nueva ideal sobre %d redes (Rf y Rtierra de 0.001 a 1000 veces |Z2|; sin suelos de medida)\n', height(Fs));
disp(G(:, {'argZ0_menos_argZ2_deg','GroupCount','min_LG_nueva_pct','mean_LG_nueva_pct','min_LLG_nueva_pct','mean_LLG_nueva_pct','mean_LG_publicado_pct','mean_LLG_publicado_pct'}));
k = Fs.argZ0_menos_argZ2_deg >= -45 & Fs.argZ0_menos_argZ2_deg <= 15;
fprintf('   con arg Z0 - arg Z2 entre -45 y 15 grados (%d redes): L-G correcta en el %.1f %% como minimo; L-L-G en el %.1f %% como minimo (media %.2f %%)\n', nnz(k), min(Fs.LG_nueva_pct(k)), min(Fs.LLG_nueva_pct(k)), mean(Fs.LLG_nueva_pct(k)));
kb = k & Fs.LLG_nueva_pct < 100; B = Fs(kb,:);
fprintf('   redes de ese grupo con algun fallo en L-L-G: %d; |Z0/Z2| de %g a %g; peor %.1f %% (arg Z2 = %d, dif = %d, |Z0/Z2| = %g)\n', height(B), min(B.absZ0_sobre_absZ2), max(B.absZ0_sobre_absZ2), ...
    min(B.LLG_nueva_pct), B.argZ2_deg(B.LLG_nueva_pct == min(B.LLG_nueva_pct) & (1:height(B))' == find(B.LLG_nueva_pct == min(B.LLG_nueva_pct), 1)), ...
    B.argZ0_menos_argZ2_deg(find(B.LLG_nueva_pct == min(B.LLG_nueva_pct), 1)), B.absZ0_sobre_absZ2(find(B.LLG_nueva_pct == min(B.LLG_nueva_pct), 1)));
k2 = k & Fs.absZ0_sobre_absZ2 <= 6.7;
fprintf('   con ademas |Z0/Z2| <= 6.7 (|V2/V0| >= 0.15 en L-G): L-L-G correcta en el %.1f %% como minimo\n', min(Fs.LLG_nueva_pct(k2)));

%% G. redes de T21 con hipotesis de puesta a tierra
N = readtable(fullfile(raiz2,'resultados','T21','eng_networks.csv'), 'TextType','char', 'Delimiter', ',');
z2n = @(s) N.R2_ohm(contains(N.network, s)) + 1j*N.X2_ohm(contains(N.network, s));
w = 2*pi*50; Vf20 = 20e3/sqrt(3);
ZMT = z2n('barras de subestacion'); ZMT5 = z2n('a 5 km de linea aerea');
ZAT1 = z2n('Sk = 10 GVA'); ZAT2 = z2n('2 GVA + 50 km'); ZMAT2 = z2n('5 GVA + 100 km'); ZBT = z2n('bornes del trafo');
lin0 = @(Zlin, km) real(Zlin) + 0.148*km + 3j*imag(Zlin);     % Z0 de linea: R1 + 3*pi^2*f*1e-4 ohm/km, X0 = 3 X1 (hipotesis)
RED = {
 'Red del paper 0.283 kV, bus de falta', 'rigida (medida)', Z2, Z0, 'Z0 medida en EMT (emt_z0.csv)'
 'BT 0.4 kV, bornes del trafo', 'rigida (TN)', ZBT, ZBT, 'hipotesis Z0 = Z2 (transformador Dyn)'
 'MT 20 kV, barras', 'rigida', ZMT, ZMT, 'hipotesis Z0 = Z2'
 'MT 20 kV, barras', 'resistencia, 500 A', ZMT, 3*Vf20/500 + 1j*imag(ZMT), 'hipotesis: Rn = Vf/500 A, mas X del transformador'
 'MT 20 kV, barras', 'resistencia, 1000 A', ZMT, 3*Vf20/1000 + 1j*imag(ZMT), 'hipotesis: Rn = Vf/1000 A'
 'MT 20 kV, barras', 'aislado, Ic = 30 A', ZMT, -1j*3*Vf20/30, 'hipotesis: corriente capacitiva de 30 A'
 'MT 20 kV, barras', 'aislado, Ic = 100 A', ZMT, -1j*3*Vf20/100, 'hipotesis: corriente capacitiva de 100 A'
 'MT 20 kV, barras', 'Petersen, Ic = 100 A, +5 %, d = 5 %', ZMT, 1/((0.05 - 0.05j)/(3*Vf20/100)), 'hipotesis: sobrecompensacion 5 %, amortiguamiento 5 %'
 'MT 20 kV, barras', 'Petersen, Ic = 100 A, sintonia exacta, d = 5 %', ZMT, 1/(0.05/(3*Vf20/100)), 'hipotesis: sintonia exacta, amortiguamiento 5 %'
 'MT 20 kV, a 5 km de linea aerea', 'resistencia, 500 A', ZMT5, 3*Vf20/500 + 1j*imag(ZMT) + lin0(ZMT5 - ZMT, 5), 'hipotesis: Rn = Vf/500 A; linea con X0 = 3 X1'
 'AT 220 kV, 10 GVA', 'rigida', ZAT1, ZAT1, 'hipotesis Z0 = Z2 (X0/X1 = 1)'
 'AT 220 kV, 10 GVA', 'rigida', ZAT1, 0.5*ZAT1, 'hipotesis X0/X1 = 0.5 (subestacion con muchos neutros a tierra)'
 'AT 220 kV, 2 GVA + 50 km', 'rigida', ZAT2, 1j*24.2 + lin0(ZAT2 - 1j*24.2, 50), 'hipotesis: fuente Z0 = Z1; linea X0 = 3 X1, R0 = R1 + 0.148 ohm/km'
 'MAT 380 kV, 5 GVA + 100 km', 'rigida', ZMAT2, 1j*28.88 + lin0(ZMAT2 - 1j*28.88, 100), 'hipotesis: fuente Z0 = Z1; linea X0 = 3 X1, R0 = R1 + 0.148 ohm/km'
 };
rows = {}; Xr = logspace(-4, 6, 4001);
for i = 1:size(RED,1)
    z2 = RED{i,3}; z0 = RED{i,4};
    okv = false(size(Xr)); lab = cell(size(Xr)); a0 = zeros(size(Xr)); a2 = a0; okn = okv;
    for q = 1:numel(Xr)
        [u0, u2] = r2_ground_rf_seqlaw('A-G', z0, z2, Xr(q), 0);
        lab{q} = r2_ground_rf_legacy(u0, u2); okv(q) = strcmp(lab{q}, 'A-G'); a0(q) = abs(u0); a2(q) = abs(u2);
        okn(q) = strcmp(r2_ground_rf_rule(u0, u2), 'A-G');
    end
    idet = find(a0 < 0.01, 1); if isempty(idet), idet = numel(Xr); end
    j = find(~okv(1:idet), 1);
    if isempty(j), rc = Inf; la = 'ninguna'; else, rc = Xr(j); la = lab{j}; end
    act = a2 >= 0.01 & a2 >= 0.15*a0;                 % regla nueva activa (suelos por defecto)
    aplic = any(act(1:10));                           % activa en falta casi franca
    jn = find(~(act & okn), 1);                       % primera Rf en que deja de estar activa o deja de acertar
    if ~aplic, rn = 0; nota = 'no aplica: |V2/V0| < 0.15'; if any(act & ~okn), nota = [nota '; activa y erronea en parte']; end %#ok<AGROW>
    elseif isempty(jn), rn = Inf; nota = 'activa y correcta en todo el barrido';
    else, rn = Xr(jn); if act(jn), nota = 'deja de acertar'; else, nota = 'deja de estar activa (suelo de |V2/V1|)'; end
    end
    rows(end+1,:) = {RED{i,1}, RED{i,2}, abs(z2), rad2deg(angle(z2)), abs(z0), rad2deg(angle(z0)), abs(z0/z2), rad2deg(angle(z0/z2)), ...
        rc, 3*rc/abs(z0), la, Xr(idet), rn, nota, all(okn(1:idet) | ~act(1:idet)), RED{i,5}}; %#ok<AGROW>
end
NT = cell2table(rows, 'VariableNames', {'red','puesta_a_tierra','absZ2_ohm','argZ2_deg','absZ0_ohm','argZ0_deg','absZ0_sobre_Z2','argZ0Z2_deg', ...
    'Rf_cruce_publicado_ohm','x_cruce_3Rf_sobre_Z0','declara_tras_cruce','Rf_deteccion_u0_1pct_ohm','Rf_max_regla_nueva_ohm','nota_regla_nueva','regla_nueva_sin_error_si_activa','hipotesis_Z0'});
writetable(NT, fullfile(res,'redes_tierra.csv'));
fprintf('\nG. Falta monofasica en el punto de falta, redes de T21 con hipotesis de Z0\n');
for i = 1:height(NT)
    fprintf('   %-36s %-46s |Z0| %8.4g ohm, %6.1f gr, |Z0/Z2| %7.3g | publicado correcto hasta Rf = %8.4g ohm (declara %s) | deteccion hasta %8.4g ohm | regla nueva hasta %8.4g ohm (%s)\n', ...
        NT.red{i}, NT.puesta_a_tierra{i}, NT.absZ0_ohm(i), NT.argZ0_deg(i), NT.absZ0_sobre_Z2(i), NT.Rf_cruce_publicado_ohm(i), NT.declara_tras_cruce{i}, NT.Rf_deteccion_u0_1pct_ohm(i), NT.Rf_max_regla_nueva_ohm(i), NT.nota_regla_nueva{i});
end
end

function lab = lg_(Z0, Z2, Rf)
[u0, u2] = r2_ground_rf_seqlaw('A-G', Z0, Z2, Rf, 0); lab = r2_ground_rf_legacy(u0, u2);
end

function lab = legacy_tau_(tau)
% etiqueta del clasificador publicado para una inclinacion compleja tau dada
% (equivale a u2 = 0 y u0 = -conj(tau)/2)
lab = r2_ground_rf_legacy(-conj(tau)/2, 0);
end
