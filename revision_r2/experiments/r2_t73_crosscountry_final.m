% r2_t73_crosscountry_final.m  (T73) Copia de r2_crosscountry_illustration.m (T71) con UN solo cambio:
% la etiqueta directa sale del clasificador FINAL de la ronda 2 (puerta de secuencia cero G5, umbral 0.05,
% via r2_ga_classifier_param con v0_mode = 'rms_win'), en vez de ga_fault_classifier_pure. Salidas en
% resultados/T73/crosscountry_final/. Sirve para comprobar que la tabla de T71 no cambia con la puerta final.
%% R2_CROSSCOUNTRY_ILLUSTRATION  (T71, peticion R4.1 del revisor 4)
% Ilustracion sintetica: que etiqueta da el clasificador de elipses cuando una
% sola medida local de tension ve a la vez dos faltas a tierra en fases
% distintas y en puntos distintos (cross-country), p. ej. A-G en el alimentador 1
% y B-G en el alimentador 2.
%
% Modelo (fasores, regimen permanente, 1 pu = pico de fase, 50 Hz):
%   - Equivalente de Thevenin pasivo con secuencias (Zs1 = Zs2, Zs0) detras
%     del bus de medida m. Dos alimentadores radiales salen de m (circuitos
%     distintos), con impedancias Zl1 y Zl2 hasta los puntos de falta f1 y f2.
%   - Falta 1: fase A a tierra en f1 con resistencia Rf1. Falta 2: fase B a
%     tierra en f2 con Rf2. Matriz de impedancias de barra Zbus en componentes
%     de fase a partir de las secuencias: Zphase = [ (Z0+2Z1)/3 diag ; (Z0-Z1)/3 resto ].
%   - Se resuelven las dos corrientes de falta y la tension en m por
%     superposicion: V_m = E - Zm1(:,a)*I1a - Zm2(:,b)*I2b.
%   - Se sintetizan ventanas de 20 muestras a 4 kHz (cuarto de ciclo, paso 5
%     muestras, 16 posiciones por ciclo) y se pasan por
%     src/fault_classifiers/ga_fault_classifier_pure.m. Es la ETIQUETA DIRECTA
%     del clasificador, no la cadena con persistencia y puertas de T03. El
%     clasificador no tiene clase 'sano': el caso (sin falta) sale con etiqueta
%     de falta y se marca aparte.
%
% LIMITACIONES de la ilustracion: red pasiva (un IBR con limitacion de
% corriente cambia Zs1, Zs2 y Zs0 vistas desde m), regimen permanente sin
% transitorios ni ruido, solo dos faltas a tierra. Sirve para razonar sobre la
% forma del efecto conjunto, no para cuantificar tasas de acierto.
%
% Salidas (resultados\T71): cross_country_cases.csv, cross_country_grid.csv,
% cross_country_equivalence.csv, cross_country_map.png, cross_country_log.txt

clear; clc;
RAIZ = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(RAIZ, 'src', 'fault_classifiers')); addpath(fullfile(RAIZ, 'experiments'));
OUT = fullfile(RAIZ, '_revision_TPWRD', 'ronda2', 'resultados', 'T73', 'crosscountry_final');
if ~exist(OUT, 'dir'), mkdir(OUT); end
if exist(fullfile(OUT, 'cross_country_log.txt'), 'file'), delete(fullfile(OUT, 'cross_country_log.txt')); end
diary(fullfile(OUT, 'cross_country_log.txt')); diary on;
fprintf('r2_t73_crosscountry_final %s  MATLAB %s\n', datestr(now,31), version);
warning('off', 'all');  % el clasificador avisa si el centro de la elipse no esta en el origen
rng(1);

fs = 4000; f0 = 50; W = 20; STRIDE = 5; w = 2*pi*f0;
a = exp(1j*2*pi/3);
E = [1; a^2; a];                                   % tension prefalta (pico 1 pu)
Zs1 = 0.02 + 0.10j;  Zs0 = 2*Zs1;                  % equivalente visto desde m
Zl1 = 0.05 + 0.15j;  Zl2 = 0.03 + 0.06j;           % alimentadores (angulos distintos)
% (Zl se usa con Z0 = 3 Z1 en lineas)
ZP = @(z0, z1) ((z0 - z1)/3) * ones(3) + z1 * eye(3);   % Z de fase desde secuencias

Zmm = ZP(Zs0, Zs1);
Zm1 = Zmm; Zm2 = Zmm;                                % transferencia m -> f = Zs (radial desde m)
Z11 = ZP(Zs0 + 3*Zl1, Zs1 + Zl1);                    % impedancia propia en f1
Z22 = ZP(Zs0 + 3*Zl2, Zs1 + Zl2);
Z12 = Zmm;                                           % f1 <-> f2 solo comparten Zs

solve_two = @(Rf1, Rf2) two_faults(E, Z11, Z22, Z12, Zm1, Zm2, Rf1, Rf2);

% ---- 1. casos concretos (tabla del manuscrito) ---------------------------
% (Rf1, Rf2) en pu; Inf = esa falta ausente
cases = [
    0.05  Inf;   % A-G sola, severa
    Inf   0.05;  % B-G sola, severa
    0.05  0.05;  % cross-country, igual severidad
    0.02  0.20;  % A 10 veces mas severa (en resistencia)
    0.20  0.02;  % B 10 veces mas severa
    0.01  0.30;  % A 30 veces mas severa
    0.30  0.01;  % B 30 veces mas severa
    0.10  0.10;  % igual severidad, moderada
    0.20  0.20;  % igual severidad, suave (cada falta sola sale mal etiquetada)
    Inf   Inf];  % sin falta (el clasificador no tiene clase 'sano')
nC = size(cases, 1);
Cs = struct([]);
for i = 1:nC
    [V, I1a, I2b] = solve_two(cases(i,1), cases(i,2));
    [lab, frac, labs] = classify_phasor(V, fs, f0, W, STRIDE);
    seq = seqcomp(V);
    Cs(i).Rf1 = cases(i,1); Cs(i).Rf2 = cases(i,2);
    Cs(i).Va = abs(V(1)); Cs(i).Vb = abs(V(2)); Cs(i).Vc = abs(V(3));
    Cs(i).V0 = abs(seq(1)); Cs(i).V1 = abs(seq(2)); Cs(i).V2 = abs(seq(3));
    Cs(i).arg_V2_over_V1_deg = rad2deg(angle(seq(3)/seq(2)));
    Cs(i).label = string(lab); Cs(i).label_fraction = frac; Cs(i).labels_seen = string(strjoin(unique(labs), ' / '));
    Cs(i).note = "";
    if isinf(cases(i,1)) && isinf(cases(i,2)), Cs(i).note = "sin falta: el clasificador no tiene clase sano"; end
    if isinf(cases(i,2)) && ~isinf(cases(i,1)), Cs(i).note = "falta unica A-G"; end
    if isinf(cases(i,1)) && ~isinf(cases(i,2)), Cs(i).note = "falta unica B-G"; end
    fprintf('Rf1=%5.2f Rf2=%5.2f |Va,Vb,Vc|=[%.2f %.2f %.2f] |V0 V1 V2|=[%.3f %.3f %.3f] label=%s (%.0f%% de 16 ventanas) vistas: %s\n', ...
        cases(i,1), cases(i,2), Cs(i).Va, Cs(i).Vb, Cs(i).Vc, Cs(i).V0, Cs(i).V1, Cs(i).V2, lab, 100*frac, Cs(i).labels_seen);
end
Tc = struct2table(Cs); writetable(Tc, fullfile(OUT, 'cross_country_cases.csv'));

% ---- 2. barrido (Rf1, Rf2) ------------------------------------------------
Rgrid = [0.01 0.02 0.03 0.05 0.1 0.2 0.3 0.5 1 2 Inf];   % Inf = esa falta ausente (falta unica)
nG = numel(Rgrid); LAB = strings(nG, nG); FR = zeros(nG, nG);
rows = [];
for i = 1:nG
    for j = 1:nG
        V = solve_two(Rgrid(i), Rgrid(j));
        [lab, frac, labs] = classify_phasor(V, fs, f0, W, STRIDE);
        LAB(i,j) = lab; FR(i,j) = frac;
        seq = seqcomp(V);
        rows = [rows; Rgrid(i), Rgrid(j), abs(V(1)), abs(V(2)), abs(V(3)), abs(seq(1)), abs(seq(2)), abs(seq(3)), frac]; %#ok<AGROW>
        grid_lab(i,j) = {lab}; grid_all(i,j) = {strjoin(unique(labs), '/')}; %#ok<AGROW>
    end
end
Tg = array2table(rows, 'VariableNames', {'Rf1_pu','Rf2_pu','Va_abs','Vb_abs','Vc_abs','V0_abs','V1_abs','V2_abs','label_fraction'});
Tg.label = reshape(LAB.', [], 1);
writetable(Tg, fullfile(OUT, 'cross_country_grid.csv'));
u = unique(LAB(:));
fprintf('\nEtiquetas en el barrido %dx%d (%d casos): ', nG, nG, nG*nG);
for k = 1:numel(u), fprintf('%s=%d  ', u(k), sum(LAB(:)==u(k))); end
fprintf('\nCasos con etiqueta no unanime entre las 16 ventanas: %d\n', sum(FR(:) < 1));
% matriz de etiquetas (filas Rf1, columnas Rf2; Inf = esa falta ausente)
fprintf('Rf1 \\ Rf2  %s\n', strjoin(compose('%6g', Rgrid), ' '));
for i = 1:nG, fprintf('%8g  %s\n', Rgrid(i), strjoin(compose('%6s', LAB(i,:)), ' ')); end
% etiqueta de cada falta unica (borde) y etiqueta del par cuando las dos, solas, se etiquetan bien
fin = isfinite(Rgrid); okA = (LAB(:, end) == "A-G"); okB = (LAB(end, :) == "B-G").';
pair = LAB(fin, fin); mask = okA(fin) & okB(fin).';
fprintf('Falta unica A-G etiquetada A-G en %d de %d niveles; B-G etiquetada B-G en %d de %d\n', ...
    sum(okA(fin)), sum(fin), sum(okB(fin)), sum(fin));
fprintf('Parejas con las dos faltas solas bien etiquetadas: %d. Etiqueta del par: ', sum(mask(:)));
up = unique(pair(mask)); for k = 1:numel(up), fprintf('%s=%d  ', up(k), sum(pair(mask) == up(k))); end; fprintf('\n');
rho = (Rgrid(fin).' ./ Rgrid(fin));   % Rf1/Rf2 (filas Rf1, columnas Rf2)
fprintf('Todos los pares finitos (%d): ', numel(pair));
ua = unique(pair(:)); for k = 1:numel(ua), m3 = pair == ua(k); fprintf('%s=%d (Rf1/Rf2 de %.3g a %.3g)  ', ua(k), sum(m3(:)), min(rho(m3)), max(rho(m3))); end; fprintf('%c', 10);
for lb = reshape(up, 1, [])
    m2 = mask & (pair == lb);
    if any(m2(:)), fprintf('  %s: Rf1/Rf2 entre %.2f y %.2f\n', lb, min(rho(m2)), max(rho(m2))); end
end

% ---- 3. equivalencia: A-B-G en el propio bus m con impedancias desiguales ---------
% Comprobacion INDEPENDIENTE de la resolucion de las dos faltas: se resuelve en m una A-B-G
% con impedancias de falta a tierra por fase Za = Rf1 + Zl1_aa y Zb = Rf2 + Zl2_bb (Zl_aa es la
% impedancia propia de fase del alimentador, (Z0l + 2 Z1l)/3 con Z0l = 3 Z1l) conectadas a las
% fases A y B de m, y V = (I + Zmm*Y)^-1 E. Se compara con la tension del caso cross-country.
Zl1aa = (3*Zl1 + 2*Zl1)/3; Zl2bb = (3*Zl2 + 2*Zl2)/3;
fin_idx = find(isfinite(Rgrid)); eq = zeros(numel(fin_idx)^2, 3); n = 0; maxerr = 0;
for i = fin_idx
    for j = fin_idx
        Vx = solve_two(Rgrid(i), Rgrid(j));
        Za = Rgrid(i) + Zl1aa; Zb = Rgrid(j) + Zl2bb;
        Y = diag([1/Za, 1/Zb, 0]);
        Ve = (eye(3) + Zmm*Y) \ E;
        e = norm(Vx - Ve); maxerr = max(maxerr, e);
        Ia = Ve(1)/Za; Ib = Ve(2)/Zb;
        n = n + 1; eq(n,:) = [e, abs(Za), abs(Zb)];
    end
end
fprintf('\nEquivalencia independiente (alimentadores radiales desde m): max |V_cross - V_ABG_en_m| = %.2e pu\n', maxerr);
writematrix([["err_pu","absZa_eq","absZb_eq"]; eq], fullfile(OUT, 'cross_country_equivalence.csv'));
% comprobacion de las formulas de secuencia del manuscrito (radial, Z_t = Zs):
%   dV0 = -Zs0 (Ia+Ib)/3, dV1 = -Zs1 (Ia + a Ib)/3, dV2 = -Zs1 (Ia + a^2 Ib)/3
errseq = 0;
for i = fin_idx
    for j = fin_idx
        [V, I1a, I2b] = solve_two(Rgrid(i), Rgrid(j));
        dV = seqcomp(V) - seqcomp(E);
        dF = -[Zs0*(I1a + I2b)/3; Zs1*(I1a + a*I2b)/3; Zs1*(I1a + a^2*I2b)/3];
        errseq = max(errseq, norm(dV - dF));
    end
end
fprintf('Formulas de secuencia del manuscrito: max |dV - formula| = %.2e pu%c', errseq, 10);

% ---- 4. mapa de etiquetas ---------------------------------------------------
labs_all = unique(LAB(:)); cmap = [0.20 0.35 0.65; 0.95 0.55 0.15; 0.55 0.75 0.55; 0.80 0.80 0.80; 0.6 0.4 0.7];
fh = figure('Visible', 'off', 'Position', [100 100 560 480]); hold on;
for i = 1:nG
    for j = 1:nG
        k = find(labs_all == LAB(i,j));
        rectangle('Position', [j-0.5 i-0.5 1 1], 'FaceColor', cmap(mod(k-1,size(cmap,1))+1,:), 'EdgeColor', 'w');
        text(j, i, char(LAB(i,j)), 'HorizontalAlignment', 'center', 'FontSize', 7);
    end
end
set(gca, 'XTick', 1:nG, 'XTickLabel', [compose('%g', Rgrid(1:end-1)) {'alone'}], 'YTick', 1:nG, 'YTickLabel', [compose('%g', Rgrid(1:end-1)) {'alone'}], 'YDir', 'normal');
xlabel('R_{f2} (B-G at f2), pu'); ylabel('R_{f1} (A-G at f1), pu'); axis equal tight;
title('Direct classifier label for simultaneous A-G (f1) and B-G (f2)');
exportgraphics(fh, fullfile(OUT, 'cross_country_map.png'), 'Resolution', 200); close(fh);
diary off;

%% ---------------- funciones locales ----------------
function [V, I1a, I2b] = two_faults(E, Z11, Z22, Z12, Zm1, Zm2, Rf1, Rf2)
    % incognitas: I1a (corriente de la falta A-G en f1), I2b (B-G en f2)
    A = zeros(2); b = zeros(2,1);
    a1 = isfinite(Rf1); a2 = isfinite(Rf2);
    if ~a1 && ~a2, V = E; I1a = 0; I2b = 0; return; end
    if a1 && ~a2
        I1a = E(1) / (Z11(1,1) + Rf1); I2b = 0;
    elseif ~a1 && a2
        I2b = E(2) / (Z22(2,2) + Rf2); I1a = 0;
    else
        % E_a - Z11(a,a) I1a - Z12(a,b) I2b = Rf1 I1a
        % E_b - Z12(b,a) I1a - Z22(b,b) I2b = Rf2 I2b
        A = [Z11(1,1) + Rf1, Z12(1,2); Z12(2,1), Z22(2,2) + Rf2]; b = [E(1); E(2)];
        x = A \ b; I1a = x(1); I2b = x(2);
    end
    V = E - Zm1(:,1)*I1a - Zm2(:,2)*I2b;
end
function s = seqcomp(V)
    a = exp(1j*2*pi/3); F = (1/3)*[1 1 1; 1 a a^2; 1 a^2 a];   % [V0; V1; V2]
    s = F * V(:);
end
function [lab, frac, labs] = classify_phasor(V, fs, f0, W, STRIDE)
    n = 80 + W;                       % un ciclo mas una ventana
    t = (0:n-1)'/fs; w = 2*pi*f0;
    Vt = real(exp(1j*w*t) * V(:).');  % [n x 3], pico 1 pu
    labs = {};
    for s = 1:STRIDE:80
        win = Vt(s:s+W-1, :);
        try, ft = r2_ga_classifier_param(win, fs, f0, r2t73_P()); catch, ft = 'Error'; end
        labs{end+1} = char(ft); %#ok<AGROW>
    end
    [uu, ~, ic] = unique(labs); cnt = accumarray(ic, 1);
    [m, k] = max(cnt); lab = uu{k}; frac = m / numel(labs);
end

function P = r2t73_P()
    P = r2_ga_classifier_param('defaults'); P.v0_mode = 'rms_win'; P.kirchhoff_thr = 0.05;
end
