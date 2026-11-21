function r2_fig_rrl()
% R2_FIG_RRL  Datos de la figura R frente a RL de la segunda revision (T94).
%
%   r2_fig_rrl()
%
% Trayectoria de las tensiones fase-tierra en bornes de Pex (bus de la falta,
% nodo 4) para la falta A-B de los registros originales:
%   R  : data/simulation_curated/DataGrid_AB.mat   (Rf 3 mOhm por fase, Lf 0)
%   RL : data/simulation_curated/DataGridRL_AB.mat (Rf 3 mOhm, Lf 0.5 mH por fase)
% (atribucion de T72, informe ronda2/informes/T72.md, tabla del apartado 2).
% Se rota cada muestra con la misma rotacion que el clasificador
% (src/fault_classifiers/ga_fault_classifier_pure.m, pasos 4-5: Rodrigues que
% lleva la normal de Kirchhoff (1,1,1)/sqrt(3) a sigma3) y se guardan las dos
% primeras coordenadas de un ciclo de prefalta [2.90, 2.92) s y de un ciclo de
% falta [3.08, 3.10) s, en p.u. del pico de fase de prefalta.
% Salidas en paper/data/: r2_rrl_{R,RL}_{pre,fault}.csv (columnas x,y) y
% r2_rrl_resumen.csv (eje mayor en grados, b/a y |v0|/|v| de cada caso);
% log en _revision_TPWRD/ronda2/resultados/T94/r2_fig_rrl_log.txt.
here = fileparts(mfilename('fullpath')); RAIZ = fullfile(here, '..');
addpath(fullfile(RAIZ, 'src', 'utils'));
DAT = fullfile(RAIZ, 'data', 'simulation_curated');
OUTD = fullfile(RAIZ, 'paper', 'data');
LOGD = fullfile(RAIZ, '_revision_TPWRD', 'ronda2', 'resultados', 'T94');
if ~exist(LOGD, 'dir'), mkdir(LOGD); end
logf = fullfile(LOGD, 'r2_fig_rrl_log.txt'); if exist(logf, 'file'), delete(logf); end
diary(logf); c0 = onCleanup(@() diary('off'));
fprintf('R2_FIG_RRL  %s  MATLAB %s\n', datestr(now, 31), version); %#ok<TNOW1,DATST>

% rotacion del clasificador: n = (1,1,1)/sqrt(3) -> z
n = [1 1 1]/sqrt(3); z = [0 0 1];
c = dot(n, z); a = cross(n, z); a = a/norm(a); s = sqrt(1 - c^2);
K = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
Rm = eye(3) + s*K + (1-c)*(K*K);

casos = {'R', 'DataGrid_AB.mat'; 'RL', 'DataGridRL_AB.mat'};
res = {};
for k = 1:size(casos, 1)
    d = load(fullfile(DAT, casos{k, 2}));
    [V, t] = extract_converter_data(d, 'Pex');
    V = double(V); t = double(t(:));
    ipre = t >= 2.90 - 1e-9 & t < 2.92 - 1e-9;
    iflt = t >= 3.08 - 1e-9 & t < 3.10 - 1e-9;
    Vpk = median(max(abs(V(ipre, :)), [], 1));        % pico de fase de prefalta
    fprintf('%s (%s): %d y %d muestras; pico de fase de prefalta %.3f V\n', casos{k, 1}, casos{k, 2}, nnz(ipre), nnz(iflt), Vpk);
    for p = {'pre', 'fault'}
        if strcmp(p{1}, 'pre'), idx = ipre; else, idx = iflt; end
        X = V(idx, :)/Vpk;
        Y = (Rm*X')';
        xy = Y(:, 1:2); xy(end+1, :) = xy(1, :); %#ok<AGROW> % cerrar la curva
        writetable(array2table(xy, 'VariableNames', {'x', 'y'}), fullfile(OUTD, sprintf('r2_rrl_%s_%s.csv', casos{k, 1}, p{1})));
        % eje mayor y redondez por covarianza (ciclo completo, elipse centrada)
        P = Y(:, 1:2); C = (P'*P)/size(P, 1);
        [E, L] = eig(C); [lv, ii] = sort(diag(L), 'descend'); E = E(:, ii);
        ang = mod(atan2d(E(2, 1), E(1, 1)), 180);
        ba = sqrt(max(lv(2), 0)/lv(1));
        v0 = norm(X*n')/norm(X, 'fro');                   % estadistico de la puerta (11)
        fprintf('   %-5s eje mayor %.2f deg, b/a %.4f, a %.3f p.u., v0 rel %.2e\n', p{1}, ang, ba, sqrt(2*lv(1)), v0);
        res(end+1, :) = {casos{k, 1}, p{1}, ang, ba, sqrt(2*lv(1)), sqrt(2*lv(2)), v0}; %#ok<AGROW>
    end
end
T = cell2table(res, 'VariableNames', {'caso', 'tramo', 'eje_mayor_deg', 'b_sobre_a', 'a_pu', 'b_pu', 'v0_rel'});
writetable(T, fullfile(OUTD, 'r2_rrl_resumen.csv'));
disp(T);
end
