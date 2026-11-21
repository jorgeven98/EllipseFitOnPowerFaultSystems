function r2_final_tiempos(out_dir)
% R2_FINAL_TIEMPOS  Tiempo de computo de escritorio por ventana (T90, R1.3 y R4.2b).
%
%   r2_final_tiempos            % salidas en ronda2\resultados\T90\tiempos
%
% Mide, en un hilo (maxNumCompThreads(1)) y sin otros procesos MATLAB del autor en
% marcha, sobre las ventanas de un cuarto de ciclo (W = 20, paso 5) de los 24
% registros del paper (t = 2.0 a 3.2 s, 957 ventanas por registro):
%   cadena_final   clasificador final (r2_ga_classifier_param con v0_mode = 'rms_win',
%                  kirchhoff_thr = 0.05; sin GA-FuL) por ventana, mas el coste
%                  amortizado por ventana de la puerta de redondez, la persistencia y
%                  las declaraciones (r2_decision_chain sobre el escaneo, dividido
%                  por el numero de ventanas)
%   cadena_pure    lo mismo con ga_fault_classifier_pure (puerta publicada)
%   ajuste_gac     solo fit_ellipse_gac_improved sobre los puntos girados al plano,
%                  con la misma construccion que ga_fault_classifier_pure (lineas
%                  20 a 82), en las ventanas en que pure llama al ajuste
%   clarke, pe3d   clarke_fault_classifier y pe3d_fault_classifier por ventana
%   gaful          ga_fault_classifier (GA-FuL) en las primeras NGA ventanas
% Cada llamada se cronometra con tic/toc; se descartan las 200 primeras
% (calentamiento del JIT). Se dan mediana, p99 y maximo.
% Salidas: tiempos_por_metodo.csv, tiempos_entorno.csv, tiempos_log.txt y
% ronda2\textos\T90_cifras_tiempo.tex (macros \TchainMed, \TchainMax, \TchainN,
% \TchainPnn, \TfitMed, \TclarkeMed, \TpeMed, \TgaLibMs).
here = fileparts(mfilename('fullpath')); RAIZ = fileparts(here);
addpath(here, fullfile(RAIZ, 'src', 'fault_classifiers'), fullfile(RAIZ, 'src', 'utils'));
R2 = fullfile(RAIZ, '_revision_TPWRD', 'ronda2');
if nargin < 1 || isempty(out_dir), out_dir = fullfile(R2, 'resultados', 'T90', 'tiempos'); end
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
logf = fullfile(out_dir, 'tiempos_log.txt'); if exist(logf, 'file'), delete(logf); end
diary(logf); cu = onCleanup(@() diary('off'));
ws = warning('off', 'all'); cw = onCleanup(@() warning(ws));
nth0 = maxNumCompThreads(1); ct = onCleanup(@() maxNumCompThreads(nth0));
[FO, info0] = r2_final_opts(); P = struct('kirchhoff_thr', FO.kirchhoff_thr, 'v0_mode', FO.v0_mode);
[~, cpu] = system('powershell -NoProfile -Command "(Get-CimInstance Win32_Processor).Name"'); cpu = strtrim(cpu);
fprintf('R2_FINAL_TIEMPOS  %s\nMATLAB %s\nCPU: %s\nHilos de calculo: %d (antes %d)\nOpciones finales: %s\nHuella: %s\n', ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), version, cpu, maxNumCompThreads, nth0, info0.text, info0.md5);
GA = getenv('GAFUL_PATH');   % carpeta del GA-FuL MATLAB Toolbox (opcional)
hayGA = exist(fullfile(GA, 'gafulInit.m'), 'file') == 2; if hayGA, addpath(GA); end
NGA = 3000; NWARM = 200;
[T, info] = r2_sim_cases(); fs = info.fs; f = info.f; W = 20; ST = 5;
% ventanas
Vw = {};
for k = 1:height(T)
    d = load(fullfile(info.dir, T.file{k}));
    for c = 1:2
        [V, t] = extract_converter_data(d, info.converters{c}); V = double(V); t = double(t(:));
        i0 = find(t >= info.t_range(1), 1); i1 = find(t <= info.t_range(2), 1, 'last');
        for s = i0:ST:(i1 - W + 1), Vw{end+1} = V(s:s+W-1, :); end %#ok<AGROW>
    end
end
nW = numel(Vw); fprintf('Ventanas: %d (24 registros)\n', nW);
cls_f = @(X) r2_ga_classifier_param(X, fs, f, P);
cls_p = @(X) ga_fault_classifier_pure(X, fs, f);
tf = tiempo_(cls_f, Vw, NWARM); tp = tiempo_(cls_p, Vw, NWARM);
% puerta, persistencia y declaraciones, amortizado por ventana (sobre los 24 registros)
td = 0; nd = 0;
for k = 1:height(T)
    d = load(fullfile(info.dir, T.file{k}));
    for c = 1:2
        [V, t] = extract_converter_data(d, info.converters{c});
        o = r2_decision_chain(V, t, FO, 't_range', info.t_range);
        for r = 1:5, t0 = tic; r2_decision_chain(o.scan, [], 't_range', info.t_range); td = td + toc(t0); end
        nd = nd + 5*height(o.win);
    end
end
tdec = td/nd; fprintf('Puerta de redondez + persistencia + declaraciones: %.3f us por ventana (amortizado, %d ventanas)\n', 1e6*tdec, nd);
% ajuste GAC solo
XT = {};
for j = 1:nW
    [xt, ok] = puntos_ajuste_(Vw{j}, fs, f); if ok, XT{end+1} = xt; end %#ok<AGROW>
end
tg = tiempo_(@(X) fit_ellipse_gac_improved(X), XT, NWARM);
tc = tiempo_(@(X) clarke_fault_classifier(X), Vw, NWARM);
te = tiempo_(@(X) pe3d_fault_classifier(X), Vw, NWARM);
if hayGA
    ga_fault_classifier(Vw{1}, fs, f);   % inicializa GA-FuL
    tgl = tiempo_(@(X) ga_fault_classifier(X, fs, f), Vw(1:min(nW, NGA + NWARM)), NWARM);
else
    tgl = NaN; fprintf('GA-FuL no encontrado en %s\n', GA);
end
M = {'cadena_final (clasificador final + puerta, persistencia y declaraciones amortizadas)', tf + tdec; ...
     'clasificador_final (r2_ga_classifier_param, G5 0.05)', tf; ...
     'cadena_pure (ga_fault_classifier_pure + puerta, persistencia y declaraciones amortizadas)', tp + tdec; ...
     'clasificador_pure (ga_fault_classifier_pure, puerta publicada)', tp; ...
     'ajuste_gac (fit_ellipse_gac_improved)', tg; 'clarke_fault_classifier', tc; 'pe3d_fault_classifier', te; ...
     'ga_fault_classifier (GA-FuL)', tgl};
rows = {};
for i = 1:size(M, 1)
    x = M{i, 2}; x = x(isfinite(x));
    rows(end+1, :) = {M{i, 1}, numel(x), 1e6*median(x), 1e6*prctile(x, 99), 1e6*max(x), 1e6*mean(x)}; %#ok<AGROW>
end
TT = cell2table(rows, 'VariableNames', {'metodo', 'n_ventanas', 'mediana_us', 'p99_us', 'max_us', 'media_us'});
writetable(TT, fullfile(out_dir, 'tiempos_por_metodo.csv')); disp(TT);
writetable(table(1e6*(tf + tdec), 1e6*(tp + tdec), 'VariableNames', {'cadena_final_us', 'cadena_pure_us'}), fullfile(out_dir, 'tiempos_por_ventana.csv'));
for lim = [125 250 1250]
    fprintf('Ventanas de la cadena por encima de %d us: final %d de %d, pure %d\n', lim, nnz(1e6*(tf + tdec) > lim), numel(tf), nnz(1e6*(tp + tdec) > lim));
end
writetable(cell2table({version, cpu, maxNumCompThreads, NWARM, 1e6*tdec, char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'))}, ...
    'VariableNames', {'matlab', 'cpu', 'hilos', 'descartadas_calentamiento', 'puerta_persistencia_us_por_ventana', 'fecha'}), fullfile(out_dir, 'tiempos_entorno.csv'));
% macros
g = @(n) TT(strcmp(TT.metodo, n), :);
a = g(M{1, 1}); fit = g(M{5, 1}); cl = g(M{6, 1}); pe = g(M{7, 1}); gl = g(M{8, 1});
fid = fopen(fullfile(R2, 'textos', 'T90_cifras_tiempo.tex'), 'w', 'n', 'UTF-8');
fprintf(fid, '%% T90_cifras_tiempo.tex: tiempos de computo de escritorio (resultados/T90/tiempos/tiempos_por_metodo.csv)\n');
fprintf(fid, '%% GENERADO por experiments/r2_final_tiempos.m (%s; %s; un hilo). No editar a mano.\n', version, cpu);
fprintf(fid, '%% Cargar ANTES de T33_manuscrito.tex: alli las macros se definen con \\providecommand.\n');
fprintf(fid, '\\newcommand{\\TchainMed}{%s~$\\mu$s}%% CIFRA: tiempos_por_metodo.csv, mediana_us, cadena_final\n', fmt_(a.mediana_us));
if a.max_us >= 1000
    fprintf(fid, '\\newcommand{\\TchainMax}{%.1f~ms}%% CIFRA: tiempos_por_metodo.csv, max_us/1000, cadena_final\n', a.max_us/1000);
else
    fprintf(fid, '\\newcommand{\\TchainMax}{%s~$\\mu$s}%% CIFRA: tiempos_por_metodo.csv, max_us, cadena_final\n', fmt_(a.max_us));
end
pu = g(M{3, 1});
fprintf(fid, '\\newcommand{\\TchainPureMed}{%s~$\\mu$s}%% CIFRA: tiempos_por_metodo.csv, mediana_us, cadena_pure (clasificador publicado, referencia)\n', fmt_(pu.mediana_us));
fprintf(fid, '\\newcommand{\\TchainOverWin}{%d}%% CIFRA: tiempos_por_ventana.csv, ventanas con cadena_final_us > 1250 (intervalo entre ventanas)\n', nnz(1e6*(tf + tdec) > 1250));
fprintf(fid, '\\newcommand{\\TchainPnn}{%s~$\\mu$s}%% CIFRA: tiempos_por_metodo.csv, p99_us, cadena_final\n', fmt_(a.p99_us));
fprintf(fid, '\\newcommand{\\TchainN}{%s}%% CIFRA: tiempos_por_metodo.csv, n_ventanas, cadena_final\n', miles_(a.n_ventanas));
fprintf(fid, '\\newcommand{\\TfitMed}{%s}%% CIFRA: tiempos_por_metodo.csv, mediana_us, ajuste_gac (microsegundos, sin unidad)\n', fmt_(fit.mediana_us));
fprintf(fid, '\\newcommand{\\TclarkeMed}{%s}%% CIFRA: tiempos_por_metodo.csv, mediana_us, clarke (microsegundos, sin unidad)\n', fmt_(cl.mediana_us));
fprintf(fid, '\\newcommand{\\TpeMed}{%s}%% CIFRA: tiempos_por_metodo.csv, mediana_us, pe3d (microsegundos, sin unidad)\n', fmt_(pe.mediana_us));
fprintf(fid, '\\newcommand{\\TgaLibMs}{%s}%% CIFRA: tiempos_por_metodo.csv, mediana_us/1000, GA-FuL (milisegundos, sin unidad)\n', fmt_(gl.mediana_us/1000));
fclose(fid);
fprintf('Escrito textos\\T90_cifras_tiempo.tex\nFIN r2_final_tiempos\n');
end

function x = tiempo_(fn, C, nwarm)
n = numel(C); x = nan(n, 1);
for j = 1:min(n, nwarm), fn(C{j}); end
for j = 1:n
    t0 = tic; fn(C{j}); x(j) = toc(t0);
end
x = x(nwarm+1:end);
end

function [x_T, ok] = puntos_ajuste_(V_abc, fs, f)
% Misma construccion que ga_fault_classifier_pure (puerta publicada) hasta la
% llamada a fit_ellipse_gac_improved; ok = false si pure no llamaria al ajuste.
x_vectors = V_abc; N = size(x_vectors, 1); q = round(fs/(4*f));
center_idx = round(N/2); idx1 = max(1, center_idx - q); idx2 = min(N, idx1 + q); idx3 = round(abs(idx1 - idx2)/2);
x1 = x_vectors(idx1, :)/norm(x_vectors(idx1, :)); x2 = x_vectors(idx2, :)/norm(x_vectors(idx2, :)); x3 = x_vectors(idx3, :)/norm(x_vectors(idx3, :));
nvec = cross(x1, x2); B = norm(nvec); nK = [1 1 1]/sqrt(3);
if all(abs([dot(x1, nK), dot(x2, nK), dot(x3, nK)]) < 1e-2), n_unit = nK; else, n_unit = nvec/B; end
z = [0 0 1]; c = dot(n_unit, z);
if c > 1 - 1e-12, R = eye(3);
elseif c < -1 + 1e-12, R = diag([1 -1 -1]);
else
    a = cross(n_unit, z); a = a/norm(a); s = sqrt(1 - c^2);
    K = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0]; R = eye(3) + s*K + (1 - c)*(K*K);
end
x_rot = (R*x_vectors')'; x_T = x_rot(:, 1:2); ok = B > 1e-3;
end

function s = fmt_(x)
if x >= 100, s = sprintf('%.0f', x); elseif x >= 10, s = sprintf('%.0f', x); else, s = sprintf('%.1f', x); end
end

function s = miles_(n)
s = sprintf('%d', n); if n >= 10000, s = sprintf('%d\\,%03d', floor(n/1000), mod(n, 1000)); end
end
