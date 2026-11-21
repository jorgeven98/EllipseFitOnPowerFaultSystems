%% R2_T15_STAR_FINAL_GATE  (T15; peticiones R1.1 y E1)
% Repite, con el CLASIFICADOR FINAL de la ronda 2 (puerta de secuencia cero G5 con
% umbral 0.05: r2_ga_classifier_param con v0_mode = 'rms_win', kirchhoff_thr = 0.05;
% ver DECISIONES_INTEGRACION.md), las comprobaciones de T11 que dependian de la puerta
% publicada (|x^.k^| < 0.01 en tres muestras):
%   A) etiquetas por ventana con referencia a tierra y con estrella flotante de ramas
%      iguales, para los nueve tipos de falta asimetrica y la monofasica en red aislada;
%   B) regimen sano con una rama de la estrella desajustada: valor del estadistico de
%      la puerta y porcentaje de ventanas con etiqueta de tierra, y desajuste al que el
%      estadistico cruza el umbral 0.05;
%   C) transformador Dy11 y Dy1 entre la falta y la medida (tensiones de fase del lado
%      de estrella y las mismas referidas al lado de la falta).
% Las ondas sinteticas son las de experiments/r2_floating_star_check.m (misma
% definicion de las faltas). Sin aleatoriedad.
% Salidas en _revision_TPWRD/ronda2/resultados/T15/.
clear; clc;
here = fileparts(mfilename('fullpath'));
RAIZ = fullfile(here, '..');
addpath(fullfile(RAIZ, 'src', 'fault_classifiers')); addpath(here);
outdir = fullfile(RAIZ, '_revision_TPWRD', 'ronda2', 'resultados', 'T15');
if ~exist(outdir, 'dir'), mkdir(outdir); end
logf = fopen(fullfile(outdir, 'star_final_gate_log.txt'), 'w');
L = @(varargin) plog(logf, varargin{:});

PF = struct('v0_mode', 'rms_win', 'kirchhoff_thr', 0.05);   % clasificador final
fs = 4000; f = 50; W = 20; STRIDE = 5; w0 = 2*pi*f;
t = (0:1/fs:0.2-1/fs)';
E = [cos(w0*t), cos(w0*t - 2*pi/3), cos(w0*t + 2*pi/3)];
one = [1 1 1]';
L('Clasificador: r2_ga_classifier_param, v0_mode = %s, kirchhoff_thr = %.2f\n', PF.v0_mode, PF.kirchhoff_thr);

cases = {};
ph = 'ABC';
for p = 1:3, cases{end+1} = struct('name', [ph(p) '-G'], 'kind', 'slg', 'p', p); end %#ok<SAGROW>
pairs = [1 2; 2 3; 3 1]; pn = {'A-B','B-C','C-A'};
for q = 1:3, cases{end+1} = struct('name', [pn{q} '-G'], 'kind', 'llg', 'p', pairs(q,:)); end %#ok<SAGROW>
for q = 1:3, cases{end+1} = struct('name', pn{q}, 'kind', 'll', 'p', pairs(q,:)); end %#ok<SAGROW>
cases{end+1} = struct('name', 'A-G', 'kind', 'slg_iso', 'p', 1);
ks = [0.25 0.5 0.75 1.0];
a = exp(2j*pi/3); A = [1 1 1; 1 a a^2; 1 a^2 a]/3;

%% A) etiquetas con ambas referencias
L('\n== A) Etiquetas: referencia a tierra frente a estrella flotante de ramas iguales ==\n');
fid = fopen(fullfile(outdir, 'star_final_gate_labels.csv'), 'w');
fprintf(fid, 'falta,red,k,referencia,V0_sobre_V1_real,estadistico_puerta_mediana,etiqueta_moda,pct_moda,otras_etiquetas,n_ventanas\n');
refs = {'tierra', 'estrella_iguales'};
for ic = 1:numel(cases)
    c = cases{ic};
    for k = ks
        V = make_fault(E, c, k);
        for ir = 1:2
            if ir == 1, Vm = V; else, Vm = ((eye(3) - one*one'/3)*V')'; end
            [lab, st] = classify_windows(Vm, fs, f, W, STRIDE, PF);
            [u, ~, j] = unique(lab); cnt = accumarray(j, 1); [cm, im] = max(cnt);
            others = strjoin(arrayfun(@(i) sprintf('%s:%d', u{i}, cnt(i)), setdiff(1:numel(u), im), 'UniformOutput', false), ' ');
            Sr = A*phasors(V, t, f);
            if strcmp(c.kind, 'slg_iso'), red = 'aislado'; else, red = 'a_tierra'; end
            fprintf(fid, '%s,%s,%.2f,%s,%.4f,%.5f,%s,%.0f,%s,%d\n', c.name, red, k, refs{ir}, ...
                abs(Sr(1))/abs(Sr(2)), median(st), u{im}, 100*cm/numel(lab), others, numel(lab));
            L('  %-6s %-9s k=%.2f %-17s -> %-8s (%3.0f %%) puerta=%.4f %s\n', c.name, red, k, refs{ir}, ...
                u{im}, 100*cm/numel(lab), median(st), others);
        end
    end
end
fclose(fid);

%% B) sano con una rama desajustada
L('\n== B) Regimen sano, estrella con una rama desajustada, Y_a = Y(1 + d) ==\n');
fid = fopen(fullfile(outdir, 'star_final_gate_healthy.csv'), 'w');
fprintf(fid, 'delta_pct,estadistico_puerta_mediana,estadistico_puerta_max,pct_ventanas_con_etiqueta_de_tierra,etiquetas\n');
for d = [0 0.01 0.02 0.03 0.05 0.10 0.14 0.15 0.16 0.17 0.18 0.20 0.30]
    Y = 1 + d*[1;0;0]; w = Y/sum(Y); Vm = ((eye(3) - one*w')*E')';
    [lab, st] = classify_windows(Vm, fs, f, W, STRIDE, PF);
    isg = endsWith(lab, '-G');
    [u, ~, j] = unique(lab); cnt = accumarray(j, 1);
    s = strjoin(arrayfun(@(i) sprintf('%s:%d', u{i}, cnt(i)), 1:numel(u), 'UniformOutput', false), ' ');
    fprintf(fid, '%.1f,%.5f,%.5f,%.1f,%s\n', 100*d, median(st), max(st), 100*mean(isg), s);
    L('    d = %4.1f %%: puerta mediana %.4f, max %.4f; %5.1f %% de ventanas con etiqueta de tierra [%s]\n', ...
        100*d, median(st), max(st), 100*mean(isg), s);
end
fclose(fid);
% desajuste al que el maximo por ventana del estadistico cruza 0.05
dd = linspace(0, 0.30, 3001); gmax = zeros(size(dd)); gmed = gmax;
for i = 1:numel(dd)
    Y = 1 + dd(i)*[1;0;0]; w = Y/sum(Y); Vm = ((eye(3) - one*w')*E')';
    st = gate_only(Vm, fs, f, W, STRIDE);
    gmax(i) = max(st); gmed(i) = median(st);
end
dmax = dd(find(gmax >= 0.05, 1)); dmed = dd(find(gmed >= 0.05, 1));
L('  el estadistico supera 0.05 en alguna ventana desde d = %.2f %% y en la ventana mediana desde d = %.2f %%\n', 100*dmax, 100*dmed);
L('  primer orden: v0_ap = d/3 de la tension de fase -> estadistico d/3 -> 0.05 con d = 15 %%\n');
fid = fopen(fullfile(outdir, 'star_final_gate_cross.csv'), 'w');
fprintf(fid, 'criterio,delta_pct\nalguna_ventana,%.2f\nventana_mediana,%.2f\nprimer_orden,15.00\n', 100*dmax, 100*dmed);
fclose(fid);

%% C) transformador Dy entre la falta y la medida
L('\n== C) Transformador Dy11 y Dy1 entre la falta (lado de red) y la medida (lado de estrella) ==\n');
% Dy11: la tension de fase del secundario adelanta 30 grados a la del primario.
% v_a'' = (v_a - v_b)/sqrt(3) y ciclicas (Dy11); v_a'' = (v_a - v_c)/sqrt(3) y ciclicas (Dy1).
T11m = [1 -1 0; 0 1 -1; -1 0 1]/sqrt(3);
T1m  = [1 0 -1; -1 1 0; 0 -1 1]/sqrt(3);
grupos = {'Dy11', T11m; 'Dy1', T1m};
sel = {cases{1}, cases{4}, cases{7}, cases{10}};          % A-G, A-B-G, A-B, A-G aislada
fid = fopen(fullfile(outdir, 'star_final_gate_dy.csv'), 'w');
fprintf(fid, 'falta,red,k,grupo,senal,etiqueta_moda,pct_moda,otras_etiquetas\n');
for ic = 1:numel(sel)
    c = sel{ic};
    if strcmp(c.kind, 'slg_iso'), red = 'aislado'; else, red = 'a_tierra'; end
    for k = ks
        V = make_fault(E, c, k);
        for ig = 1:2
            Tm = grupos{ig,2};
            Vs = (Tm*V')';                                  % tensiones de fase del lado de estrella
            % referidas al lado de la falta: inversa del grupo sobre el plano de Kirchhoff
            Vr = (pinv(Tm)*Vs')';
            sig = {'estrella_sin_referir', Vs; 'referida_al_lado_de_la_falta', Vr};
            for is = 1:2
                lab = classify_windows(sig{is,2}, fs, f, W, STRIDE, PF);
                [u, ~, j] = unique(lab); cnt = accumarray(j, 1); [cm, im] = max(cnt);
                others = strjoin(arrayfun(@(i) sprintf('%s:%d', u{i}, cnt(i)), setdiff(1:numel(u), im), 'UniformOutput', false), ' ');
                fprintf(fid, '%s,%s,%.2f,%s,%s,%s,%.0f,%s\n', c.name, red, k, grupos{ig,1}, sig{is,1}, u{im}, 100*cm/numel(lab), others);
                if k == 0.5
                    L('  %-6s %-9s k=%.2f %-5s %-29s -> %-8s (%3.0f %%) %s\n', c.name, red, k, grupos{ig,1}, sig{is,1}, u{im}, 100*cm/numel(lab), others);
                end
            end
        end
    end
end
fclose(fid);
fclose(logf);

%% ============================================================ funciones
function plog(fid, varargin)
fprintf(varargin{:}); fprintf(fid, varargin{:});
end

function V = make_fault(E, c, k)
% Misma definicion que experiments/r2_floating_star_check.m
V = E;
switch c.kind
    case 'slg'
        V(:, c.p) = (1 - k)*E(:, c.p);
    case 'llg'
        V(:, c.p(1)) = (1 - k)*E(:, c.p(1));
        V(:, c.p(2)) = (1 - k)*E(:, c.p(2));
    case 'll'
        d = E(:, c.p(1)) - E(:, c.p(2));
        V(:, c.p(1)) = E(:, c.p(1)) - k*d/2;
        V(:, c.p(2)) = E(:, c.p(2)) + k*d/2;
    case 'slg_iso'
        V = E - k*E(:, c.p)*[1 1 1];
end
end

function [labels, stat] = classify_windows(V, fs, f, W, STRIDE, PF)
labels = {}; stat = [];
for s = 1:STRIDE:(size(V,1) - W + 1)
    try
        [ft, p] = r2_ga_classifier_param(V(s:s+W-1, :), fs, f, PF);
        g = p.v0stat;
    catch
        ft = 'Error'; g = NaN;
    end
    labels{end+1} = char(ft); stat(end+1) = g; %#ok<AGROW>
end
end

function stat = gate_only(V, fs, f, W, STRIDE)
% Estadistico G5 de cada ventana: ||X k^|| / ||X||_F
stat = [];
kh = [1 1 1]'/sqrt(3);
for s = 1:STRIDE:(size(V,1) - W + 1)
    X = V(s:s+W-1, :);
    stat(end+1) = norm(X*kh)/norm(X, 'fro'); %#ok<AGROW>
end
end

function Ph = phasors(V, t, f)
X = [cos(2*pi*f*t), sin(2*pi*f*t)];
c = X \ V; Ph = (c(1,:) - 1j*c(2,:)).';
end
