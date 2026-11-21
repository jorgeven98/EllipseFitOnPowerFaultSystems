%% R2_FLOATING_STAR_CHECK  (T11, pasos 1-3 y parte del 5; peticiones R1.1 y E1)
% Medida de tension contra una estrella artificial flotante frente a medida
% contra tierra. Comprueba numericamente la derivacion de T11:
%   v_n = sum(Y_k v_k)/sum(Y_k),   v_m = (I - 1 w') v,   w_k = Y_k/sum(Y_j)
% 1) formula frente a una resolucion nodal independiente (MNA) del circuito,
%    con ramas iguales y con desajuste del 1 al 5 %, y con admitancias complejas;
% 2) anulacion de la componente comun con cualquier admitancia;
% 3) fuga de primer orden con desajuste: v0_aparente ~ -(1/3) sum((d_k - dmedia) v_k);
% 4) etiquetas del clasificador (ga_fault_classifier_pure, ventana de 1/4 de ciclo)
%    con referencia a tierra y con estrella flotante, para los nueve tipos de
%    falta asimetrica y la falta monofasica en red de neutro aislado de la ronda 1;
% 5) identidad de reconstruccion: tensiones compuestas + 3V0 (triangulo abierto).
% Salidas en _revision_TPWRD/ronda2/resultados/T11/. Sin aleatoriedad salvo el
% barrido de desajustes (semilla fija).
clear; clc; rng(20261119);
here = fileparts(mfilename('fullpath'));
RAIZ = fullfile(here, '..');
addpath(fullfile(RAIZ, 'src', 'fault_classifiers'));
outdir = fullfile(RAIZ, '_revision_TPWRD', 'ronda2', 'resultados', 'T11');
if ~exist(outdir, 'dir'), mkdir(outdir); end
logf = fopen(fullfile(outdir, 'floating_star_log.txt'), 'w');
P = @(varargin) plog(logf, varargin{:});

fs = 4000; f = 50; W = 20; STRIDE = 5; w0 = 2*pi*f;
t = (0:1/fs:0.2-1/fs)';                         % 10 ciclos en regimen permanente
E = [cos(w0*t), cos(w0*t - 2*pi/3), cos(w0*t + 2*pi/3)];   % 1 pu de pico
one = [1 1 1]';
Pk = eye(3) - one*one'/3;                       % proyector ortogonal sobre el plano de Kirchhoff
khat = one/sqrt(3);

%% ------------------------------------------------------------------ casos
% Faltas sinteticas con V0 conocido. k = severidad (1 = franca).
% Red con neutro a tierra (Z0 = Z1 = Z2): la tension de la fase en falta cae a (1-k).
% Red de neutro aislado (caso de la ronda 1): v = e - k*e_a*1 (solo desplazamiento de neutro).
cases = {};
ph = 'ABC';
for p = 1:3                                     % monofasicas a tierra
    cases{end+1} = struct('name', [ph(p) '-G'], 'kind', 'slg', 'p', p); %#ok<SAGROW>
end
pairs = [1 2; 2 3; 3 1]; pn = {'A-B','B-C','C-A'};
for q = 1:3                                     % bifasicas a tierra y bifasicas
    cases{end+1} = struct('name', [pn{q} '-G'], 'kind', 'llg', 'p', pairs(q,:)); %#ok<SAGROW>
end
for q = 1:3
    cases{end+1} = struct('name', pn{q}, 'kind', 'll', 'p', pairs(q,:)); %#ok<SAGROW>
end
cases{end+1} = struct('name', 'A-G', 'kind', 'slg_iso', 'p', 1);
ks = [0.25 0.5 0.75 1.0];

%% ----------------------------------------- 1) formula frente a resolucion nodal
P('== 1) Formula v_m = (I - 1 w'')v frente a resolucion nodal (MNA) ==\n');
Rnom = 1e6;                                     % ohm, valor nominal de cada rama (irrelevante)
deltas = [0 0.01 0.02 0.03 0.05];
maxrel = 0; maxcm = 0;
Vtest = make_fault(E, cases{4}, 0.6);           % A-B-G, k = 0.6 (contiene V0, V1 y V2)
for d = deltas
    for trial = 1:20
        dk = d*(2*rand(3,1) - 1);               % desajuste de admitancia por rama, |d_k| <= d
        Y = (1/Rnom)*(1 + dk);
        w = Y/sum(Y);
        Vm_formula = (eye(3) - one*w')*Vtest';  % 3 x N
        Vm_mna = mna_star(Vtest', Y);           % resolucion nodal independiente
        rel = norm(Vm_formula - Vm_mna, 'fro')/norm(Vm_mna, 'fro');
        maxrel = max(maxrel, rel);
        % componente comun pura: v = v0(t)*1 debe dar medida nula
        Vcm = one*cos(w0*t' + 0.3);
        maxcm = max(maxcm, max(abs(mna_star(Vcm, Y)), [], 'all'));
    end
end
P('  error relativo maximo formula-MNA (ramas iguales y desajuste hasta 5 %%, 100 casos): %.3g\n', maxrel);
P('  medida maxima ante componente comun pura de 1 pu (100 casos, con desajuste): %.3g pu\n', maxcm);
% admitancias complejas (R en paralelo con C parasita distinta por rama), dominio fasorial a 50 Hz
maxcmc = 0; maxrelc = 0;
for trial = 1:100
    Yc = (1/Rnom)*(1 + 0.05*(2*rand(3,1)-1)) + 1j*w0*1e-9*(1 + 0.5*(2*rand(3,1)-1));
    wc = Yc/sum(Yc);
    Vph = [1; exp(-2j*pi/3); exp(2j*pi/3)]*(0.7 + 0.1j) + 0.4*exp(0.9j)*one;   % fasores con V0 = 0.4 pu
    maxrelc = max(maxrelc, norm((eye(3) - one*wc.')*Vph - mna_star(Vph, Yc))/norm(Vph));
    maxcmc = max(maxcmc, max(abs(mna_star(one*(0.4*exp(0.9j)), Yc))));
end
P('  admitancias complejas (R || C parasita): error relativo maximo %.3g; medida ante V0 puro %.3g pu\n', maxrelc, maxcmc);
% ramas iguales: la medida es la proyeccion ortogonal
Vm_eq = mna_star(Vtest', ones(3,1)/Rnom);
P('  ramas iguales: ||v_m - P_K v||/||v|| = %.3g ; max|suma de las tres medidas| = %.3g pu\n', ...
    norm(Vm_eq - Pk*Vtest', 'fro')/norm(Vtest, 'fro'), max(abs(sum(Vm_eq, 1))));

%% ------------------------------------------------ 2) fuga de primer orden
P('\n== 2) Fuga con desajuste: V0 aparente exacto frente a primer orden ==\n');
P('  convencion: Y_k = Y(1 + d_k); v0_ap = media de las tres medidas = -(w - 1/3)''v\n');
P('  primer orden: v0_ap ~ -(1/3) sum_k (d_k - dmedia) v_k   (con sum v = 0 equivale a -(1/3) sum d_k v_k)\n');
fid = fopen(fullfile(outdir, 'floating_star_leakage.csv'), 'w');
fprintf(fid, 'patron,delta_pct,condicion,v0ap_rms_exacto_pu,v0ap_rms_primer_orden_pu,dif_rel_pct,gate_max_abs_xk,gate_frac_muestras_sobre_0p01\n');
patterns = {'[d 0 0]', [1 0 0]; '[d -d 0]', [1 -1 0]; '[d d 0]', [1 1 0]};
conds = {'sano', E; 'A-B k=0.6 (V0=0)', make_fault(E, cases{7}, 0.6); 'A-G k=0.6 (V0=-0.2 e_a)', make_fault(E, cases{1}, 0.6)};
worst = 0;
for ip = 1:size(patterns,1)
    for d = [0.005 0.01 0.02 0.03 0.05]
        dk = d*patterns{ip,2}';
        Y = 1 + dk; w = Y/sum(Y);
        for ic = 1:size(conds,1)
            V = conds{ic,2}';
            Vm = (eye(3) - one*w')*V;
            v0ex = mean(Vm, 1);
            v0fo = -(1/3)*((dk - mean(dk))'*V);
            difrel = 100*norm(v0ex - v0fo)/norm(v0ex);
            worst = max(worst, difrel);
            g = abs(khat'*(Vm./vecnorm(Vm)));    % metrica de la puerta de Kirchhoff, muestra a muestra
            fprintf(fid, '%s,%.1f,%s,%.6f,%.6f,%.3f,%.5f,%.4f\n', patterns{ip,1}, 100*d, conds{ic,1}, ...
                rms(v0ex), rms(v0fo), difrel, max(g), mean(g > 0.01));
        end
    end
end
fclose(fid);
P('  diferencia relativa maxima exacto-primer orden (desajuste hasta 5 %%): %.2f %%\n', worst);
% desajuste al que el V0 aparente cruza el umbral 0.01 de la puerta en regimen sano, patron [d 0 0]
dd = linspace(0, 0.06, 6001); gmax = zeros(size(dd));
for i = 1:numel(dd)
    Y = 1 + dd(i)*[1;0;0]; w = Y/sum(Y); Vm = (eye(3) - one*w')*E';
    gmax(i) = max(abs(khat'*(Vm./vecnorm(Vm))));
end
dcross = dd(find(gmax > 0.01, 1));
P('  sano, una rama desajustada: max|x^.k^| supera 0.01 a partir de d = %.2f %% (analitico: 0.01*3/sqrt(2) = %.2f %%)\n', 100*dcross, 100*0.01*3/sqrt(2));

%% --------------------------------- 3) clasificador con ambas referencias
P('\n== 3) Etiquetas del clasificador: referencia a tierra frente a estrella flotante ==\n');
refs = {'tierra', [0 0 0]'; 'estrella_iguales', [0 0 0]'; 'estrella_d1pct', 0.01*[1 0 0]'; ...
        'estrella_d2pct', 0.02*[1 0 0]'; 'estrella_d5pct', 0.05*[1 0 0]'};
fid = fopen(fullfile(outdir, 'floating_star_labels.csv'), 'w');
fprintf(fid, 'falta,red,k,referencia,V0_sobre_V1_real,V0_sobre_V1_medido,V2_sobre_V1,theta_medio_deg,redondez_media,etiqueta_moda,pct_moda,otras_etiquetas\n');
a = exp(2j*pi/3); A = [1 1 1; 1 a a^2; 1 a^2 a]/3;
for ic = 1:numel(cases)
    c = cases{ic};
    for k = ks
        V = make_fault(E, c, k);
        for ir = 1:size(refs,1)
            if ir == 1
                Vm = V;
            else
                Y = 1 + refs{ir,2}; w = Y/sum(Y); Vm = ((eye(3) - one*w')*V')';
            end
            [lab, th, rd] = classify_windows(Vm, fs, f, W, STRIDE);
            [u, ~, j] = unique(lab); cnt = accumarray(j, 1);
            [cm, im] = max(cnt);
            others = strjoin(arrayfun(@(i) sprintf('%s:%d', u{i}, cnt(i)), setdiff(1:numel(u), im), 'UniformOutput', false), ' ');
            Sr = A*phasors(V, t, f); Sm = A*phasors(Vm, t, f);
            if strcmp(c.kind, 'slg_iso'), red = 'aislado'; else, red = 'a_tierra'; end
            fprintf(fid, '%s,%s,%.2f,%s,%.4f,%.4f,%.4f,%.1f,%.3f,%s,%.0f,%s\n', c.name, red, k, refs{ir,1}, ...
                abs(Sr(1))/abs(Sr(2)), abs(Sm(1))/abs(Sm(2)), abs(Sm(3))/abs(Sm(2)), circ_mean180(th), mean(rd), ...
                u{im}, 100*cm/numel(lab), others);
            if ismember(ir, [1 2]) && k == 0.5
                P('  %-6s %-9s k=%.2f %-17s -> %-8s (%3.0f %%) theta=%6.1f r=%.2f %s\n', c.name, red, k, refs{ir,1}, ...
                    u{im}, 100*cm/numel(lab), circ_mean180(th), mean(rd), others);
            end
        end
    end
end
fclose(fid);

% sano con desajuste de la estrella: etiquetas en regimen sano
P('\n  Regimen sano, estrella con una rama desajustada (etiquetas por ventana):\n');
fid = fopen(fullfile(outdir, 'floating_star_healthy.csv'), 'w');
fprintf(fid, 'delta_pct,pct_ventanas_con_etiqueta_de_tierra,etiquetas\n');
for d = [0 0.01 0.02 0.021 0.022 0.025 0.03 0.05]
    Y = 1 + d*[1;0;0]; w = Y/sum(Y); Vm = ((eye(3) - one*w')*E')';
    lab = classify_windows(Vm, fs, f, W, STRIDE);
    isg = endsWith(lab, '-G');
    [u, ~, j] = unique(lab); cnt = accumarray(j, 1);
    s = strjoin(arrayfun(@(i) sprintf('%s:%d', u{i}, cnt(i)), 1:numel(u), 'UniformOutput', false), ' ');
    fprintf(fid, '%.1f,%.1f,%s\n', 100*d, 100*mean(isg), s);
    P('    d = %4.1f %%: %5.1f %% de ventanas con etiqueta de tierra  [%s]\n', 100*d, 100*mean(isg), s);
end
fclose(fid);

%% ------------------------- 4) reconstruccion con compuestas + triangulo abierto
P('\n== 4) Reconstruccion fase-tierra = f(compuestas, 3V0) ==\n');
V = make_fault(E, cases{4}, 0.6);               % A-B-G
vab = V(:,1) - V(:,2); vbc = V(:,2) - V(:,3); vca = V(:,3) - V(:,1);
v3v0 = sum(V, 2);                               % salida ideal del triangulo abierto
Vrec = [(vab - vca)/3, (vbc - vab)/3, (vca - vbc)/3] + v3v0/3;
P('  A-B-G k=0.6: error maximo de reconstruccion |v - v_rec| = %.3g pu\n', max(abs(V - Vrec), [], 'all'));
labg = classify_windows(V, fs, f, W, STRIDE); labr = classify_windows(Vrec, fs, f, W, STRIDE);
P('  etiquetas identicas con la senal reconstruida: %d de %d ventanas\n', sum(strcmp(labg, labr)), numel(labg));

P('\nCRITERIOS T11: formula-MNA < 1e-9: %d | componente comun anulada con desajuste: %d | fuga 1er orden dentro del 5 %%: %d\n', ...
    maxrel < 1e-9, max(maxcm, maxcmc) < 1e-9, worst < 5);
fclose(logf);

%% ============================================================ funciones
function plog(fid, varargin)
% Escribe a la vez en consola y en el registro.
fprintf(varargin{:}); fprintf(fid, varargin{:});
end

function V = make_fault(E, c, k)
% Devuelve las tensiones fase-tierra (N x 3) de la falta sintetica.
V = E;
switch c.kind
    case 'slg'                                  % neutro a tierra, Z0 = Z1 = Z2
        V(:, c.p) = (1 - k)*E(:, c.p);
    case 'llg'
        V(:, c.p(1)) = (1 - k)*E(:, c.p(1));
        V(:, c.p(2)) = (1 - k)*E(:, c.p(2));
    case 'll'                                   % las dos fases se acercan a su punto medio; V0 = 0
        d = E(:, c.p(1)) - E(:, c.p(2));
        V(:, c.p(1)) = E(:, c.p(1)) - k*d/2;
        V(:, c.p(2)) = E(:, c.p(2)) + k*d/2;
    case 'slg_iso'                              % neutro aislado: v = e - k e_p (caso de la ronda 1)
        V = E - k*E(:, c.p)*[1 1 1];
end
end

function Vm = mna_star(V, Y)
% Resolucion nodal del circuito: tres fuentes ideales fase-tierra (columnas de V, 3 x N),
% tres ramas de admitancia Y(k) entre cada fase y el nudo flotante n. Se ensambla la
% matriz de admitancias nodal completa (nudos a, b, c, n) por sellado de ramas y se
% resuelve el nudo desconocido. Devuelve v_k - v_n (3 x N).
Yn = zeros(4);
for k = 1:3
    Yn([k 4], [k 4]) = Yn([k 4], [k 4]) + Y(k)*[1 -1; -1 1];
end
vn = -(Yn(4,4) \ (Yn(4,1:3)*V));               % KCL en n: Y_nn v_n + Y_nk v_k = 0
Vm = V - [1;1;1]*vn;
end

function [labels, theta, rnd] = classify_windows(V, fs, f, W, STRIDE)
labels = {}; theta = []; rnd = [];
for s = 1:STRIDE:(size(V,1) - W + 1)
    try
        [ft, p] = ga_fault_classifier_pure(V(s:s+W-1, :), fs, f);
    catch
        ft = 'Error'; p.angle_deg = NaN; p.roundness = NaN;
    end
    labels{end+1} = char(ft); theta(end+1) = p.angle_deg; rnd(end+1) = p.roundness; %#ok<AGROW>
end
end

function Ph = phasors(V, t, f)
X = [cos(2*pi*f*t), sin(2*pi*f*t)];
c = X \ V; Ph = (c(1,:) - 1j*c(2,:)).';
end

function m = circ_mean180(th)
th = th(isfinite(th));
m = mod(angle(mean(exp(2j*deg2rad(th))))/2*180/pi, 180);
end
