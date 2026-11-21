%% R2_DY_TRANSFORMER_CHECK  (T11, paso 5; peticiones R1.1 y E1)
% Que ve el metodo cuando hay un transformador con triangulo entre el punto de
% falta (lado de red, triangulo) y el punto de medida (lado del convertidor, estrella).
% Transformador ideal 1:1 (tension compuesta a compuesta), sin caida interna:
%   Dy11: v_s = T v_p,  T = (I - C)/sqrt(3)   (v_sa = (v_pA - v_pB)/sqrt(3))
%   Dy1 : v_s = T' v_p, T' = (I - C^-1)/sqrt(3)
%   Dy5 = -T (Dy11 con polaridad invertida), Dy7 = -T'.
% En componentes simetricas: V0 se anula, V1 gira +30 h' y V2 gira -30 h', de modo
% que arg(V2/V1) cambia 60 grados y el eje mayor de la elipse gira 30 grados.
% Comprueba con ondas sinteticas y el clasificador de referencia:
%  (a) faltas a tierra del lado de red (neutro a tierra y neutro aislado),
%  (b) faltas bifasicas con angulo residual arg(V2/V1) no nulo (impedancia de falta),
%  (c) el efecto de referir la medida al lado de red (tensiones compuestas/sqrt(3)).
% Salidas en _revision_TPWRD/ronda2/resultados/T11/. Determinista.
clear; clc;
here = fileparts(mfilename('fullpath'));
RAIZ = fullfile(here, '..');
addpath(fullfile(RAIZ, 'src', 'fault_classifiers'));
outdir = fullfile(RAIZ, '_revision_TPWRD', 'ronda2', 'resultados', 'T11');
if ~exist(outdir, 'dir'), mkdir(outdir); end
logf = fopen(fullfile(outdir, 'dy_transformer_log.txt'), 'w');
P = @(varargin) plog(logf, varargin{:});

fs = 4000; f = 50; W = 20; STRIDE = 5; w0 = 2*pi*f;
t = (0:1/fs:0.2-1/fs)';
a = exp(2j*pi/3); A = [1 1 1; 1 a a^2; 1 a^2 a]/3; Ainv = [1 1 1; 1 a^2 a; 1 a a^2];
C = [0 1 0; 0 0 1; 1 0 0];                      % (C v)_a = v_b
T11 = (eye(3) - C)/sqrt(3); T1 = (eye(3) - C')/sqrt(3);
% {nombre, matriz al secundario, matriz que refiere el secundario al lado de red}
trafos = {'sin_trafo', eye(3), eye(3); 'Dy11', T11, T1; 'Dy1', T1, T11; 'Dy5', -T11, -T1; 'Dy7', -T1, -T11};
synth = @(S) real((Ainv*S(:))*exp(1j*w0*t')).'; % S = [V0; V1; V2] referidos a la fase a -> N x 3

% comprobacion de las ganancias de secuencia de cada matriz
P('== Ganancias de secuencia de cada conexion (fasor de salida / fasor de entrada) ==\n');
for i = 1:size(trafos,1)
    G = A*trafos{i,2}*Ainv;                     % diagonal en componentes simetricas
    P('  %-9s  V0: %.3f   V1: %.3f ang %+6.1f   V2: %.3f ang %+6.1f   giro de theta = %+5.1f grados (mod 180)\n', trafos{i,1}, ...
        abs(G(1,1)), abs(G(2,2)), angle(G(2,2))*180/pi, abs(G(3,3)), angle(G(3,3))*180/pi, ...
        wrap90(-0.5*(angle(G(3,3)) - angle(G(2,2)))*180/pi));
end

fid = fopen(fullfile(outdir, 'dy_transformer_labels.csv'), 'w');
fprintf(fid, 'falta_lado_red,red,severidad,rho_V2V1,phi_residual_deg,conexion,senal,V0_sobre_V1,V2_sobre_V1,theta_teorico_deg,theta_ajuste_deg,margen_a_frontera_deg,redondez,etiqueta_moda,pct_moda,otras\n');

%% (a) faltas a tierra del lado de red
P('\n== (a) Faltas a tierra en el lado de red, vistas en el lado de estrella ==\n');
E = synth([0; 1; 0]);
gcases = {'A-G','a_tierra',@(k) setcol(E, 1, 1-k); 'B-G','a_tierra',@(k) setcol(E, 2, 1-k); ...
          'C-G','a_tierra',@(k) setcol(E, 3, 1-k); ...
          'A-B-G','a_tierra',@(k) setcol(E, [1 2], 1-k); 'B-C-G','a_tierra',@(k) setcol(E, [2 3], 1-k); ...
          'C-A-G','a_tierra',@(k) setcol(E, [3 1], 1-k); ...
          'A-G','aislado',@(k) E - k*E(:,1)*[1 1 1]};
for ic = 1:size(gcases,1)
    for k = [0.25 0.5 0.75 1.0]
        Vp = gcases{ic,3}(k);
        for it = 1:size(trafos,1)
            sigs = {'secundario', (trafos{it,2}*Vp')'};
            if it > 1, sigs(2,:) = {'referida_a_red', (trafos{it,3}*trafos{it,2}*Vp')'}; end
            for is = 1:size(sigs,1)
                r = analyse(sigs{is,2}, t, fs, f, W, STRIDE, A);
                fprintf(fid, '%s,%s,%.2f,,,%s,%s,%.4f,%.4f,%.1f,%.1f,%.1f,%.3f,%s,%.0f,%s\n', gcases{ic,1}, gcases{ic,2}, k, ...
                    trafos{it,1}, sigs{is,1}, r.v0, r.v2, r.th_theory, r.th_fit, r.margin, r.rnd, r.mode, r.pct, r.others);
                if k == 0.5 && ismember(it, [1 2 3]) && is == 1
                    P('  %-6s %-9s k=0.50  %-9s -> %-7s (%3.0f %%)  V0/V1=%.3f V2/V1=%.3f theta=%6.1f (teorico %6.1f) %s\n', gcases{ic,1}, ...
                        gcases{ic,2}, trafos{it,1}, r.mode, r.pct, r.v0, r.v2, r.th_fit, r.th_theory, r.others);
                end
            end
        end
    end
end

%% (b) faltas bifasicas con angulo residual
P('\n== (b) Faltas bifasicas en el lado de red con angulo residual phi = arg(V2/V1) - valor ideal ==\n');
P('  (V2/V1 = Z2/(Z2 + Zf): phi > 0 con falta resistiva y red inductiva; phi = 0 en la falta franca)\n');
pairs = {'A-B', a^-1; 'B-C', 1; 'C-A', a};      % giro de referencia de V2/V1 para cada par
phis = [-20 -10 -5 -2 0 2 5 10 20];
summary = cell(numel(phis), 3, size(trafos,1), 2);
for ip = 1:3
    for rho = [1.0 0.6 0.3]
        for iph = 1:numel(phis)
            S = [0; 1; rho*exp(1j*deg2rad(phis(iph)))*pairs{ip,2}]/(1 + rho);
            Vp = synth(S);
            for it = 1:size(trafos,1)
                sigs = {'secundario', (trafos{it,2}*Vp')'};
                if it > 1, sigs(2,:) = {'referida_a_red', (trafos{it,3}*trafos{it,2}*Vp')'}; end
                for is = 1:size(sigs,1)
                    r = analyse(sigs{is,2}, t, fs, f, W, STRIDE, A);
                    fprintf(fid, '%s,,,%.2f,%d,%s,%s,%.4f,%.4f,%.1f,%.1f,%.1f,%.3f,%s,%.0f,%s\n', pairs{ip,1}, rho, phis(iph), ...
                        trafos{it,1}, sigs{is,1}, r.v0, r.v2, r.th_theory, r.th_fit, r.margin, r.rnd, r.mode, r.pct, r.others);
                    if rho == 0.6
                        summary{iph, ip, it, is} = sprintf('%s(%.0f)', r.mode, r.th_fit);
                    end
                end
            end
        end
    end
end
fclose(fid);
for ip = 1:3
    P('\n  Falta %s en el lado de red, rho = |V2/V1| = 0.6: etiqueta de moda (theta ajustado)\n', pairs{ip,1});
    P('  %6s | %-12s | %-12s %-12s | %-12s %-12s | %-12s %-12s\n', 'phi', 'sin trafo', 'Dy11 sec.', 'Dy11 ref.', 'Dy1 sec.', 'Dy1 ref.', 'Dy5 sec.', 'Dy7 sec.');
    for iph = 1:numel(phis)
        P('  %+6d | %-12s | %-12s %-12s | %-12s %-12s | %-12s %-12s\n', phis(iph), summary{iph,ip,1,1}, summary{iph,ip,2,1}, ...
            summary{iph,ip,2,2}, summary{iph,ip,3,1}, summary{iph,ip,3,2}, summary{iph,ip,4,1}, summary{iph,ip,5,1});
    end
end
P('\nLectura: tras un Dy el eje mayor de una bifasica franca cae en una frontera de sector (15, 75 o 135 grados);\n');
P('la etiqueta depende del signo del angulo residual y del indice horario. Referida al lado de red recupera el centro.\n');
fclose(logf);

%% ============================================================ funciones
function plog(fid, varargin)
fprintf(varargin{:}); fprintf(fid, varargin{:});
end

function V = setcol(E, cols, g)
V = E; V(:, cols) = g*E(:, cols);
end

function y = wrap90(x)
y = mod(x + 90, 180) - 90;
end

function r = analyse(V, t, fs, f, W, STRIDE, A)
% Fasores de secuencia, angulo teorico theta = 165 - arg(V2/V1)/2 (mod 180) y etiquetas por ventana.
X = [cos(2*pi*f*t), sin(2*pi*f*t)]; c = X \ V; S = A*(c(1,:) - 1j*c(2,:)).';
r.v0 = abs(S(1))/abs(S(2)); r.v2 = abs(S(3))/abs(S(2));
r.th_theory = mod(165 - 0.5*angle(S(3)/S(2))*180/pi, 180);
if r.v2 < 1e-6, r.th_theory = NaN; end
labels = {}; th = []; rd = [];
for s = 1:STRIDE:(size(V,1) - W + 1)
    [ft, p] = ga_fault_classifier_pure(V(s:s+W-1, :), fs, f);
    labels{end+1} = char(ft); th(end+1) = p.angle_deg; rd(end+1) = p.roundness; %#ok<AGROW>
end
[u, ~, j] = unique(labels); cnt = accumarray(j, 1); [cm, im] = max(cnt);
r.mode = u{im}; r.pct = 100*cm/numel(labels);
r.others = strjoin(arrayfun(@(i) sprintf('%s:%d', u{i}, cnt(i)), setdiff(1:numel(u), im), 'UniformOutput', false), ' ');
r.th_fit = mod(angle(mean(exp(2j*deg2rad(th))))/2*180/pi, 180);
r.rnd = mean(rd);
b = [15 75 135];                                % fronteras de sector
d = abs(mod(r.th_fit - b + 90, 180) - 90); r.margin = min(d);
end
