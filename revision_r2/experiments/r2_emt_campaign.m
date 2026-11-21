function varargout = r2_emt_campaign(accion, conjunto, w, nw)
% R2_EMT_CAMPAIGN  Campana EMT de la ronda 2 (T23): define los casos, los
% simula con r2_emt_run_case (varios procesos en paralelo, cada uno con su
% carpeta de trabajo) y construye el indice ronda2\emt\indice_casos.csv.
%
%   C = r2_emt_campaign('definir', conjunto)     tabla de casos del conjunto
%   r2_emt_campaign('ejecutar', conjunto, w, nw) trabajador w de nw (salta los hechos)
%   r2_emt_campaign('indice')                    reescribe indice_casos.csv
%
% Conjuntos (ver ronda2\informes\T23.md):
%   S2a  barrido de resistencia de falta (A-G, A-B, A-B-G), incepcion en t = 3.000 s
%   S2d  refinamiento del barrido de Rf alrededor de los suelos vistos en S2a
%   S2b  angulo de incidencia 0/45/90/135 grados respecto a la tension del bucle
%        en falta, para Rf franca, media y alta
%   S2c  A-B-G con la resistencia a tierra barrida (Rf de fase franca)
%   S1a  L-L (A-B, B-C, C-A), un factor cada vez: Rf, X/R y fortaleza de la red,
%        punto de operacion, aportacion remota, frecuencia
%   S1b  interacciones Rf x X/R y Rf x fortaleza
%   S1c  frecuencia de red con desviaciones pequenas (+-0.1 y +-0.2 Hz)
%   S1t  los 69 casos propuestos por el analisis teorico T21 (resultados\T21\emt_cases_proposed.csv)
%   XC   prueba de la segunda falta simultanea en el nudo de union (cross-country, para T71)
%
% Convenciones del indice:
%   Rf       resistencia por fase en falta (ohm): Ron del interruptor + R de la rama serie.
%            En L-L la resistencia total entre fases es 2*Rf; en L-G es Rf + Rtierra.
%   Rtierra  resistencia a tierra de la falta (ohm). No confundir con red_R (equivalente de red).
%   ang_inc  angulo de incidencia pedido (grados): fase de la tension del bucle en falta en el
%            bus de la falta (v_a en A-G; v_a - v_b en A-B y A-B-G) con v = V*sin(ang);
%            0 = paso por cero ascendente, 90 = cresta. NaN = no controlado (t = 3.000 s).
%
% Ronda 2 TPWRD, tarea T23.

raiz  = fileparts(fileparts(mfilename('fullpath')));
raiz2 = fullfile(raiz,'_revision_TPWRD','ronda2');
dcasos = fullfile(raiz2,'emt','casos');
dlogs  = fullfile(raiz2,'emt','logs');
addpath(fullfile(raiz,'experiments'), fullfile(raiz,'src','utils'));
todos = {'S2a','S2d','S2b','S2c','S1a','S1b','S1c','S1t','XC'};

switch lower(accion)
    case 'definir'
        varargout{1} = definir(conjunto, raiz2);

    case 'ejecutar'
        C = definir(conjunto, raiz2);
        if nargin < 4, w = 1; nw = 1; end
        wd = fullfile(raiz2,'emt','work',sprintf('w%d',w));
        flog = fullfile(dlogs, sprintf('%s_w%d.csv', conjunto, w));
        for k = w:nw:numel(C)
            c = C(k); arch = fullfile(dcasos, [c.id '.mat']);
            if exist(arch,'file'), continue, end
            p = rmfield(c, {'conjunto','factor','ang_inc','nota'});
            p.dir_trabajo = wd; p.archivo = arch;
            p.meta = struct('conjunto',c.conjunto,'factor',c.factor,'ang_inc_pedido',c.ang_inc,'nota',c.nota);
            try
                [G, info] = r2_emt_run_case(p);
                v = [G.Pex.Vabc.Data(:); G.DG1.Vac.Data(:); G.DG2.Vac.Data(:)];
                if any(~isfinite(v))
                    est = 'no_finito'; msg = 'valores no finitos en las tensiones';
                else
                    est = 'ok'; msg = '';
                end
                linea(flog, c.id, est, info.t_sim, msg);
            catch ME
                msg = regexprep(ME.message, '[\r\n,]+', ' ');
                linea(flog, c.id, 'error', NaN, msg(1:min(end,400)));
                % tras un error, recargar el modelo limpio para el caso siguiente
                try, bdclose('all'); catch, end
            end
        end
        try, bdclose('all'); catch, end

    case 'indice'
        filas = {};
        cab = {'id','conjunto','factor','estado','tipo','tipo2','Rf_ohm','Rtierra_ohm','Lf_H','R_bucle_ohm', ...
            't_falta_s','dur_falta_s','t_stop_s','ang_inc_pedido_deg','Id_ref_pu','dg1_dP_W', ...
            'red_R_ohm','red_L_H','red_Z_ohm','red_XR','red_f_Hz','dg1_on','dg2_on','modelo', ...
            't_sim_s','archivo','nota','causa_fallo','Ron_ohm'};
        L = leer_logs(dlogs);
        for q = 1:numel(todos)
            C = definir(todos{q}, raiz2);
            for k = 1:numel(C)
                c = C(k); arch = fullfile(dcasos, [c.id '.mat']);
                est = 'pendiente'; ts = NaN; causa = '';
                if isKey(L, c.id), e = L(c.id); est = e{1}; ts = e{2}; causa = e{3}; end
                if exist(arch,'file') && ~strcmp(est,'no_finito'), est = 'ok'; end
                ron = min(0.002, c.Rf/2);
                if exist(arch,'file')      % valores efectivos guardados con la onda
                    try
                        S = load(arch,'caso'); ts = S.caso.t_sim_s; c.t_falta = S.caso.t_falta; ron = S.caso.Ron;
                    catch
                    end
                end
                t = upper(strrep(c.tipo,'-','')); nf = sum(ismember(t,'ABC'));
                if any(t=='G') && nf == 1, Rb = c.Rf + c.Rtierra;       % bucle fase-tierra
                elseif nf >= 2, Rb = 2*c.Rf;                             % bucle entre fases
                else, Rb = NaN; end
                X = 2*pi*c.red_f*c.red_L;
                filas(end+1,:) = {c.id, c.conjunto, c.factor, est, c.tipo, c.tipo2, c.Rf, c.Rtierra, c.Lf, Rb, ...
                    round(c.t_falta,5), c.dur_falta, round(c.t_falta + c.dur_falta + 0.14,5), c.ang_inc, c.Id_ref_pu, c.dg1_dP, ...
                    c.red_R, c.red_L, hypot(c.red_R,X), X/c.red_R, c.red_f, double(c.dg1_on), double(c.dg2_on), c.modelo, ...
                    ts, ['casos\' c.id '.mat'], c.nota, causa, ron}; %#ok<AGROW>
            end
        end
        T = cell2table(filas, 'VariableNames', cab);
        writetable(T, fullfile(raiz2,'emt','indice_casos.csv'));
        varargout{1} = T;

    otherwise
        error('Accion desconocida: %s', accion);
end
end

% =========================================================================
function C = definir(conjunto, raiz2)
% Caso base = estado guardado de GFM.slx (el que reprodujo DataGrid_ABG.mat).
base = struct('id','', 'conjunto',conjunto, 'factor','', 'modelo','GFM', 'tipo','A-B', ...
    'Rf',0.003, 'Rtierra',0.01, 'Lf',0, 't_falta',3.0, 'dur_falta',0.16, 'ang_inc',NaN, ...
    'Id_ref_pu',-0.15, 'dg1_dP',800, 'red_R',0.8929, 'red_L',16.58e-3, 'red_f',50, ...
    'dg1_on',true, 'dg2_on',true, 'tipo2','NONE', 'Rf2',0.003, 'Rtierra2',0.01, 'nota','');
C = base([]);
Zg = hypot(0.8929, 2*pi*50*16.58e-3);          % 5.285 ohm
xr0 = 2*pi*50*16.58e-3/0.8929;                 % 5.83
redxr = @(xr, esc) deal(esc*Zg/sqrt(1+xr^2), esc*Zg*xr/sqrt(1+xr^2)/(2*pi*50));   % [R, L] para X/R y escala de |Z|

switch conjunto
    case 'S2a'
        Rfs = [0.003 0.01 0.1 1 3 10 30 100 300 1000];
        for tp = {'A-G','A-B','A-B-G'}
            for r = Rfs
                c = base; c.tipo = tp{1}; c.Rf = r; c.factor = 'Rf';
                C(end+1) = c; %#ok<AGROW>
            end
        end

    case 'S2d'
        % Refinamiento tras el piloto S2a: la etiqueta deja de ser correcta entre 1 y 3 ohm
        % (A-G, A-B; A-B-G en DG1) y entre 0.01 y 0.1 ohm (A-B-G en Pex).
        for tp = {'A-G','A-B','A-B-G'}
            for r = [0.03 0.05 0.3 1.5 2]
                c = base; c.tipo = tp{1}; c.Rf = r; c.factor = 'Rf';
                C(end+1) = c; %#ok<AGROW>
            end
        end

    case 'S2b'
        % Rf franca, media y alta por tipo (elegidas con el piloto S2a, ver T23.md)
        niveles = struct('AG',[0.003 1 30], 'AB',[0.003 1 30], 'ABG',[0.003 1 30]);
        cal = fullfile(raiz2,'resultados','T22','humo_ABG_GFM.mat');   % misma prefalta que el caso base
        for tp = {'A-G','A-B','A-B-G'}
            phi0 = fase_bucle(cal, tp{1}, 3.0, 50);
            for r = niveles.(strrep(tp{1},'-',''))
                for ang = [0 45 90 135]
                    c = base; c.tipo = tp{1}; c.Rf = r; c.factor = 'incidencia'; c.ang_inc = ang;
                    c.t_falta = 3.0 + mod(ang - phi0, 360)/360/50;
                    c.nota = sprintf('fase del bucle en t=3.000 s: %.2f grados', phi0);
                    C(end+1) = c; %#ok<AGROW>
                end
            end
        end

    case 'S2c'
        for rt = [1 10 100 1000]
            c = base; c.tipo = 'A-B-G'; c.Rf = 0.003; c.Rtierra = rt; c.factor = 'Rtierra';
            c.nota = 'A-B-G con fases francas y tierra resistiva: tiende a A-B';
            C(end+1) = c; %#ok<AGROW>
        end

    case 'S1a'
        Rf0 = 0.71;     % Rf de referencia de los factores de red: 2*Rf/|Z2| = 1.5 con |Z2| = 0.95 ohm (piloto S2a)
        for tp = {'A-B','B-C','C-A'}
            b = base; b.tipo = tp{1}; b.Rf = Rf0;
            % Rf: 2*Rf/|Z2| de 0 a ~5 (|Z2| medido en el piloto, ver T23.md)
            for r = [0.003 0.24 0.47 0.71 0.95 1.19 1.42 1.9 2.37]     % 2*Rf/|Z2| = 0, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5
                c = b; c.Rf = r; c.factor = 'Rf'; C(end+1) = c; %#ok<AGROW>
            end
            % X/R de la red a |Zg| constante
            for xr = [1 3 10 20]
                c = b; [c.red_R, c.red_L] = redxr(xr, 1); c.factor = 'red_XR'; C(end+1) = c; %#ok<AGROW>
            end
            % fortaleza de la red: escala de |Zg| a X/R constante
            for esc = [0.25 0.5 2 4]
                c = b; [c.red_R, c.red_L] = redxr(xr0, esc); c.factor = 'red_fortaleza'; C(end+1) = c; %#ok<AGROW>
            end
            % punto de operacion prefalta: consigna de Pex y escalon de DG1
            for id = [-0.05 0 0.05 0.15]
                c = b; c.Id_ref_pu = id; c.factor = 'Id_ref_Pex'; C(end+1) = c; %#ok<AGROW>
            end
            for dp = [0 1600]
                c = b; c.dg1_dP = dp; c.factor = 'dP_DG1'; C(end+1) = c; %#ok<AGROW>
            end
            % aportacion remota
            c = b; c.dg1_on = false; c.factor = 'remota'; c.nota = 'DG1 desconectada'; C(end+1) = c; %#ok<AGROW>
            c = b; c.dg2_on = false; c.factor = 'remota'; c.nota = 'DG2 desconectada'; C(end+1) = c; %#ok<AGROW>
            c = b; c.dg1_on = false; c.dg2_on = false; c.factor = 'remota'; c.nota = 'DG1 y DG2 desconectadas'; C(end+1) = c; %#ok<AGROW>
            % inductancia serie de la falta (impedancia de falta no resistiva)
            for lf = [0.5e-3 3e-3]
                c = b; c.Lf = lf; c.factor = 'Lf'; C(end+1) = c; %#ok<AGROW>
            end
            % frecuencia de la red
            for fr = [49.5 50.5]
                c = b; c.red_f = fr; c.factor = 'frecuencia'; C(end+1) = c; %#ok<AGROW>
            end
        end

    case 'S1b'
        for tp = {'A-B','B-C','C-A'}
            for r = [0.24 0.71 1.42]      % 2*Rf/|Z2| = 0.5, 1.5, 3
                for xr = [1 20]
                    c = base; c.tipo = tp{1}; c.Rf = r; [c.red_R, c.red_L] = redxr(xr, 1);
                    c.factor = 'Rf_x_XR'; C(end+1) = c; %#ok<AGROW>
                end
                for esc = [0.25 4]
                    c = base; c.tipo = tp{1}; c.Rf = r; [c.red_R, c.red_L] = redxr(xr0, esc);
                    c.factor = 'Rf_x_fortaleza'; C(end+1) = c; %#ok<AGROW>
                end
            end
        end

    case 'S1t'
        % Lista final de T21, con su prediccion por caso (incluye los que atraviesan la frontera de sector).
        % Se lee una copia congelada (la del momento de simular), para que indice y ondas no dejen de
        % corresponderse si T21 regenera su CSV.
        ft21 = fullfile(raiz2,'emt','T21_emt_cases_congelado.csv');
        if ~exist(ft21,'file'), copyfile(fullfile(raiz2,'resultados','T21','emt_cases_proposed.csv'), ft21); end
        if exist(ft21,'file')
            P = readtable(ft21, 'TextType','char');
            for k = 1:height(P)
                c = base; c.tipo = P.tipo{k}; c.Rf = P.Rf_ohm_por_fase(k); c.Lf = P.Lf_H(k);
                c.red_R = P.red_R_ohm(k); c.red_L = P.red_L_H(k); c.red_f = P.red_f_Hz(k);
                c.Id_ref_pu = P.Id_ref_pu(k); c.dg1_on = logical(P.dg1_on(k)); c.dg2_on = logical(P.dg2_on(k));
                c.factor = ['T21_' P.factor{k}];
                c.nota = sprintf('T21 %s; theta prevista en Pex %.1f; etiqueta prevista %s', P.id{k}, P.pred_theta_Pex_deg(k), P.pred_label_Pex{k});
                C(end+1) = c; %#ok<AGROW>
            end
        end

    case 'XC'
        % Dos faltas simultaneas en fases distintas y nudos distintos: principal en el bus de Pex,
        % segunda en el nudo de union de lineas (bloques 'Three-Phase Fault' + 'Fault Load1').
        for t2 = {'B-G','C-G'}
            c = base; c.tipo = 'A-G'; c.tipo2 = t2{1}; c.factor = 'cross_country';
            c.nota = ['segunda falta ' t2{1} ' en el nudo de union; la etiqueta de un solo bucle no aplica'];
            C(end+1) = c; %#ok<AGROW>
        end

    case 'S1c'
        % Frecuencia con desviaciones pequenas. En S1a, +-0.5 Hz saca a las DG de su regimen
        % (la consigna del droop sigue en 50 Hz y cada DG desplaza unos 3 kW), asi que se repite
        % con +-0.1 y +-0.2 Hz. (T21.md no existia al lanzar S1: no hay casos tomados de T21.)
        for tp = {'A-B','B-C','C-A'}
            for fr = [49.8 49.9 50.1 50.2]
                c = base; c.tipo = tp{1}; c.Rf = 0.71; c.red_f = fr; c.factor = 'frecuencia';
                C(end+1) = c; %#ok<AGROW>
            end
        end

    otherwise
        error('Conjunto desconocido: %s', conjunto);
end
for k = 1:numel(C)
    C(k).id = sprintf('%s_%s_%03d', conjunto, strrep(C(k).tipo,'-',''), k);
end
end

% =========================================================================
function phi = fase_bucle(archivo, tipo, tref, f)
% Fase (grados, convenio seno) de la tension del bucle en falta en el bus de
% Pex en el instante tref, a partir de los 5 ciclos anteriores.
S = load(archivo); [V, t] = extract_converter_data(S, 'Pex');
tp = upper(strrep(tipo,'-','')); ph = find(ismember('ABC', tp));
if numel(ph) == 1, v = V(:,ph); else, v = V(:,ph(1)) - V(:,ph(2)); end
m = t >= tref - 5/f - 1e-9 & t < tref - 1e-9;
P = 2*mean(v(m).*exp(-1j*2*pi*f*(t(m) - tref)));
phi = mod(rad2deg(angle(P)) + 90, 360);
end

function linea(flog, id, est, ts, msg)
fid = fopen(flog,'a'); fprintf(fid,'%s,%s,%.1f,%s\n', id, est, ts, msg); fclose(fid);
end

function L = leer_logs(dlogs)
L = containers.Map();
d = dir(fullfile(dlogs,'S*_w*.csv'));
for k = 1:numel(d)
    tx = strsplit(fileread(fullfile(d(k).folder,d(k).name)), newline);
    for q = 1:numel(tx)
        p = strsplit(strtrim(tx{q}), ',');
        if numel(p) >= 3, L(p{1}) = {p{2}, str2double(p{3}), strjoin(p(4:end),' ')}; end
    end
end
end
