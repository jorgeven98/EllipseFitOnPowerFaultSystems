function [GridConvs, info] = r2_emt_run_case(caso)
% R2_EMT_RUN_CASE  Simula sin interfaz un caso de falta parametrizado con el
% modelo EMT que genero el dataset del paper (models\GFM.slx, MATLAB R2024b)
% y guarda V_abc e I_abc en bornes de DG1, DG2 y Pex a 4 kHz en el mismo
% formato que los .mat originales (variable GridConvs, legible con
% src\utils\extract_converter_data.m).
%
%   [GridConvs, info] = r2_emt_run_case(caso)
%
% No modifica los .slx: carga una COPIA del modelo desde caso.dir_trabajo
% (si no existe, la copia desde caso.dir_modelo), cambia parametros con
% set_param y nunca guarda el modelo. Todos los parametros se escriben en
% cada llamada, de modo que nada se arrastra del caso anterior aunque el
% modelo siga cargado en la sesion.
%
% Campos de caso (todos opcionales; por defecto, el estado guardado de GFM.slx,
% que es el que genero DataGrid_ABG.mat):
%   id            cadena identificadora                       ('caso')
%   modelo        'GFM' | 'GFM_activeDamping'                  ('GFM')
%   dir_modelo    carpeta con las copias de partida            (ronda2\emt\modelo)
%   dir_trabajo   carpeta de trabajo propia del proceso        (ronda2\emt\work\w0)
%   archivo       .mat de salida                               (dir_trabajo\<id>.mat)
%  Falta principal (bus de Pex, bloques 'Three-Phase Fault1' + 'Fault Load'):
%   tipo          'A-G','B-G','C-G','A-B','B-C','C-A','A-B-G','B-C-G','C-A-G',
%                 'A-B-C','A-B-C-G' o 'NONE' (sin falta)        ('A-B-G')
%   Rf            resistencia POR FASE en falta, ohm: suma de la Ron del bloque
%                 de falta y de la R de la rama serie 'Fault Load'   (0.003)
%   Ron           parte de Rf que se deja en el interruptor     (min(0.002, Rf/2))
%                 OJO: con Rf = 0.003 el reparto por defecto es 1.5 + 1.5 mohm. El modelo
%                 original reparte 2 mohm (interruptor) + 1 mohm (rama); para reproducirlo
%                 bit a bit hay que pasar Ron = 0.002. Los dos repartos son electricamente
%                 equivalentes (misma resistencia serie) y solo cambia el rizado de conmutacion.
%   Rtierra       resistencia a tierra de la falta (GroundResistance), ohm (0.01)
%   Lf            inductancia serie por fase de la rama de falta, H     (0)
%   t_falta       instante de incepcion, s (resolucion = paso 1e-5 s)   (3.0)
%   dur_falta     duracion de la falta, s                                (0.16)
%   t_stop        fin de la simulacion, s             (t_falta + dur_falta + 0.14)
%  Segunda falta opcional (nudo de union de lineas, bloques 'Three-Phase Fault'
%  + 'Fault Load1', solo en GFM): tipo2 ('NONE'), Rf2, Ron2, Rtierra2, Lf2,
%   t_falta2, dur_falta2 (por defecto, los de la falta principal).
%  Punto de operacion:
%   Id_ref_pu     consigna de corriente activa de Pex en p.u. de Ibase = 40 A
%                 (escalon 'ref. Export  current1', negativo = importa)  (-0.15)
%   dg1_dP        escalon de potencia de DG1 en t = 1.5 s, W             (800)
%  Equivalente de red (bloque 'Three-Phase Source'):
%   red_R, red_L  resistencia (ohm) e inductancia (H) serie      (0.8929, 16.58e-3)
%   red_f         frecuencia de la fuente, Hz                    (50)
%   red_Vll       tension de linea eficaz de la fuente, V        (Vbase*sqrt(3)/sqrt(2))
%  Aportacion remota (se comentan las lineas; el control no se toca):
%   dg1_on        true = DG1 conectada (lineas 'Line' y 'Line3') (true)
%   dg2_on        true = DG2 conectada (linea 'Line2')           (true)
%
% Salida: GridConvs (struct de timeseries: DG1.Vac, DG1.Io, DG2.Vac, DG2.Io,
% Pex.Vabc, Pex.Io, Pex.Vdc, Pex.Idc) e info (parametros aplicados, tiempo de
% simulacion, archivo). El .mat guardado contiene GridConvs y la variable
% 'caso' con todos los parametros efectivos.
%
% Nota sobre el equivalente de red: en la mascara, la casilla "Specify
% short-circuit level parameters" esta desmarcada (SpecifyImpedance = off) y
% en ese estado el bloque usa directamente Resistance e Inductance (medido
% con un circuito minimo, ver ronda2\informes\T22.md). Los campos
% ShortCircuitLevel = 100e6 y XRratio = 7 NO estan activos.
%
% Ronda 2 TPWRD, tareas T22 y T23.

if nargin < 1, caso = struct(); end
raiz2 = [fileparts(fileparts(mfilename('fullpath'))) '\_revision_TPWRD\ronda2'];
d = struct( ...
    'id','caso', 'modelo','GFM', ...
    'dir_modelo', fullfile(raiz2,'emt','modelo'), ...
    'dir_trabajo', fullfile(raiz2,'emt','work','w0'), ...
    'archivo','', ...
    'tipo','A-B-G', 'Rf',0.003, 'Ron',[], 'Rtierra',0.01, 'Lf',0, ...
    't_falta',3.0, 'dur_falta',0.16, 't_stop',[], ...
    'tipo2','NONE', 'Rf2',0.003, 'Ron2',[], 'Rtierra2',0.01, 'Lf2',0, ...
    't_falta2',[], 'dur_falta2',[], ...
    'Id_ref_pu',-0.15, 'dg1_dP',800, ...
    'red_R',0.8929, 'red_L',16.58e-3, 'red_f',50, 'red_Vll',[], ...
    'dg1_on',true, 'dg2_on',true);
fn = fieldnames(d);
for k = 1:numel(fn)
    if ~isfield(caso, fn{k}) || (isempty(caso.(fn{k})) && ~isempty(d.(fn{k})))
        caso.(fn{k}) = d.(fn{k});
    end
end
extra = setdiff(fieldnames(caso), [fn; {'meta'}]);
if ~isempty(extra)
    error('r2_emt_run_case:CampoDesconocido', 'Campo no reconocido: %s', strjoin(extra', ', '));
end
if isempty(caso.Ron),  caso.Ron  = min(0.002, caso.Rf/2);  end
if isempty(caso.Ron2), caso.Ron2 = min(0.002, caso.Rf2/2); end
if isempty(caso.t_stop), caso.t_stop = caso.t_falta + caso.dur_falta + 0.14; end
if isempty(caso.t_falta2), caso.t_falta2 = caso.t_falta; end
if isempty(caso.dur_falta2), caso.dur_falta2 = caso.dur_falta; end
if isempty(caso.archivo), caso.archivo = fullfile(caso.dir_trabajo, [caso.id '.mat']); end
paso = 1e-5;                                    % paso fijo del modelo (ode3, powergui discreto)
caso.t_falta  = round(caso.t_falta/paso)*paso;  % el interruptor solo conmuta en multiplos del paso
caso.t_falta2 = round(caso.t_falta2/paso)*paso;

mdl = caso.modelo;

% --- carpeta de trabajo con copia propia del modelo (slprj no se comparte)
if ~exist(caso.dir_trabajo,'dir'), mkdir(caso.dir_trabajo); end
if ~exist(fullfile(caso.dir_trabajo,[mdl '.slx']),'file')
    copyfile(fullfile(caso.dir_modelo,'*'), caso.dir_trabajo);
end
dir0 = pwd; vuelve = onCleanup(@() cd(dir0));
cd(caso.dir_trabajo);

% --- carga (PreLoadFcn = GridDesign hace clc y clear all en el workspace base;
%     las variables locales de esta funcion no se ven afectadas)
if ~bdIsLoaded(mdl)
    evalc('load_system(mdl)');
end
if ~strcmpi(fileparts(get_param(mdl,'FileName')), caso.dir_trabajo)
    error('r2_emt_run_case:ModeloAjeno', ...
        'El modelo cargado (%s) no es la copia de %s', get_param(mdl,'FileName'), caso.dir_trabajo);
end
if evalin('base','~exist(''Vbase'',''var'')')   % por si otra carga limpio el workspace base
    evalc('evalin(''base'',''GridDesign'')');
end

% --- bloques
b.src   = [mdl '/Three-Phase Source'];
b.f1    = [mdl '/Three-Phase Fault1'];
b.fl1   = [mdl '/Fault Load'];
b.f2    = [mdl '/Three-Phase Fault'];
b.fl2   = [mdl '/Fault Load1'];
b.tofile = char(find_system(mdl,'SearchDepth',1,'BlockType','ToFile'));
b.idref = find_system(mdl,'LookUnderMasks','all','Regexp','on','BlockType','Step', ...
    'Name','^ref\. Export\s+current1$');
b.dg1step = find_system(mdl,'LookUnderMasks','all','Regexp','on','BlockType','Step', ...
    'Parent','DG1.*Primary control$','Name','^Step$');
if numel(b.idref) ~= 1, error('r2_emt_run_case:Bloque','No se localiza el escalon de Id de Pex'); end
if numel(b.dg1step) ~= 1, error('r2_emt_run_case:Bloque','No se localiza el escalon de potencia de DG1'); end
hay_f2 = getSimulinkBlockHandle(b.f2) > 0 && getSimulinkBlockHandle(b.fl2) > 0;

% --- falta principal
[fa, fb, fc, fg, hay] = fases(caso.tipo);
set_param(b.f1, 'Commented','off'); set_param(b.fl1, 'Commented','off');
if hay
    sw = [caso.t_falta, caso.t_falta + caso.dur_falta];
else
    fa = 'on'; fb = 'on'; fc = 'off'; fg = 'off';      % bloque presente pero sin conmutar
    sw = caso.t_stop + [10 11];
end
set_param(b.f1, 'InitialStates','0', 'FaultA',fa, 'FaultB',fb, 'FaultC',fc, 'GroundFault',fg, ...
    'External','off', 'SwitchTimes', sprintf('[%.5f %.5f]', sw(1), sw(2)), ...
    'FaultResistance', num2str(caso.Ron,12), 'GroundResistance', num2str(caso.Rtierra,12));
rama(b.fl1, caso.Rf - caso.Ron, caso.Lf);

% --- segunda falta (nudo de union)
[fa2, fb2, fc2, fg2, hay2] = fases(caso.tipo2);
if hay2 && ~hay_f2
    error('r2_emt_run_case:SinSegundaFalta','El modelo %s no tiene los bloques de la segunda falta', mdl);
end
if hay_f2
    if hay2
        set_param(b.f2, 'Commented','off'); set_param(b.fl2, 'Commented','off');
        sw2 = [caso.t_falta2, caso.t_falta2 + caso.dur_falta2];
        set_param(b.f2, 'InitialStates','0', 'FaultA',fa2, 'FaultB',fb2, 'FaultC',fc2, 'GroundFault',fg2, ...
            'External','off', 'SwitchTimes', sprintf('[%.5f %.5f]', sw2(1), sw2(2)), ...
            'FaultResistance', num2str(caso.Ron2,12), 'GroundResistance', num2str(caso.Rtierra2,12));
        rama(b.fl2, caso.Rf2 - caso.Ron2, caso.Lf2);
    else
        set_param(b.f2, 'Commented','on'); set_param(b.fl2, 'Commented','on');
    end
end

% --- punto de operacion
set_param(b.idref{1}, 'Time','1', 'Before','Ibase*0.0', 'After', sprintf('Ibase*(%.10g)', caso.Id_ref_pu));
set_param(b.dg1step{1}, 'Time','1.5', 'Before','0', 'After', num2str(caso.dg1_dP,12));

% --- equivalente de red (SpecifyImpedance = off usa Resistance e Inductance)
if isempty(caso.red_Vll), vll = 'Vbase*sqrt(3)/sqrt(2)'; else, vll = num2str(caso.red_Vll,12); end
set_param(b.src, 'SpecifyImpedance','off', 'NonIdealSource','on', ...
    'Resistance', num2str(caso.red_R,12), 'Inductance', num2str(caso.red_L,12), ...
    'Frequency', num2str(caso.red_f,12), 'Voltage', vll);

% --- aportacion remota
set_param([mdl '/Line2'], 'Commented', onoff(~caso.dg2_on));
set_param([mdl '/Line'],  'Commented', onoff(~caso.dg1_on));
set_param([mdl '/Line3'], 'Commented', onoff(~caso.dg1_on));

% --- simulacion
tmp = fullfile(caso.dir_trabajo, 'r2tmp_DataGrid.mat');
if exist(tmp,'file'), delete(tmp); end
set_param(b.tofile, 'Filename', tmp);
set_param(mdl, 'StopTime', num2str(caso.t_stop,12));
t0 = tic;
evalc('sim(mdl)');
info.t_sim = toc(t0);

% --- salida en el formato del dataset + parametros
S = load(tmp);
GridConvs = S.GridConvs;
caso.paso_s = paso;
caso.fs_registro = 1/250e-6;
caso.matlab = version;
caso.t_sim_s = info.t_sim;
caso.fecha = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
dsal = fileparts(caso.archivo);
if ~isempty(dsal) && ~exist(dsal,'dir'), mkdir(dsal); end
movefile(tmp, caso.archivo, 'f');
save(caso.archivo, 'caso', '-append');
info.caso = caso;
info.archivo = caso.archivo;
end

% =========================================================================
function [fa, fb, fc, fg, hay] = fases(tipo)
t = upper(strrep(strrep(char(tipo),'-',''),' ',''));
hay = ~any(strcmp(t, {'NONE','','SANO'}));
if ~hay, fa = 'off'; fb = 'off'; fc = 'off'; fg = 'off'; return, end
if ~isempty(regexp(t,'[^ABCG]','once')), error('r2_emt_run_case:Tipo','Tipo de falta no valido: %s', tipo); end
fa = onoff(any(t=='A')); fb = onoff(any(t=='B')); fc = onoff(any(t=='C')); fg = onoff(any(t=='G'));
end

function s = onoff(tf)
if tf, s = 'on'; else, s = 'off'; end
end

function rama(blk, R, L)
% Rama serie trifasica de la falta: R por fase (>= 0) y L por fase.
if R < 0, error('r2_emt_run_case:Rf','Rf debe ser >= Ron'); end
R = max(R, 1e-6);       % la rama no admite R = 0 con BranchType R
if L > 0
    set_param(blk, 'BranchType','RL', 'Resistance', num2str(R,12), 'Inductance', num2str(L,12));
else
    set_param(blk, 'BranchType','R', 'Resistance', num2str(R,12));
end
end
