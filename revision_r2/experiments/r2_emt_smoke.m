function r2_emt_smoke(solo_comparar)
% R2_EMT_SMOKE  Prueba de humo del banco EMT (T22): re-simula con
% r2_emt_run_case los casos DataGrid_ABG (GFM.slx tal cual) y DataGridRL_AB
% (parametros deducidos en los informes de procedencia) y mide la diferencia
% con los registros originales: instantanea, de fasores de 50 Hz por ciclo y
% de etiquetas ventana a ventana. Tambien compara con la ejecucion del modelo
% sin tocar (humo0) para comprobar determinismo y que set_param no altera nada.
%
% Salidas en ronda2\resultados\T22\: humo_*.mat, humo_diferencias.csv, humo_log.txt
if nargin < 1, solo_comparar = false; end
raiz  = fileparts(fileparts(mfilename('fullpath')));
raiz2 = fullfile(raiz,'_revision_TPWRD','ronda2');
res   = fullfile(raiz2,'resultados','T22');
addpath(fullfile(raiz,'experiments'), fullfile(raiz,'src','fault_classifiers'), fullfile(raiz,'src','utils'));
wd = fullfile(raiz2,'emt','work','humo');

casos = struct('id',{},'ref',{},'c',{});
% 1) ABG: estado guardado de GFM.slx
casos(end+1) = struct('id','humo_ABG_GFM', 'ref','DataGrid_ABG.mat', ...
    'c', struct('modelo','GFM','tipo','A-B-G','Rf',0.003,'Rtierra',0.01,'Lf',0,'Id_ref_pu',-0.15,'t_stop',3.4));
% 2) RL_AB con GFM_activeDamping (modelo deducido), Fault Load RL 1 mohm + 0.5 mH, Id = -0.05*Ibase
casos(end+1) = struct('id','humo_RLAB_AD', 'ref','DataGridRL_AB.mat', ...
    'c', struct('modelo','GFM_activeDamping','tipo','A-B','Rf',0.003,'Lf',0.5e-3,'Id_ref_pu',-0.05,'t_stop',3.4));
% 3) RL_AB con GFM (para decidir el modelo de campana)
casos(end+1) = struct('id','humo_RLAB_GFM', 'ref','DataGridRL_AB.mat', ...
    'c', struct('modelo','GFM','tipo','A-B','Rf',0.003,'Lf',0.5e-3,'Id_ref_pu',-0.05,'t_stop',3.4));
% 4) ABG repetido (determinismo dentro de la misma sesion, tras otros casos)
casos(end+1) = struct('id','humo_ABG_GFM_rep', 'ref','DataGrid_ABG.mat', ...
    'c', struct('modelo','GFM','tipo','A-B-G','Rf',0.003,'Rtierra',0.01,'Lf',0,'Id_ref_pu',-0.15,'t_stop',3.4));

% 5 y 6) Reparto fiel al modelo original: 2 mohm en el interruptor (Ron) y 1 mohm en la rama serie.
%    Por defecto r2_emt_run_case reparte Rf = 3 mohm en 1.5 + 1.5 mohm (electricamente equivalente).
casos(end+1) = struct('id','humo_ABG_GFM_Ron2', 'ref','DataGrid_ABG.mat', ...
    'c', struct('modelo','GFM','tipo','A-B-G','Rf',0.003,'Ron',0.002,'Rtierra',0.01,'Lf',0,'Id_ref_pu',-0.15,'t_stop',3.4));
casos(end+1) = struct('id','humo_RLAB_AD_Ron2', 'ref','DataGridRL_AB.mat', ...
    'c', struct('modelo','GFM_activeDamping','tipo','A-B','Rf',0.003,'Ron',0.002,'Lf',0.5e-3,'Id_ref_pu',-0.05,'t_stop',3.4));

flog = fopen(fullfile(res,'humo_log.txt'),'a');   % el CSV de diferencias se reescribe entero en cada llamada
fprintf(flog,'\n==== r2_emt_smoke %s (MATLAB %s)\n', char(datetime('now')), version);
if ~solo_comparar
    for k = 1:numel(casos)
        c = casos(k).c; c.id = casos(k).id; c.dir_trabajo = wd; c.archivo = fullfile(res,[casos(k).id '.mat']);
        if exist(c.archivo,'file'), continue, end      % no se repite lo ya simulado
        try
            [~, info] = r2_emt_run_case(c);
            fprintf(flog,'%s: simulado en %.1f s\n', casos(k).id, info.t_sim);
        catch ME
            fprintf(flog,'%s: ERROR %s\n', casos(k).id, ME.message);
        end
    end
end

% ---- comparaciones
pares = {};
for k = 1:numel(casos)
    pares(end+1,:) = {casos(k).id, fullfile(res,[casos(k).id '.mat']), fullfile(raiz,'data','simulation_curated',casos(k).ref), casos(k).ref}; %#ok<AGROW>
end
pares(end+1,:) = {'humo_ABG_GFM', fullfile(res,'humo_ABG_GFM.mat'), fullfile(res,'humo0_GFM_sin_tocar.mat'), 'humo0 (GFM.slx sin set_param)'};
pares(end+1,:) = {'humo_ABG_GFM_rep', fullfile(res,'humo_ABG_GFM_rep.mat'), fullfile(res,'humo_ABG_GFM.mat'), 'humo_ABG_GFM (repeticion)'};
pares(end+1,:) = {'S2a_ABG_021 (otra sesion)', fullfile(raiz2,'emt','casos','S2a_ABG_021.mat'), fullfile(res,'humo_ABG_GFM.mat'), 'humo_ABG_GFM (mismo caso; otra sesion de MATLAB)'};
pares(end+1,:) = {'humo_ABG_GFM_Ron2', fullfile(res,'humo_ABG_GFM_Ron2.mat'), fullfile(res,'humo0_GFM_sin_tocar.mat'), 'humo0 (GFM.slx sin set_param)'};
pares(end+1,:) = {'humo0', fullfile(res,'humo0_GFM_sin_tocar.mat'), fullfile(raiz,'data','simulation_curated','DataGrid_ABG.mat'), 'DataGrid_ABG.mat'};

fs = 4000; f = 50; W = 20; PASO = 5; Vn = 231;          % Vn: pico de fase nominal (Vbase), base del p.u.
fid = fopen(fullfile(res,'humo_diferencias.csv'),'w');
fprintf(fid,['nuevo,referencia,convertidor,tramo,rms_dif_V,rms_dif_pu_pct,max_dif_V,' ...
    'rms_dif_filtrada_V,rms_dif_filtrada_pu_pct,dif_fasor_V1_pct,dif_fasor_max_pct,' ...
    'n_ventanas,etiquetas_iguales_pct,etq_mayoritaria_nuevo,etq_mayoritaria_ref\n']);
for k = 1:size(pares,1)
    if ~exist(pares{k,2},'file') || ~exist(pares{k,3},'file'), continue, end
    A = load(pares{k,2}); B = load(pares{k,3});
    for cn = {'DG1','DG2','Pex'}
        [va, ta] = extract_converter_data(A, cn{1}); [vb, ~] = extract_converter_data(B, cn{1});
        n = min(size(va,1), size(vb,1)); va = va(1:n,:); vb = vb(1:n,:); t = ta(1:n);
        tramos = {'prefalta', 2.5, 3.0; 'falta', 3.0, 3.16};
        for q = 1:2
            m = t >= tramos{q,2} & t < tramos{q,3} - 1e-9;
            dd = va(m,:) - vb(m,:);
            % diferencia tras quitar el rizado de conmutacion: media movil de 1 ms (4 muestras) en ambos
            fa_ = movmean(va,4,1); fb_ = movmean(vb,4,1); df = fa_(m,:) - fb_(m,:);
            % fasores de 50 Hz por ciclo (80 muestras)
            i0 = find(m,1); nc = floor(nnz(m)/80); dfas = zeros(nc,1); dv1 = zeros(nc,1);
            e = exp(-1j*2*pi*f*t);
            for c = 1:nc
                ii = i0 + (c-1)*80 + (0:79);
                Pa = 2*mean(va(ii,:).*e(ii),1); Pb = 2*mean(vb(ii,:).*e(ii),1);
                dfas(c) = max(abs(Pa - Pb))/Vn;
                a = exp(2j*pi/3); s1 = [1 a a^2]/3;
                dv1(c) = abs(s1*Pa.' - s1*Pb.')/Vn;
            end
            % etiquetas ventana a ventana (clasificador puro, voltios crudos como en la cadena nominal)
            idx = find(m); la = {}; lb = {};
            for s = idx(1):PASO:(idx(end)-W+1)
                la{end+1} = etq(va(s:s+W-1,:), fs, f); lb{end+1} = etq(vb(s:s+W-1,:), fs, f); %#ok<AGROW>
            end
            fprintf(fid,'%s,%s,%s,%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%d,%.1f,%s,%s\n', pares{k,1}, pares{k,4}, cn{1}, tramos{q,1}, ...
                sqrt(mean(dd.^2,'all')), 100*sqrt(mean(dd.^2,'all'))/(Vn/sqrt(2)), max(abs(dd),[],'all'), ...
                sqrt(mean(df.^2,'all')), 100*sqrt(mean(df.^2,'all'))/(Vn/sqrt(2)), 100*mean(dv1), 100*max(dfas), ...
                numel(la), 100*mean(strcmp(la,lb)), moda(la), moda(lb));
        end
    end
end
fclose(fid); fclose(flog);
end

function l = etq(V, fs, f)
ws = warning('off','all'); c = onCleanup(@() warning(ws));
try, l = char(ga_fault_classifier_pure(V, fs, f)); catch, l = 'Error'; end
end

function m = moda(l)
[u,~,j] = unique(l); [~,i] = max(accumarray(j(:),1)); m = u{i};
end
