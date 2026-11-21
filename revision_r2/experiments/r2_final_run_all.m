function r2_final_run_all(varargin)
% R2_FINAL_RUN_ALL  Regenera TODOS los resultados finales de la ronda 2 (T90) con el
% clasificador final (r2_final_opts: puerta de secuencia cero G5, umbral 0.05).
%
%   r2_final_run_all                                % todo, salvo el bloque B de T41 (unos 48 min)
%   r2_final_run_all('bloques', {'nominal','robustez'})
%   r2_final_run_all('persistencia_larga', true)    % incluye el bloque B de T41 (2e7 ventanas)
%
% Bloques, en orden (salidas en _revision_TPWRD\ronda2\resultados\...):
%   nominal          r2_final_nominal                         T90\nominal
%   cca              compare_classifier_accuracy.m (sin modificar) en una copia aislada,
%                    con GA-FuL, para Clarke y PE3D            T90\nominal\cca
%   robustez         r2_final_robustez                         T90\robustez
%   instrumentacion  r2_final_instrumentacion                  T90\instrumentacion
%   emt              r2_final_emt (base, ruido, tex, resumen)  T90\emt
%   lab              r2_final_lab (reutiliza T20\G5_005 si la huella coincide)  T90\lab
%   lab_fig          r2_final_lab_fig                          T90\figuras
%   sintef           r2_t20_run_all('tag','G5_005', 'chain', r2_final_opts())  T20\G5_005
%                    (solo si 'repetir_sintef' es true; si no, se comprueba la huella del log)
%   degenerada       r2_degenerate_security(opts,'final',{'A','C','D'}) [+ {'B'}, 2e7 ventanas]
%                    y r2_degenerate_guard_test(opts,'final')  T41\final
%   puerta           r2_gate_sensitivity(opts,'final')         T61\final
%   macros           python -I r2_t42_t62_macros.py final y copia a T42_cifras.tex y T62_cifras.tex
%   tiempos          r2_final_tiempos (un hilo; no lanzar con otros MATLAB en marcha)  T90\tiempos
%   cifras           python -I r2_final_cifras.py: cifras_finales.csv, diferencias.csv,
%                    cifras_en_textos.csv                       T90
here = fileparts(mfilename('fullpath')); RAIZ = fileparts(here);
R2 = fullfile(RAIZ, '_revision_TPWRD', 'ronda2'); RES = fullfile(R2, 'resultados');
addpath(here, fullfile(RAIZ, 'src', 'fault_classifiers'), fullfile(RAIZ, 'src', 'utils'));
% En la copia publica no estan r2_t42_t62_macros.py ni r2_final_cifras.py (escriben texto del
% manuscrito), asi que 'macros' y 'cifras' no van en la lista por defecto.
p.bloques = {'nominal', 'cca', 'robustez', 'instrumentacion', 'emt', 'lab', 'lab_fig', 'sintef', 'degenerada', 'puerta', 'tiempos'};
p.persistencia_larga = false; p.repetir_sintef = false;
for k = 1:2:numel(varargin)
    if ~isfield(p, varargin{k}), error('r2_final_run_all:opt', 'Opcion desconocida: %s', varargin{k}); end
    p.(varargin{k}) = varargin{k+1};
end
[FO, info] = r2_final_opts();
fprintf('R2_FINAL_RUN_ALL  %s\nOpciones finales: %s\nHuella: %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), info.text, info.md5);
B = cellstr(p.bloques); t0 = tic;
for b = B
    fprintf('\n===== bloque %s (%.0f s)\n', b{1}, toc(t0));
    switch b{1}
        case 'nominal', r2_final_nominal();
        case 'cca'
            C = fullfile(RES, 'T90', 'nominal', 'cca');
            for d = {fullfile(C, 'experiments'), fullfile(C, 'results', 'latex'), fullfile(C, 'results', 'figures'), fullfile(C, 'data')}
                if ~exist(d{1}, 'dir'), mkdir(d{1}); end
            end
            copyfile(fullfile(here, 'compare_classifier_accuracy.m'), fullfile(C, 'experiments'));
            if ~exist(fullfile(C, 'data', 'simulation_curated'), 'dir')
                system(sprintf('cmd /c mklink /J "%s" "%s"', fullfile(C, 'data', 'simulation_curated'), fullfile(RAIZ, 'data', 'simulation_curated')));
            end
            cmd = sprintf(['"%s" -batch "addpath(''' strrep(getenv('GAFUL_PATH'), '\', '\\') '''); addpath(''%s''); addpath(''%s''); ' ...
                'cd(''%s''); compare_classifier_accuracy" > "%s" 2>&1'], fullfile(matlabroot, 'bin', 'matlab.exe'), fullfile(RAIZ, 'src', 'fault_classifiers'), ...
                fullfile(RAIZ, 'src', 'utils'), fullfile(C, 'experiments'), fullfile(RES, 'T90', 'nominal', 'cca_stdout.txt'));
            st = system(cmd); fprintf('compare_classifier_accuracy: estado %d\n', st);
        case 'robustez', correr_script_('r2_final_robustez');   % es un script: se ejecuta en un espacio de trabajo propio
        case 'instrumentacion', r2_final_instrumentacion();
        case 'emt', r2_final_emt();
        case 'lab', r2_final_lab();
        case 'lab_fig', r2_final_lab_fig();
        case 'sintef'
            if p.repetir_sintef
                r2_t20_run_all('tag', 'G5_005', 'chain', FO);
            else
                lg = fileread(fullfile(RES, 'T20', 'G5_005', 'sintef_log.txt'));
                ok = contains(lg, info.md5) && contains(lg, '"v0_mode":"rms_win","kirchhoff_thr":0.05');
                fprintf('T20\\G5_005: huella y opciones iguales a las finales: %d (si es 0, repetir con ''repetir_sintef'', true)\n', ok);
            end
        case 'degenerada'
            r2_degenerate_security(FO, 'final', {'A', 'C', 'D'});
            if p.persistencia_larga, r2_degenerate_security(FO, 'final', {'B'}, 2e7); end
            r2_degenerate_guard_test(FO, 'final');
        case 'puerta', r2_gate_sensitivity(FO, 'final');
        case 'macros'
            system(sprintf('python -I "%s" final', fullfile(here, 'r2_t42_t62_macros.py')));
            copyfile(fullfile(R2, 'textos', 'T42_cifras_final.tex'), fullfile(R2, 'textos', 'T42_cifras.tex'));
            copyfile(fullfile(R2, 'textos', 'T62_cifras_final.tex'), fullfile(R2, 'textos', 'T62_cifras.tex'));
        case 'tiempos', r2_final_tiempos();
        case 'cifras', system(sprintf('python -I "%s"', fullfile(here, 'r2_final_cifras.py')));
        otherwise, error('r2_final_run_all:bloque', 'Bloque desconocido: %s', b{1});
    end
end
fprintf('\nFIN r2_final_run_all (%.0f s). Huella al terminar: %s\n', toc(t0), r2_t20_md5());
end

function correr_script_(nombre)
% Ejecuta un script en el espacio de trabajo de esta subfuncion (no pisa las variables de r2_final_run_all).
eval(nombre);
end
