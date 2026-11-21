function r2_t73_xc()
% R2_T73_XC  Dobles faltas simultaneas de la campana EMT (casos XC de T23) con el clasificador FINAL.
%   A-G en el bus de Pex y B-G (XC_AG_001) o C-G (XC_AG_002) a la vez en el nudo de union de lineas.
%   La decision sale de r2_decision_chain con r2_final_opts (puerta G5, umbral 0.05), como en r2_final_emt.m;
%   aqui no se reimplementa nada. La "verdad" que se pasa es la union de fases en falta (A-B-G o C-A-G):
%   es solo una etiqueta de referencia para las metricas, no una clase que el metodo deba dar.
%   Salida: _revision_TPWRD/ronda2/resultados/T73/xc_final.csv y xc_final_log.txt.
%   Uso: matlab -batch "addpath('...\experiments'); r2_t73_xc"
here = fileparts(mfilename('fullpath')); addpath(here);
addpath(fullfile(here, '..', 'src', 'fault_classifiers')); addpath(fullfile(here, '..', 'src', 'utils'));
R2 = fullfile(here, '..', '_revision_TPWRD', 'ronda2');
OUT = fullfile(R2, 'resultados', 'T73'); if ~exist(OUT, 'dir'), mkdir(OUT); end
diary(fullfile(OUT, 'xc_final_log.txt')); diary on
[o, info] = r2_final_opts();
fprintf('R2_T73_XC %s\nOpciones finales: %s\n', datestr(now, 31), info.text);
casos = {'XC_AG_001', 'A-B-G'; 'XC_AG_002', 'C-A-G'};
TERM = {'Pex', 'DG1'};
rows = {};
for k = 1:size(casos, 1)
    S = load(fullfile(R2, 'emt', 'casos', [casos{k, 1} '.mat'])); c = S.caso;
    t_on = c.t_falta; t_off = t_on + c.dur_falta; f = c.red_f;
    for q = 1:2
        [V, t] = extract_converter_data(S, TERM{q}); V = double(V); t = double(t(:));
        out = r2_decision_chain(V, t, o, 'truth', casos{k, 2}, 't_on', t_on, 't_off', t_off, 'f', f, ...
            't_range', [t_on - 1.0, min(t(end), t_off + 0.12)]);
        e = out.event;
        r = struct('id', casos{k, 1}, 'terminal', TERM{q}, 'union_fases', casos{k, 2});
        fn = fieldnames(e);
        for j = 1:numel(fn)
            v = e.(fn{j});
            if (isnumeric(v) || islogical(v)) && isscalar(v), r.(fn{j}) = double(v);
            elseif ischar(v) || isstring(v), r.(fn{j}) = char(v);
            end
        end
        rows{end+1} = r; %#ok<AGROW>
        fprintf('%s %s: ', casos{k, 1}, TERM{q}); disp(r);
    end
end
T = struct2table([rows{:}]);
writetable(T, fullfile(OUT, 'xc_final.csv'));
fprintf('FIN\n'); diary off
end
