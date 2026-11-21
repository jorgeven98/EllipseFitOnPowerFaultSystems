function [opts, info] = r2_final_opts()
% R2_FINAL_OPTS  Opciones FINALES de la cadena de decision de la ronda 2 (decision D2).
%
%   opts = r2_final_opts();
%   out  = r2_decision_chain(V, t, opts, 'truth', 'A-B', 't_on', 3.0, 't_off', 3.16);
%   [opts, info] = r2_final_opts();   % info: opciones resueltas, huella MD5 y texto
%
% Clasificador final = el publicado con la puerta de secuencia cero G5 de
% umbral 0.05 (T52): la ventana va a la rama sin tierra si
% ||X k||/||X||_F < 0.05 (X ventana completa N x 3, k = (1,1,1)/sqrt(3)).
% Todo lo demas es lo publicado y sale de los valores por defecto de
% r2_decision_chain: M = 10, Delta = 0.01, K = 80 ventanas, puerta adaptativa
% de redondez con nivel ambiente causal, clase de tierra de dos muestras
% (gnd_mode = '2s'), marcas temporales al final de la ventana y scale = 1
% (voltios en simulacion y EMT, pu en laboratorio, como en la ronda 1).
%
% opts solo lleva lo que difiere de la cadena de referencia, para que los
% scripts que barren Delta, K o M (r2_gate_sensitivity) puedan hacerlo. La
% funcion comprueba que los valores por defecto de la cadena siguen siendo
% los publicados y da error si alguien los ha cambiado.
%
% TODO script r2_final_* toma de aqui las opciones. No se escriben en ningun
% otro sitio.
here = fileparts(mfilename('fullpath')); addpath(here);
addpath(fullfile(here, '..', 'src', 'fault_classifiers')); addpath(fullfile(here, '..', 'src', 'utils'));
opts = struct('v0_mode', 'rms_win', 'kirchhoff_thr', 0.05);
d = r2_decision_chain('defaults');
chk = {'M', 10; 'Delta', 0.01; 'K', 80; 'thr_cap', 0.99; 'scale', 1};
for k = 1:size(chk, 1)
    if ~isequal(d.(chk{k, 1}), chk{k, 2})
        error('r2_final_opts:defaults', 'El valor por defecto de %s en r2_decision_chain ya no es el publicado.', chk{k, 1});
    end
end
txt = {'gate', 'adaptive'; 'amb_mode', 'causal'; 'gnd_mode', '2s'; 'healthy_label', 'A-B-C'; 'amb_init', 'assume_healthy'; 'amb_feed', 'healthy_label'};
for k = 1:size(txt, 1)
    if ~strcmp(char(d.(txt{k, 1})), txt{k, 2})
        error('r2_final_opts:defaults', 'El valor por defecto de %s en r2_decision_chain ya no es el publicado.', txt{k, 1});
    end
end
if d.amb_holdoff ~= 0, error('r2_final_opts:defaults', 'amb_holdoff por defecto ya no es 0.'); end
if nargout > 1
    r = d; fn = fieldnames(opts);
    for k = 1:numel(fn), r.(fn{k}) = opts.(fn{k}); end
    info.resolved = r;
    info.md5 = r2_t20_md5();
    info.text = sprintf('v0_mode=%s kirchhoff_thr=%.3g gnd_mode=%s gate=%s amb_mode=%s M=%d Delta=%.3g K=%d scale=%g', ...
        r.v0_mode, r.kirchhoff_thr, r.gnd_mode, r.gate, r.amb_mode, r.M, r.Delta, r.K, r.scale);
end
end
