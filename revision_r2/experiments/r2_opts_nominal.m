function o = r2_opts_nominal(opts, vpk)
% R2_OPTS_NOMINAL  Ajusta v0_nominal de las opciones de r2_decision_chain a cada senal.
%
%   o = r2_opts_nominal(opts, vpk)
%
% Algunas variantes de la puerta de secuencia cero de T52 (v0_mode 'rms_nom' y
% 'ph_nom') necesitan el pico de fase nominal en las unidades de V_abc. Los
% experimentos de T41 y T61 mezclan senales en voltios (simulacion) y en pu
% (eventos sinteticos, laboratorio), asi que el valor no puede ser unico:
% esta funcion lo fija para la senal que se va a clasificar. Con las opciones
% de referencia no cambia nada.
%   vpk  pico de fase sano de la senal, en las unidades de V_abc (antes de
%        dividir por `scale`; la cadena ya divide v0_nominal por scale).
o = opts;
if isfield(o, 'v0_mode') && any(strcmp(o.v0_mode, {'rms_nom', 'ph_nom'})), o.v0_nominal = vpk; end
end
