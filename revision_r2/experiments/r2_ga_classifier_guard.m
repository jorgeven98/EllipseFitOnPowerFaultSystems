function [fault_type, params] = r2_ga_classifier_guard(V_abc, fs, f, G)
% R2_GA_CLASSIFIER_GUARD  Envoltorio de prueba para la decision D5 (T41, ronda 2).
%
% Llama al clasificador base SIN modificarlo y, despues, aplica condiciones
% de amplitud medidas sobre la propia ventana por ajuste sinusoidal de
% minimos cuadrados (v = c cos wt + s sin wt por fase; semieje mayor = mayor
% valor singular de [c s]). Si una condicion bloquea la ventana, la etiqueta
% pasa a 'Unclassified', que la cadena no declara y que rompe la racha
% ('A-B-C' no vale: la puerta de redondez de la cadena la recalcularia).
%
%   [etq, p] = r2_ga_classifier_guard(Vw, fs, f, G)
%
% Campos de G (NaN o ausente = condicion desactivada):
%   base        clasificador base, [etq, p] = base(Vw, fs, f)
%               (por defecto @ga_fault_classifier_pure)
%   v_nominal   semieje mayor nominal en las unidades de Vw: sqrt(3/2) por el
%               pico de fase sano (1.2247 en pu; unos 274 V en Pex)
%   a_min_line  SEMIEJE MAYOR EN LA RAMA DEGENERADA (lo que describia el texto
%               del paper): si la ventana se trato como linea, la etiqueta no
%               es de tierra y el semieje mayor medido es menor que
%               a_min_line * v_nominal, se bloquea
%   v_min       TENSION MINIMA (linea sin tension): cualquier etiqueta se
%               bloquea si el semieje mayor medido es menor que v_min * v_nominal
%   uv_start    ARRANQUE POR SUBTENSION: una etiqueta de falta solo se acepta
%               si la amplitud de alguna fase es menor que uv_start por la
%               nominal de fase (v_nominal/sqrt(3/2))
% Salida: los campos del clasificador base y, ademas, a_fit, b_fit,
% amp_min_rel, a_fit_rel y guard (texto: '', 'a_min_line', 'v_min', 'uv_start').
% Uso con la cadena:
%   G = struct('v_nominal', sqrt(3/2)*Npk, 'v_min', 0.1);
%   o = r2_decision_chain(V, t, 'classifier', @(Vw, fs, f) r2_ga_classifier_guard(Vw, fs, f, G));
if ~isfield(G, 'base') || isempty(G.base), G.base = @ga_fault_classifier_pure; end
[fault_type, params] = G.base(V_abc, fs, f);
fault_type = char(fault_type);
N = size(V_abc, 1); tw = (0:N-1)'/fs; X = [cos(2*pi*f*tw), sin(2*pi*f*tw)];
C = X \ V_abc; C(~isfinite(C)) = 0;                      % 2 x 3
s = svd(C); amp = sqrt(sum(C.^2, 1));
params.a_fit = s(1); params.b_fit = s(2);
params.a_fit_rel = s(1)/G.v_nominal; params.amp_min_rel = min(amp)/(G.v_nominal/sqrt(3/2));
params.guard = '';
GROUND = {'A-G','B-G','C-G','A-B-G','B-C-G','C-A-G'};
is_fault = ~ismember(fault_type, {'A-B-C', 'Error', 'Unclassified'});
if ~is_fault, return, end
on = @(name) isfield(G, name) && ~isempty(G.(name)) && isfinite(G.(name));
if on('v_min') && params.a_fit_rel < G.v_min
    params.guard = 'v_min';
elseif on('a_min_line') && isfield(params, 'is_line') && params.is_line && ~ismember(fault_type, GROUND) && params.a_fit_rel < G.a_min_line
    params.guard = 'a_min_line';
elseif on('uv_start') && params.amp_min_rel >= G.uv_start
    params.guard = 'uv_start';
end
if ~isempty(params.guard), fault_type = 'Unclassified'; end
end
