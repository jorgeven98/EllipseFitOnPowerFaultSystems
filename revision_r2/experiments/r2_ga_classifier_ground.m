function [fault_type, params] = r2_ga_classifier_ground(V_abc, fs, f, P)
% R2_GA_CLASSIFIER_GROUND  Clasificador GA con la regla corregida de fase en la
% rama de tierra (T23b). Envoltura de ga_fault_classifier_pure, que no se
% modifica: la deteccion de tierra (puerta de Kirchhoff), la rama entre fases
% y todos los parametros devueltos son los del clasificador publicado. Solo
% cambia la etiqueta cuando la ventana va por la rama de tierra y la
% componente inversa es medible y proporcionada a la homopolar.
%
%   [etq, p] = r2_ga_classifier_ground(Vw, fs, f)        % valores por defecto
%   [etq, p] = r2_ga_classifier_ground(Vw, fs, f, P)     % P struct de opciones
%   P = r2_ga_classifier_ground('defaults')
%
% Opciones:
%   u2_min    0.01   |V2/V1| minimo estimado en la ventana para aplicar la
%                    regla nueva (suelo de medida, no umbral de clase).
%   v2v0_min  0.15   |V2/V0| minimo. En el punto de falta, |V2/V0| = |Z2/Z0|
%                    en una L-G con cualquier Rf (0.22 en la red del paper).
%                    Por debajo, la tension inversa no informa de la red
%                    (bornes de un convertidor que la regula, o red con
%                    |Z0| > 6.7 |Z2|: neutro aislado, resonante o por
%                    resistencia alta) y se devuelve la etiqueta publicada.
%                    El valor 0.15 se ha elegido con los registros de este
%                    proyecto (resultados\T23b\orig_v2v0.csv): en DG1 el
%                    cociente llega a 0.20 y en Pex su mediana es 0.22.
% Condicion de red de la regla nueva (no comprobable en la ventana):
% |arg Z0 - arg Z2| < 60 grados. Con Z0 resistiva la regla puede activarse y
% fallar donde la publicada acierta (informes\T23b.md, apartado 3.3).
%   metodo    'conj' estimador de u0 y u2 (r2_ground_rf_estimate)
% Campos anadidos a p: label_legacy, regla ('nueva','legacy','no_tierra'),
% u0, u2, alfa_deg, beta_deg, rK, paso.
% Compatible con r2_decision_chain(..., 'classifier', @r2_ga_classifier_ground).
D = struct('u2_min', 0.01, 'v2v0_min', 0.15, 'metodo', 'conj');
if (ischar(V_abc) || isstring(V_abc)) && strcmpi(V_abc, 'defaults')
    fault_type = D; params = []; return
end
if nargin >= 4 && ~isempty(P)
    fn = fieldnames(P); for i = 1:numel(fn), D.(fn{i}) = P.(fn{i}); end
end
if nargin < 3 || isempty(f), f = 50; end
[fault_type, params] = ga_fault_classifier_pure(V_abc, fs, f);
fault_type = char(fault_type);
params.label_legacy = fault_type; params.regla = 'no_tierra';
params.u0 = NaN; params.u2 = NaN; params.alfa_deg = NaN; params.beta_deg = NaN; params.rK = NaN; params.paso = '';
if ~ismember(fault_type, {'A-G','B-G','C-G','A-B-G','B-C-G','C-A-G'}), return, end
E = r2_ground_rf_estimate(V_abc, fs, f, D.metodo);
params.u0 = E.u0; params.u2 = E.u2; params.rK = E.rK;
params.regla = 'legacy';
if isfinite(E.u2) && isfinite(E.u0) && abs(E.u2) >= D.u2_min && abs(E.u2) >= D.v2v0_min*abs(E.u0)
    [lab, info] = r2_ground_rf_rule(E.u0, E.u2);
    if ~isempty(lab)
        fault_type = lab; params.regla = 'nueva';
        params.alfa_deg = info.alfa_deg; params.beta_deg = info.beta_deg; params.paso = info.paso;
    end
end
end
