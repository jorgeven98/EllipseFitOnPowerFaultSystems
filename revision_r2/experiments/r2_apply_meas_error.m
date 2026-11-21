function [Ve, info] = r2_apply_meas_error(V, t, varargin)
% R2_APPLY_MEAS_ERROR  Errores de medida sobre la senal COMPLETA (antes de segmentar).
%
%   Ve = r2_apply_meas_error(V, t, 'phase_deg', [0.5 0 0])      % retardo en A
%   Ve = r2_apply_meas_error(V, t, 'gain', [0.01 0 0])          % +1 % en A
%   Ve = r2_apply_meas_error(V, t, 'snr_db', 30, 'seed', 1)     % ruido blanco
%
% Opciones (se aplican en este orden: fase, ganancia, ruido):
%   phase_deg [1x3]  error de fase por canal en grados a la frecuencia f.
%                    Signo positivo = RETARDO temporal de phase_deg/(360 f),
%                    por interpolacion spline del canal (convenio de Codex y
%                    de robustness_sensor_errors_r12.m).
%   gain      [1x3]  error de ganancia relativo por canal (0.01 = +1 %).
%   snr_db    escalar; ruido blanco gaussiano independiente por canal con
%                    sigma = v_peak/sqrt(2) * 10^(-snr_db/20) (convenio de
%                    DIAG_R2_Paco.m: SNR referida al valor eficaz de fase en
%                    regimen sano).
%   v_peak    pico de fase sano; por defecto max|V| en t_ref.
%   t_ref     [t0 t1] tramo sano para v_peak; por defecto [2.4 2.8) si existe,
%                    si no todo el registro.
%   seed      semilla (obligatoria si hay ruido); se usa un RandStream propio,
%                    no altera el generador global.
%   f         50 Hz.
% info: struct con v_peak, sigma y los parametros aplicados.
p = struct('phase_deg', [0 0 0], 'gain', [0 0 0], 'snr_db', [], 'v_peak', [], ...
    't_ref', [2.4 2.8], 'seed', [], 'f', 50);
for k = 1:2:numel(varargin)
    if ~isfield(p, varargin{k}), error('r2_apply_meas_error:opt', 'Opcion desconocida: %s', varargin{k}); end
    p.(varargin{k}) = varargin{k+1};
end
t = t(:); Ve = V;
if isempty(p.v_peak)
    ref = t >= p.t_ref(1) & t < p.t_ref(2);
    if ~any(ref), ref = true(size(t)); end
    p.v_peak = max(abs(V(ref, :)), [], 'all');
end
for ph = 1:3
    if p.phase_deg(ph) ~= 0
        dt = p.phase_deg(ph)/360/p.f;
        Ve(:, ph) = interp1(t, V(:, ph), t - dt, 'spline', 'extrap');
    end
end
Ve = Ve .* (1 + p.gain(:)');
sigma = 0;
if ~isempty(p.snr_db) && isfinite(p.snr_db)
    if isempty(p.seed), error('r2_apply_meas_error:seed', 'El ruido exige una semilla.'); end
    rs = RandStream('mt19937ar', 'Seed', p.seed);
    sigma = p.v_peak/sqrt(2) * 10^(-p.snr_db/20);
    Ve = Ve + sigma*randn(rs, size(Ve));
end
info = p; info.sigma = sigma;
end
