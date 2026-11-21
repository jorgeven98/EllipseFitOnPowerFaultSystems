function r2_t20_run_all(varargin)
% R2_T20_RUN_ALL  Regenera toda la tarea T20 con una sola llamada.
%
%   r2_t20_run_all()                                             % cadena de referencia
%   r2_t20_run_all('tag', 'G4', 'chain', struct('kirchhoff_thr', 0.02))
%   r2_t20_run_all('tag', 'M8', 'chain', struct('M', 8))         % reutiliza la clasificacion
%
% Los argumentos son los de r2_t20_defaults. 'chain' se pasa integro a
% r2_decision_chain (T03_API.md). Con tag vacio escribe en resultados\T20\ y en
% textos\T20_tabla.tex y T20_cifras.tex; con tag no vacio, en resultados\T20\<tag>\ y
% textos\T20_<tag>_*.tex, sin tocar la linea base.
%
% Pasos: 1) dataset EMT de SINTEF (r2_sintef_dataset, unos 12 min la primera vez);
%        2) registros de laboratorio (r2_lab_autokontaktor, unos 3 min);
%        3) tabla y cifras en LaTeX a partir de los CSV (r2_t20_tex).
% Ejecucion por linea de comandos, desde resultados\T20:
%   "C:\Program Files\MATLAB\R2024b\bin\matlab.exe" -batch "addpath('<raiz del proyecto>\experiments'); r2_t20_run_all"
opts = r2_t20_defaults(varargin{:});
r2_sintef_dataset(opts);
r2_lab_autokontaktor(opts);
r2_t20_tex(opts);
end
