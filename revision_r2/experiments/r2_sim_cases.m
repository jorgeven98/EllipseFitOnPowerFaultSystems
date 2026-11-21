function [T, info] = r2_sim_cases()
% R2_SIM_CASES  Los 12 registros de falta asimetrica de la simulacion del paper.
%
%   [T, info] = r2_sim_cases()
%
% T: tabla con archivo (en data/simulation_curated), clase verdadera con la
%    nomenclatura del clasificador y familia ('LL_RL', 'LG', 'LLG', 'LL').
%    Son los 12 archivos de experiments/compare_classifier_accuracy.m.
% info: constantes comunes comprobadas en la forma de onda (T03, nov-2026):
%    fs = 4000 Hz, f = 50 Hz; la falta se aplica en t_on = 3.000 s (la muestra
%    12001, t = 3.000, es la ultima de prefalta; la primera alterada es
%    t = 3.00025) y la orden de despeje llega en t_off = 3.160 s (el despeje
%    real ocurre entre 3.161 y 3.170 s, en el paso por cero de la corriente).
%    Regimen estacionario desde t = 0.1 s aprox.; t_range recomendado
%    [2.0 3.2]. Convertidores: 'Pex' (GFL) y 'DG1' (GFM). Unidades: voltios.
T = cell2table({
    'DataGridRL_AB.mat', 'A-B',   'LL_RL'
    'DataGridRL_AC.mat', 'C-A',   'LL_RL'
    'DataGridRL_BC.mat', 'B-C',   'LL_RL'
    'DataGrid_AG.mat',   'A-G',   'LG'
    'DataGrid_BG.mat',   'B-G',   'LG'
    'DataGrid_CG.mat',   'C-G',   'LG'
    'DataGrid_ABG.mat',  'A-B-G', 'LLG'
    'DataGrid_ACG.mat',  'C-A-G', 'LLG'
    'DataGrid_BCG.mat',  'B-C-G', 'LLG'
    'DataGrid_AB.mat',   'A-B',   'LL'
    'DataGrid_BC.mat',   'B-C',   'LL'
    'DataGrid_AC.mat',   'C-A',   'LL'}, 'VariableNames', {'file', 'truth', 'family'});
here = fileparts(mfilename('fullpath'));
info = struct('fs', 4000, 'f', 50, 't_on', 3.0, 't_off', 3.16, 't_range', [2.0 3.2], ...
    't_healthy', [2.4 2.8], 'converters', {{'Pex', 'DG1'}}, ...
    'dir', fullfile(here, '..', 'data', 'simulation_curated'));
end
