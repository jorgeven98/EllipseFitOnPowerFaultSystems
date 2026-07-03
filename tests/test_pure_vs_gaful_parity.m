%% TEST: paridad ga_fault_classifier_pure vs ga_fault_classifier (GA-FuL)
% Compara etiqueta a etiqueta sobre el protocolo completo de la tabla de
% precisión del paper (12 ficheros, DG1+Pex, ventanas 1/4 y 1 ciclo dentro
% de la falta). Requiere GA-FuL instalado (si no, el original ya delega en
% el port puro y la comparación sería trivial).
% Resultado de referencia (jul-2026, R2024b): 3696/3696 ventanas idénticas.
addpath(fullfile('..', 'classifiers'));
addpath(fullfile('..', 'utils'));
assert(~isempty(which('gafulInit')), 'Este test requiere GA-FuL instalado.');
fs = 4000; f = 50;
files = {'DataGridRL_AB.mat','DataGridRL_AC.mat','DataGridRL_BC.mat','DataGrid_AG.mat', ...
         'DataGrid_BG.mat','DataGrid_CG.mat','DataGrid_ABG.mat','DataGrid_ACG.mat', ...
         'DataGrid_BCG.mat','DataGrid_AB.mat','DataGrid_BC.mat','DataGrid_AC.mat'};
n_tot = 0; n_diff = 0;
for fi = 1:numel(files)
    data = load(fullfile('..', 'data', files{fi}));
    for conv = {'DG1','Pex'}
        [V, t] = extract_converter_data(data, conv{1});
        i0 = find(t >= 3.0, 1); i1 = find(t <= 3.16, 1, 'last');
        for W = [20 80]
            stride = round(W*0.25);
            for s = i0:stride:(i1 - W + 1)
                Vw = V(s:s+W-1,:);
                l1 = ga_fault_classifier(Vw, fs, f);
                l2 = ga_fault_classifier_pure(Vw, fs, f);
                n_tot = n_tot + 1;
                if ~strcmp(l1, l2), n_diff = n_diff + 1; end
            end
        end
    end
end
fprintf('PARIDAD: %d/%d ventanas idénticas, %d discrepancias\n', n_tot-n_diff, n_tot, n_diff);
assert(n_diff == 0, 'El port puro difiere del original en %d ventanas', n_diff);
fprintf('test_pure_vs_gaful_parity: OK\n');
