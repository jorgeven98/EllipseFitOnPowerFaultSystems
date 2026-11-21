%% V2 sostenida con gate estricto M=10 (cierre de la tabla de robustez)
clear; clc;
addpath(fullfile('..', 'src', 'fault_classifiers'));
addpath(fullfile('..', 'src', 'utils'));
fs=4000; f=50; W=20; STRIDE=5; M=10;
base = load(fullfile('..','data','simulation_curated','DataGridRL_AB.mat'));
[Vpex,tpex] = extract_converter_data(base,'Pex');
seg = tpex>=2.4 & tpex<2.8; Vh = Vpex(seg,:); th = tpex(seg); th=th-th(1);
Npeak = max(abs(Vh(:)));
for k2 = [0.005 0.01 0.02]
    v2 = k2*Npeak*[cos(2*pi*f*th), cos(2*pi*f*th+2*pi/3), cos(2*pi*f*th-2*pi/3)];
    lab = clsloc(Vh+v2, W, STRIDE, fs, f);
    fprintf('V2=%.1f%%: estricto M=%d -> %.2f%% falsa operación\n', 100*k2, M, 100*prate_strict(lab,M));
end
function labels = clsloc(V, W, stride, fs, f)
    n = size(V,1); labels = {};
    for s = 1:stride:(n - W + 1)
        try, ft = ga_fault_classifier(V(s:s+W-1,:), fs, f); catch, ft='Error'; end
        labels{end+1} = ft; %#ok<AGROW>
    end
end
function r = prate_strict(labels,M)
    n=numel(labels); fire=false(1,n);
    for i=M:n
        w=labels(i-M+1:i);
        if all(~strcmp(w,'A-B-C')) && numel(unique(w))==1, fire(i)=true; end
    end
    r = mean(fire);
end
