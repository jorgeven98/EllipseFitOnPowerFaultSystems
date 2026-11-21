%% Chequeos residuales de la suite de robustez (900 Hz M=12; declaración con DC)
clear; clc;
addpath(fullfile('..', 'src', 'fault_classifiers'));
addpath(fullfile('..', 'src', 'utils'));
fs=4000; f=50; W=20; STRIDE=5;
base = load(fullfile('..','data','simulation_curated','DataGridRL_AB.mat'));
[Vpex,tpex] = extract_converter_data(base,'Pex');
seg = tpex>=2.4 & tpex<2.8; Vh = Vpex(seg,:); th = tpex(seg); th = th-th(1);
Npeak = max(abs(Vh(:)));
segF = tpex>=3.00 & tpex<3.14; Vf = Vpex(segF,:); tf_ = tpex(segF); tf_=tf_-tf_(1);

% 1) Residuo 900 Hz A=0.1 con M=8/10/12
tau_osc=0.010; t0=0.05;
osc = 0.10*Npeak*exp(-max(th-t0,0)/tau_osc).*sin(2*pi*900*(th-t0)).*(th>=t0);
[~,i0_]=min(abs(th-t0)); w0 = Vh(i0_,:)/Npeak;
Vcs = Vh; for ph=1:3, Vcs(:,ph)=Vcs(:,ph)+osc.*w0(ph); end
lab = clsloc(Vcs, W, STRIDE, fs, f);
for M = [8 10 12]
    fprintf('900Hz A=0.1: persistencia M=%d -> %.2f%%\n', M, prate(lab,M)*100);
end

% 2) Falta real A-B con DC 0.3 pu bajo persistencia M=8
tau=0.030; dcf = 0.30*Npeak*exp(-tf_/tau);
Vfdc = Vf; Vfdc(:,1)=Vfdc(:,1)+dcf;
labf = clsloc(Vfdc, W, STRIDE, fs, f);
starts = 1:STRIDE:(size(Vfdc,1)-W+1); tcent = (starts+(W-1)/2)/fs;
[det, detlab] = firstfire(labf, tcent, 8);
fprintf('Falta A-B + DC 0.3pu: declaración M=8 a %.2f ms, etiqueta %s\n', 1e3*det, detlab);
labf0 = clsloc(Vf, W, STRIDE, fs, f);
[det0, detlab0] = firstfire(labf0, tcent, 8);
fprintf('Falta A-B sin DC:     declaración M=8 a %.2f ms, etiqueta %s\n', 1e3*det0, detlab0);

function labels = clsloc(V, W, stride, fs, f)
    n = size(V,1); labels = {};
    for s = 1:stride:(n - W + 1)
        try, ft = ga_fault_classifier(V(s:s+W-1,:), fs, f); catch, ft='Error'; end
        labels{end+1} = ft; %#ok<AGROW>
    end
end
function r = prate(labels,M)
    n=numel(labels); fire=false(1,n);
    for i=M:n
        w=labels(i-M+1:i);
        if all(~strcmp(w,'A-B-C')) && numel(unique(w))==1, fire(i)=true; end
    end
    r = mean(fire);
end
function [det, detlab] = firstfire(labels, tcent, M)
    det = NaN; detlab = '';
    for i=M:numel(labels)
        w = labels(i-M+1:i);
        if all(~strcmp(w,'A-B-C')) && numel(unique(w))==1
            det = tcent(i); detlab = w{1}; return
        end
    end
end
