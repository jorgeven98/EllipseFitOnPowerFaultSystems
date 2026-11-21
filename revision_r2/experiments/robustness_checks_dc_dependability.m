%% Dependabilidad con DC: defensa combinada (LS-hist + persistencia) y gate relajado
% Gate estricto: M etiquetas de falta IDENTICAS consecutivas.
% Gate relajado: M consecutivas no-normales; etiqueta = voto mayoritario del bloque.
% Seguridad: se re-mide sobre control sano, cap-switching y DC sano.
clear; clc;
addpath(fullfile('..', 'src', 'fault_classifiers'));
addpath(fullfile('..', 'src', 'utils'));
fs=4000; f=50; W=20; STRIDE=5; M=10;
base = load(fullfile('..','data','simulation_curated','DataGridRL_AB.mat'));
[Vpex,tpex] = extract_converter_data(base,'Pex');
seg = tpex>=2.4 & tpex<2.8; Vh = Vpex(seg,:); th = tpex(seg); th=th-th(1);
Npeak = max(abs(Vh(:)));
segF = tpex>=3.00 & tpex<3.14; Vf = Vpex(segF,:); tf_=tpex(segF); tf_=tf_-tf_(1);
starts = 1:STRIDE:(size(Vf,1)-W+1); tcent = (starts+(W-1)/2)/fs;

fprintf('=== DEPENDABILIDAD: falta A-B real + DC, defensa LS-hist + persistencia M=%d ===\n', M);
for dc = [0 0.05 0.10 0.20 0.30]
    dcf = dc*Npeak*exp(-tf_/0.030);
    Vfdc = Vf; Vfdc(:,1)=Vfdc(:,1)+dcf;
    lab_raw = clsloc(Vfdc, W, STRIDE, fs, f, false);
    lab_def = clsloc(Vfdc, W, STRIDE, fs, f, true);
    [dS_raw, lS_raw] = fire_strict(lab_raw, tcent, M);
    [dS_def, lS_def] = fire_strict(lab_def, tcent, M);
    [dR_def, lR_def] = fire_majority(lab_def, tcent, M);
    fprintf('DC=%.2f: estricto crudo %6.2f ms (%s) | estricto+LS %6.2f ms (%s) | mayoritario+LS %6.2f ms (%s)\n', ...
        dc, 1e3*dS_raw, lS_raw, 1e3*dS_def, lS_def, 1e3*dR_def, lR_def);
end

fprintf('\n=== SEGURIDAD del gate mayoritario (M=%d) sobre perturbaciones sanas ===\n', M);
% control
lab0 = clsloc(Vh, W, STRIDE, fs, f, false);
fprintf('Control sano:            mayoritario %5.2f%% | estricto %5.2f%%\n', 100*prate_major(lab0,M), 100*prate_strict(lab0,M));
% DC sano 0.3
dcH = 0.30*Npeak*exp(-th/0.030); Vdc = Vh; Vdc(:,1)=Vdc(:,1)+dcH;
labD = clsloc(Vdc, W, STRIDE, fs, f, false);
fprintf('DC 0.3 pu sano:          mayoritario %5.2f%% | estricto %5.2f%%\n', 100*prate_major(labD,M), 100*prate_strict(labD,M));
% cap switching peor caso (400 Hz, 0.3)
tau_osc=0.010; t0=0.05;
osc = 0.30*Npeak*exp(-max(th-t0,0)/tau_osc).*sin(2*pi*400*(th-t0)).*(th>=t0);
[~,i0_]=min(abs(th-t0)); w0 = Vh(i0_,:)/Npeak;
Vcs = Vh; for ph=1:3, Vcs(:,ph)=Vcs(:,ph)+osc.*w0(ph); end
labC = clsloc(Vcs, W, STRIDE, fs, f, false);
fprintf('Cap-sw 400Hz 0.3 pu:     mayoritario %5.2f%% | estricto %5.2f%%\n', 100*prate_major(labC,M), 100*prate_strict(labC,M));
% V2 0.5%
v2 = 0.005*Npeak*[cos(2*pi*f*th), cos(2*pi*f*th+2*pi/3), cos(2*pi*f*th-2*pi/3)];
labV = clsloc(Vh+v2, W, STRIDE, fs, f, false);
fprintf('V2 0.5%% sostenida:       mayoritario %5.2f%% | estricto %5.2f%%\n', 100*prate_major(labV,M), 100*prate_strict(labV,M));

function labels = clsloc(V, W, stride, fs, f, dcrem)
    n = size(V,1); labels = {};
    H = 80; tH = (0:H-1)'/fs;
    X = [cos(2*pi*f*tH), sin(2*pi*f*tH), ones(H,1)]; pinvX = pinv(X);
    for s = 1:stride:(n - W + 1)
        e = s + W - 1; Vw = V(s:e,:);
        if dcrem && (e - H + 1) >= 1
            hist = V(e-H+1:e,:);
            for ph = 1:3
                coef = pinvX * hist(:,ph);
                Vw(:,ph) = Vw(:,ph) - coef(3);
            end
        end
        try, ft = ga_fault_classifier(Vw, fs, f); catch, ft='Error'; end
        labels{end+1} = ft; %#ok<AGROW>
    end
end
function [det, detlab] = fire_strict(labels, tcent, M)
    det = NaN; detlab = '--';
    for i=M:numel(labels)
        w = labels(i-M+1:i);
        if all(~strcmp(w,'A-B-C')) && numel(unique(w))==1, det=tcent(i); detlab=w{1}; return, end
    end
end
function [det, detlab] = fire_majority(labels, tcent, M)
    det = NaN; detlab = '--';
    for i=M:numel(labels)
        w = labels(i-M+1:i);
        if all(~strcmp(w,'A-B-C'))
            u = unique(w); cnt = cellfun(@(x) sum(strcmp(w,x)), u);
            [cmax, im] = max(cnt);
            if cmax > M/2, det=tcent(i); detlab=u{im}; return, end
        end
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
function r = prate_major(labels,M)
    n=numel(labels); fire=false(1,n);
    for i=M:n
        w=labels(i-M+1:i);
        if all(~strcmp(w,'A-B-C'))
            u = unique(w); cnt = cellfun(@(x) sum(strcmp(w,x)), u);
            if max(cnt) > M/2, fire(i)=true; end
        end
    end
    r = mean(fire);
end
