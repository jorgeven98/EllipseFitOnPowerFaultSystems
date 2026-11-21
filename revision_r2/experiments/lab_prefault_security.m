%% Seguridad pre-falta en los registros del lab: ¿dispara el gate M=10 en falso?
clear; clc;
addpath(fullfile('..', 'src', 'fault_classifiers'));
labdir = fullfile('..', 'data', 'lab', 'Fault_emulator_Outkontaktor');
files = dir(fullfile(labdir, '*.mat'));
fsL=10000; f=50; W=50; STRIDE=12; M=10; cyc=200; PERM=[3 2 1];
tot_win = 0; tot_fire = 0;
for k = 1:numel(files)
    fn = files(k).name;
    pref = regexp(fn, '^[a-z]+', 'match', 'once');
    if isempty(pref) || strcmp(pref,'abc'), continue; end
    d = load(fullfile(labdir, fn));
    R = d.results; t = R(1,:)'; t = t-t(1); V = R(59:61,:)'; V = V(:,PERM);
    nw = floor(numel(t)/cyc); dev = zeros(nw,1);
    for i=1:nw
        seg=(i-1)*cyc+1:i*cyc; r3=sqrt(mean(V(seg,:).^2,1)); dev(i)=max(r3)-min(r3);
    end
    base=median(dev); a0=find(dev>5*base+0.02,1);
    iEnd = max(1, (a0-3)*cyc);                    % hasta 3 ciclos antes de la falta
    iStart = max(1, iEnd - 10000);                % 1 s de pre-falta
    Vp = V(iStart:iEnd,:);
    lab = {};
    for s=1:STRIDE:(size(Vp,1)-W+1)
        try, ft=ga_fault_classifier(Vp(s:s+W-1,:),fsL,f); catch, ft='Err'; end
        lab{end+1}=ft; %#ok<AGROW>
    end
    fire = 0; n=numel(lab);
    for i=M:n
        w=lab(i-M+1:i);
        if all(~strcmp(w,'A-B-C')) && numel(unique(w))==1, fire=fire+1; end
    end
    tot_win = tot_win + n; tot_fire = tot_fire + fire;
    fprintf('%-28s pre-falta %4d ventanas, disparos falsos gate M=%d: %d\n', fn, n, M, fire);
end
fprintf('\nTOTAL: %d ventanas pre-falta, %d disparos falsos (%.2f%%)\n', tot_win, tot_fire, 100*tot_fire/max(tot_win,1));
