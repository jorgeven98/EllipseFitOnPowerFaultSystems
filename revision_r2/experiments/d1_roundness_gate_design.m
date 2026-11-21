%% D1 - Diseno del gate de roundness (rama L-L vs normal) ante desbalance ambiente V2
% Decision D1 de la revision TPWRD.
%
% Problema: ga_fault_classifier desvia a la rama "circulo/normal" solo si
% roundness = b/a > 0.99. Con desbalance ambiente sostenido V2 >= ~0.5% eso
% produce etiquetas L-L falsas persistentes (V2=1% -> ~10.7% de falsa operacion
% asegurada con M=10; V2=2% -> ~55.8%). Bajar el umbral fijo a ~0.91 mata la
% dependabilidad en DG1 (grid-forming): sus faltas L-L reales tienen r ~ 0.95.
%
% Este script:
%   FASE 1  una sola pasada del clasificador por senal (cacheada) guardando por
%           ventana: etiqueta original, roundness, angle_deg, semiejes.
%   PARTE A caracterizacion de la distribucion de roundness (sano vs falta).
%   PARTE B barrido de gate FIJO r_c.
%   PARTE C gate ADAPTATIVO: thr = min(r_amb - Delta, 0.99), con r_amb = mediana
%           del roundness del regimen sano del sitio (K=80 ventanas).
%   PARTE D suelo de detectabilidad del diseno ganador (faltas A-B sinteticas).
%
% La etiqueta re-gateada se reconstruye desde (roundness, angle_deg): las
% etiquetas con tierra {A-G,B-G,C-G,A-B-G,B-C-G,C-A-G} no se tocan; para el
% resto: roundness < thr -> L-L por sector del angulo (A-B 45+-30, C-A 105+-30,
% B-C 165+-30, mod 180); si no -> A-B-C. (No sirve la etiqueta devuelta por el
% clasificador para la rama L-L/normal: ya lleva aplicado su 0.99 interno.)
%
% Salidas: ../_revision_TPWRD/robustness/d1_parte{A,B,C,D}_*.csv
% Cache:   ../_revision_TPWRD/robustness/d1_scan_cache.mat (borrar para re-escanear)

clear; clc;
addpath(fullfile('..','src','fault_classifiers'));
addpath(fullfile('..','src','utils'));

fs=4000; f=50; W=20; STRIDE=5; M=10; K_AMB=80;
OUTDIR = fullfile('..','_revision_TPWRD','robustness');
CACHE  = fullfile(OUTDIR,'d1_scan_cache.mat');
GROUND = {'A-G','B-G','C-G','A-B-G','B-C-G','C-A-G'};
K2H    = [0 0.005 0.01 0.02];      % V2 superpuesta al sano
K2F    = [0 0.01];                 % V2 superpuesta a la falta
CONVS  = {'Pex','DG1'};
FFILES = {'DataGridRL_AB.mat','A-B'; 'DataGridRL_BC.mat','B-C'; 'DataGridRL_AC.mat','C-A'};
RC     = [0.90 0.92 0.94 0.95 0.96 0.97 0.98 0.99 0.995];   % gates fijos
DELTAS = [0.01 0.02 0.03 0.05];                             % gates adaptativos
KSEV   = [0.01 0.02 0.03 0.04 0.05 0.10 0.15 0.20 0.30];    % severidades sinteticas (extendida por debajo de 0.05 para localizar el suelo)
T_ON=3.0; T_OFF=3.16;              % falta real
T0S=0.5;                           % onset sintetico

if exist(CACHE,'file'), tmp=load(CACHE); scans=tmp.scans; else, scans=struct(); end

%% ================= FASE 1: escaneos (1 pasada del clasificador por senal) ====
tic;
% --- (1) sano x V2: segmento t=[2.4,2.8) del fichero AB, Pex y DG1
dAB = load(fullfile('..','data','simulation_curated','DataGridRL_AB.mat'));
for c=1:numel(CONVS)
    cv=CONVS{c};
    [V,t]=extract_converter_data(dAB,cv);
    seg = t>=2.4 & t<2.8; Vh=V(seg,:); th=t(seg);
    npk = max(abs(Vh),[],'all');
    for k2=K2H
        key=sprintf('H_%s_k%03d',cv,round(1000*k2));
        if ~isfield(scans,key)
            scans.(key)=scanwin(Vh+v2wave(k2,npk,f,th), th, W,STRIDE,fs,f);
            save(CACHE,'scans'); fprintf('[scan] %-16s %3d ventanas (t=%.0fs)\n',key,numel(scans.(key).label),toc);
        end
    end
end

% --- (2) faltas reales: ventanas con centro en [3.0,3.16], con y sin V2=1%
for c=1:numel(CONVS)
    cv=CONVS{c};
    for q=1:size(FFILES,1)
        need=false;
        for k2=K2F
            if ~isfield(scans,fkey(cv,FFILES{q,2},k2)), need=true; end
        end
        if ~need, continue; end
        dq = load(fullfile('..','data','simulation_curated',FFILES{q,1}));
        [V,t]=extract_converter_data(dq,cv);
        npk = max(abs(V(t>=2.4&t<2.8,:)),[],'all');
        seg = t>=2.99 & t<=3.18; Vf=V(seg,:); tf_=t(seg);
        for k2=K2F
            key=fkey(cv,FFILES{q,2},k2);
            if isfield(scans,key), continue; end
            S=scanwin(Vf+v2wave(k2,npk,f,tf_), tf_, W,STRIDE,fs,f);
            scans.(key)=subsel(S, S.tc>=T_ON & S.tc<=T_OFF);
            save(CACHE,'scans'); fprintf('[scan] %-16s %3d ventanas (t=%.0fs)\n',key,numel(scans.(key).label),toc);
        end
    end
end

% --- (3) faltas A-B sinteticas de severidad k sobre fondo V2
tsyn=(0:4399)'/fs;   % 1.1 s; onset en t0=0.5 s
vbase=[cos(2*pi*f*tsyn), cos(2*pi*f*tsyn-2*pi/3), cos(2*pi*f*tsyn+2*pi/3)]; % secuencia abc (como los datos)
for k=KSEV
    for k2=K2F
        key=sprintf('S_k%03d_v%03d',round(1000*k),round(1000*k2));
        if isfield(scans,key), continue; end
        V=vbase; dv=(k/2)*(V(:,1)-V(:,2)); msk=tsyn>=T0S;
        V(msk,1)=V(msk,1)-dv(msk); V(msk,2)=V(msk,2)+dv(msk);
        scans.(key)=scanwin(V+v2wave(k2,1,f,tsyn), tsyn, W,STRIDE,fs,f);
        save(CACHE,'scans'); fprintf('[scan] %-16s %3d ventanas (t=%.0fs)\n',key,numel(scans.(key).label),toc);
    end
end
fprintf('FASE 1 completa (%.0f s)\n\n',toc);

% Registro de trazas
HT={};                             % sanas: {key, cv, k2}
for c=1:numel(CONVS), for k2=K2H, HT(end+1,:)={sprintf('H_%s_k%03d',CONVS{c},round(1000*k2)),CONVS{c},k2}; end, end %#ok<*SAGROW,*AGROW>
FT={};                             % faltas: {key, cv, truth, k2}
for c=1:numel(CONVS), for q=1:size(FFILES,1), for k2=K2F, FT(end+1,:)={fkey(CONVS{c},FFILES{q,2},k2),CONVS{c},FFILES{q,2},k2}; end, end, end

%% ================= PARTE A: caracterizacion del roundness =====================
fprintf('=== PARTE A: distribucion de roundness por ventana ===\n');
fprintf('%-18s %-5s %5s | %7s %7s %7s %7s | ang_med  %%ground\n','traza','tipo','nwin','mediana','p5','p95','min');
rowsA={};
for i=1:size(HT,1)
    S=scans.(HT{i,1}); r=S.roundness;
    rowsA(end+1,:)=stataRow(HT{i,1},'sano',HT{i,2},HT{i,3},'',S,true(numel(r),1),GROUND);
end
for i=1:size(FT,1)
    S=scans.(FT{i,1});
    full = S.tstart>=T_ON & S.tend<=T_OFF;    % ventanas totalmente dentro de la falta
    rowsA(end+1,:)=stataRow(FT{i,1},'falta',FT{i,2},FT{i,4},FT{i,3},S,full,GROUND);
end
for i=1:size(rowsA,1)
    fprintf('%-18s %-5s %5d | %7.4f %7.4f %7.4f %7.4f | %7.1f  %5.1f\n',rowsA{i,1},rowsA{i,2},rowsA{i,7},rowsA{i,8},rowsA{i,9},rowsA{i,10},rowsA{i,11},rowsA{i,12},rowsA{i,13});
end
TA=cell2table(rowsA,'VariableNames',{'traza','tipo','convertidor','v2_pct','truth','pendiente','n_ventanas','r_mediana','r_p5','r_p95','r_min','angulo_mediana','pct_ventanas_ground'});
TA.pendiente=[]; writetable(TA,fullfile(OUTDIR,'d1_parteA_roundness.csv'));

%% ================= PARTE B: gate FIJO r_c =====================================
fprintf('\n=== PARTE B: gate fijo (L-L si roundness < r_c) | falsa op. asegurada M=%d ===\n',M);
rowsB={};
for rc=RC
    fo=zeros(1,size(HT,1));
    for i=1:size(HT,1)
        S=scans.(HT{i,1});
        lab=regate(S,rc*ones(numel(S.roundness),1),GROUND);
        fo(i)=100*assured_rate(lab,M);
        rowsB(end+1,:)={rc,'sano',HT{i,1},HT{i,2},100*HT{i,3},'',fo(i),NaN,NaN,NaN,NaN};
    end
    nok=0; latmax=0; retmin=100;
    for i=1:size(FT,1)
        S=scans.(FT{i,1});
        lab=regate(S,rc*ones(numel(S.roundness),1),GROUND);
        [dec,corr,lat,~]=declare(lab,S.tend,M,FT{i,3},T_ON);
        full=S.tstart>=T_ON & S.tend<=T_OFF;
        ret=100*mean(strcmp(lab(full),FT{i,3}));
        nok=nok+double(corr); if corr, latmax=max(latmax,lat); end; retmin=min(retmin,ret);
        rowsB(end+1,:)={rc,'falta',FT{i,1},FT{i,2},100*FT{i,4},FT{i,3},NaN,double(dec),double(corr),lat,ret};
    end
    fprintf('r_c=%.3f | FO%% sano [Pex 0/0.5/1/2%%: %5.1f %5.1f %5.1f %5.1f | DG1: %5.1f %5.1f %5.1f %5.1f] | faltas OK %2d/12 latmax %5.1f ms retmin %5.1f%%\n',...
        rc,fo(1),fo(2),fo(3),fo(4),fo(5),fo(6),fo(7),fo(8),nok,latmax,retmin);
end
TB=cell2table(rowsB,'VariableNames',{'r_c','tipo','traza','convertidor','v2_pct','truth','falsa_op_pct','declarada','correcta','latencia_ms','retencion_pct'});
writetable(TB,fullfile(OUTDIR,'d1_parteB_gate_fijo.csv'));

%% ================= PARTE C: gate ADAPTATIVO thr = min(r_amb - D, 0.99) ========
fprintf('\n=== PARTE C: gate adaptativo (r_amb = mediana roundness sano del sitio, K=%d) ===\n',K_AMB);
rowsC={}; resC=struct();
for dd=1:numel(DELTAS)
    D=DELTAS(dd);
    fo=zeros(1,size(HT,1)); marg=zeros(1,size(HT,1));
    for i=1:size(HT,1)
        S=scans.(HT{i,1});
        ramb=median(S.roundness,'omitnan'); thr=min(ramb-D,0.99);
        lab=regate(S,thr*ones(numel(S.roundness),1),GROUND);
        fo(i)=100*assured_rate(lab,M);
        marg(i)=pct(S.roundness,5)-thr;   % margen del p5 ambiente sobre el umbral
        rowsC(end+1,:)={D,'sano',HT{i,1},HT{i,2},100*HT{i,3},'',ramb,thr,marg(i),fo(i),NaN,NaN,NaN,NaN};
    end
    nok=0; latmax=0; latsum=0; retmin=100;
    for i=1:size(FT,1)
        S=scans.(FT{i,1});
        hk=sprintf('H_%s_k%03d',FT{i,2},round(1000*FT{i,4}));   % sano del mismo convertidor y misma V2
        ramb=median(scans.(hk).roundness,'omitnan'); thr=min(ramb-D,0.99);
        lab=regate(S,thr*ones(numel(S.roundness),1),GROUND);
        [dec,corr,lat,~]=declare(lab,S.tend,M,FT{i,3},T_ON);
        full=S.tstart>=T_ON & S.tend<=T_OFF;
        ret=100*mean(strcmp(lab(full),FT{i,3}));
        nok=nok+double(corr); if corr, latmax=max(latmax,lat); latsum=latsum+lat; end; retmin=min(retmin,ret);
        rowsC(end+1,:)={D,'falta',FT{i,1},FT{i,2},100*FT{i,4},FT{i,3},ramb,thr,NaN,NaN,double(dec),double(corr),lat,ret};
    end
    resC(dd).D=D; resC(dd).foMax=max(fo); resC(dd).nok=nok; resC(dd).latmax=latmax;
    resC(dd).latmean=latsum/max(nok,1); resC(dd).margmin=min(marg); resC(dd).retmin=retmin;
    fprintf('Delta=%.2f | FO%% max sano %5.1f | margen p5-thr min %+.4f | faltas OK %2d/12 latmax %5.1f ms latmed %5.1f ms retmin %5.1f%%\n',...
        D,max(fo),min(marg),nok,latmax,latsum/max(nok,1),retmin);
end
TC=cell2table(rowsC,'VariableNames',{'delta','tipo','traza','convertidor','v2_pct','truth','r_amb','thr','margen_p5','falsa_op_pct','declarada','correcta','latencia_ms','retencion_pct'});
writetable(TC,fullfile(OUTDIR,'d1_parteC_gate_adaptativo.csv'));

% Eleccion de Delta*: seguridad perfecta (FO=0 hasta V2=2%) + 12/12 faltas correctas;
% entre candidatos, menor latencia maxima; empate -> mayor margen ambiente.
ok=find([resC.foMax]==0 & [resC.nok]==12);
if isempty(ok)
    [~,ib]=min([resC.foMax]*1e3-[resC.nok]); Dstar=resC(ib).D;
    fprintf('AVISO: ningun Delta cumple seguridad+dependabilidad perfectas; mejor compromiso Delta=%.2f\n',Dstar);
else
    cand=resC(ok); [~,ord]=sortrows([[cand.latmax]' , -[cand.margmin]'],[1 2]);
    Dstar=cand(ord(1)).D;
    fprintf('Delta* elegido = %.2f (candidatos perfectos: %s)\n',Dstar,mat2str([cand.D]));
end

%% ================= PARTE D: suelo de detectabilidad (diseno ganador) ==========
fprintf('\n=== PARTE D: faltas A-B sinteticas, gate adaptativo (todas las Deltas; ganador Delta*=%.2f) ===\n',Dstar);
rowsD={};
for dd=1:numel(DELTAS)
    D=DELTAS(dd);
    for k2=K2F
        kmin=NaN;
        for k=KSEV
            S=scans.(sprintf('S_k%03d_v%03d',round(1000*k),round(1000*k2)));
            pre=find(S.tend<=T0S); pre=pre(max(1,end-K_AMB+1):end);   % ultimas K ventanas sanas
            ramb=median(S.roundness(pre),'omitnan'); thr=min(ramb-D,0.99);
            lab=regate(S,thr*ones(numel(S.roundness),1),GROUND);
            [dec,corr,lat,~]=declare(lab,S.tend,M,'A-B',T0S);
            full=S.tstart>=T0S;
            ret=100*mean(strcmp(lab(full),'A-B'));
            rmed=median(S.roundness(full),'omitnan');
            if corr && isnan(kmin), kmin=k; end
            rowsD(end+1,:)={D,double(D==Dstar),k,100*k2,ramb,thr,rmed,double(dec),double(corr),lat,ret};
            if D==Dstar
                fprintf('Delta*=%.2f V2=%.0f%% k=%.2f | r_amb=%.4f thr=%.4f r_falta=%.4f | %s lat %6.1f ms ret %5.1f%%\n',...
                    D,100*k2,k,ramb,thr,rmed,okstr(dec,corr),lat,ret);
            end
        end
        if D==Dstar, fprintf('  -> suelo detectable (V2=%.0f%%): k_min = %s\n',100*k2,fmtk(kmin)); end
    end
end
TD=cell2table(rowsD,'VariableNames',{'delta','es_ganador','k_severidad','v2_pct','r_amb','thr','r_falta_mediana','declarada','correcta','latencia_ms','retencion_pct'});
writetable(TD,fullfile(OUTDIR,'d1_parteD_detectabilidad.csv'));
fprintf('\nCSVs escritos en %s\n',OUTDIR);

%% ================= ANEXO: limite DG1 con V2=2% (fuera del criterio D1, para acotar el claim)
% El suelo de la parte D predice que la falta L-L de un GFM (NS ~2.7%) se pierde con
% ambiente V2=2%. Se comprueba con las faltas reales de DG1 + V2=2% superpuesta.
fprintf('\n=== ANEXO: faltas DG1 con V2=2%% ambiente, gate adaptativo Delta*=%.2f ===\n',Dstar);
for q=1:size(FFILES,1)
    key=fkey('DG1',FFILES{q,2},0.02);
    if ~isfield(scans,key)
        dq=load(fullfile('..','data','simulation_curated',FFILES{q,1}));
        [V,t]=extract_converter_data(dq,'DG1');
        npk=max(abs(V(t>=2.4&t<2.8,:)),[],'all');
        seg=t>=2.99 & t<=3.18; Vf=V(seg,:); tf_=t(seg);
        S=scanwin(Vf+v2wave(0.02,npk,f,tf_), tf_, W,STRIDE,fs,f);
        scans.(key)=subsel(S, S.tc>=T_ON & S.tc<=T_OFF);
        save(CACHE,'scans');
    end
    S=scans.(key);
    ramb=median(scans.H_DG1_k020.roundness,'omitnan'); thr=min(ramb-Dstar,0.99);
    lab=regate(S,thr*ones(numel(S.roundness),1),GROUND);
    [dec,corr,lat,~]=declare(lab,S.tend,M,FFILES{q,2},T_ON);
    full=S.tstart>=T_ON & S.tend<=T_OFF;
    fprintf('DG1 %s + V2=2%%: r_falta med=%.4f vs thr=%.4f | %s lat %6.1f ms ret %5.1f%%\n',...
        FFILES{q,2},median(S.roundness(full),'omitnan'),thr,okstr(dec,corr),lat,100*mean(strcmp(lab(full),FFILES{q,2})));
end

%% ================= FUNCIONES LOCALES ==========================================
function key=fkey(cv,truth,k2)
    key=sprintf('F_%s_%s_k%03d',cv,strrep(truth,'-',''),round(1000*k2));
end

function v2=v2wave(k2,npk,f,t)
% Secuencia inversa sostenida de amplitud k2*npk (t absoluto; fase irrelevante al ser sostenida)
    v2=k2*npk*[cos(2*pi*f*t), cos(2*pi*f*t+2*pi/3), cos(2*pi*f*t-2*pi/3)];
end

function S=scanwin(V,t,W,STRIDE,fs,f)
% Una pasada del clasificador con ventana deslizante; guarda trazas por ventana.
    n=size(V,1); starts=1:STRIDE:(n-W+1); m=numel(starts);
    S.tstart=zeros(m,1); S.tc=zeros(m,1); S.tend=zeros(m,1);
    S.label=cell(m,1); S.roundness=zeros(m,1); S.angle=zeros(m,1); S.a=zeros(m,1); S.b=zeros(m,1);
    for j=1:m
        s=starts(j); e=s+W-1;
        S.tstart(j)=t(s); S.tend(j)=t(e); S.tc(j)=(t(s)+t(e))/2;
        try
            [ft,p]=ga_fault_classifier(V(s:e,:),fs,f);
            S.label{j}=char(ft); S.roundness(j)=p.roundness; S.angle(j)=p.angle_deg;
            S.a(j)=p.semi_major; S.b(j)=p.semi_minor;
        catch
            S.label{j}='Error'; S.roundness(j)=NaN; S.angle(j)=NaN; S.a(j)=NaN; S.b(j)=NaN;
        end
    end
end

function S2=subsel(S,keep)
    fn=fieldnames(S);
    for q=1:numel(fn), S2.(fn{q})=S.(fn{q})(keep); end
end

function lab=regate(S,thr,GROUND)
% Reconstruye la etiqueta con el gate thr (vector por ventana):
% ground -> intocable; roundness<thr -> L-L por sector; si no -> A-B-C.
    m=numel(S.label); lab=cell(m,1);
    for j=1:m
        L=S.label{j};
        if any(strcmp(L,GROUND)) || strcmp(L,'Error')
            lab{j}=L;
        elseif S.roundness(j)<thr(j)
            lab{j}=sector_ll(S.angle(j));
        else
            lab{j}='A-B-C';
        end
    end
end

function s=sector_ll(ang)
% Sector L-L por angulo (mod 180): A-B 45, C-A 105, B-C 165 (+-30 tesela [0,180))
    cts=[45 105 165]; names={'A-B','C-A','B-C'};
    d=mod(ang-cts,180); d=min(d,180-d);
    [~,ix]=min(d); s=names{ix};
end

function r=assured_rate(lab,M)
% Fraccion de ventanas con declaracion asegurada (M etiquetas de falta identicas)
    n=numel(lab); fire=false(1,n);
    for i=M:n
        w=lab(i-M+1:i);
        if all(~strcmp(w,'A-B-C')) && all(~strcmp(w,'Error')) && numel(unique(w))==1, fire(i)=true; end
    end
    r=mean(fire);
end

function [dec,corr,lat_ms,declab]=declare(lab,tend,M,truth,t_on)
% Primera declaracion asegurada; latencia = fin de la M-esima ventana - onset
    dec=false; corr=false; lat_ms=NaN; declab='--';
    for i=M:numel(lab)
        w=lab(i-M+1:i);
        if all(~strcmp(w,'A-B-C')) && all(~strcmp(w,'Error')) && numel(unique(w))==1
            dec=true; declab=w{1}; corr=strcmp(declab,truth); lat_ms=1e3*(tend(i)-t_on); return
        end
    end
end

function v=pct(x,p)
% Percentil sin toolbox (interpolacion lineal sobre orden)
    x=sort(x(~isnan(x))); n=numel(x);
    if n==0, v=NaN; return, elseif n==1, v=x; return, end
    v=interp1((0:n-1)/(n-1),x,p/100,'linear');
end

function row=stataRow(name,tipo,cv,k2,truth,S,mask,GROUND)
    r=S.roundness(mask); a=S.angle(mask); L=S.label(mask);
    g=100*mean(ismember(L,GROUND));
    row={name,tipo,cv,100*k2,truth,'',numel(r),median(r,'omitnan'),pct(r,5),pct(r,95),min(r),median(a,'omitnan'),g};
end

function s=okstr(dec,corr)
    if ~dec, s='NO-DECL   '; elseif corr, s='OK        '; else, s='DECL-ERR  '; end
end

function s=fmtk(k)
    if isnan(k), s='ninguna (>0.30)'; else, s=sprintf('%.2f',k); end
end
