%% LAB_LATENCY_R3M1 — latencia de detección medida en registros del lab (OPAL)
% Para cada registro de falta L-L del banco (Fault_emulator_Outkontaktor):
%  1) localiza la región de falta (RMS por ciclo),
%  2) estima la incepción por cambio-punto: primer instante en que la señal
%     se desvía >5% del valor extrapolado del ajuste sinusoidal pre-falta,
%  3) clasifica con ga_fault_classifier (10 kHz, ventana 1/4 ciclo = 50 m.),
%  4) reporta: primera ventana correcta (fin de ventana) y declaración
%     asegurada (M=10 etiquetas idénticas consecutivas), relativas a la incepción.
% Canales: U1abc del sincronizador (filas 26-28 de results, 1-indexed 25:27+1).
clear; clc;
addpath(fullfile('..', 'src', 'fault_classifiers'));

labdir = fullfile('..', 'data', 'lab', 'Fault_emulator_Outkontaktor');
files = dir(fullfile(labdir, '*.mat'));
fsL = 10000; f = 50; W = 50; STRIDE = 12; M = 10;
CH = 59:61;                                   % V_abc_meas_pu del bus conv1
% El cableado del banco intercambia las fases a y c respecto al orden de
% columnas del clasificador (verificado con los 3 tipos de falta L-L):
PERM = [3 2 1];
expected_map = containers.Map({'ab','bc','ac','abc'}, {'A-B','B-C','C-A','A-B-C-G'});

fprintf('%-28s %-6s %-9s %-11s %-11s\n', 'registro', 'falta', 't_incep', '1a correcta', 'asegurada');
lat1 = []; latM = [];
for k = 1:numel(files)
    fn = files(k).name;
    pref = regexp(fn, '^[a-z]+', 'match', 'once');
    if isempty(pref) || ~isKey(expected_map, pref), continue; end
    expected = expected_map(pref);
    if strcmp(pref, 'abc'), continue; end     % simétrica: convención de etiqueta distinta, fuera del alcance L-L
    d = load(fullfile(labdir, fn));
    R = d.results; t = R(1,:)'; t = t - t(1); V = R(CH,:)'; V = V(:, PERM);
    Vdet = V;             % bus del convertidor: hundimiento claro, para LOCALIZAR la falta

    % 1) región de falta por RMS de 1 ciclo (sobre el bus del convertidor)
    cyc = fsL/f; nw = floor(numel(t)/cyc);
    dev = zeros(nw,1);
    for i=1:nw
        seg = (i-1)*cyc+1 : i*cyc;
        r3 = sqrt(mean(Vdet(seg,:).^2,1));
        dev(i) = max(r3) - min(r3);
    end
    base = median(dev); anom = find(dev > 5*base + 0.02, 1);
    if isempty(anom), fprintf('%-28s sin falta detectable\n', fn); continue; end
    i_region = (anom-2)*cyc;                  % empezar 2 ciclos antes

    % 2) incepción por cambio-punto con extrapolación sinusoidal
    % ajuste [cos,sin,1] sobre los 2 ciclos ANTERIORES a la región
    idx_fit = i_region-2*cyc : i_region-1;
    tt = t(idx_fit); X = [cos(2*pi*f*tt), sin(2*pi*f*tt), ones(numel(tt),1)];
    coefs = X \ V(idx_fit,:);
    idx_test = i_region : min(i_region+3*cyc, numel(t));
    Xt = [cos(2*pi*f*t(idx_test)), sin(2*pi*f*t(idx_test)), ones(numel(idx_test),1)];
    resid = abs(V(idx_test,:) - Xt*coefs);
    % primer instante con residuo >5% del pico nominal en cualquier fase,
    % sostenido 5 muestras (0.5 ms) para inmunidad a ruido
    Vpk = max(abs(V(idx_fit,:)),[],'all');
    hit = all(movmin(max(resid,[],2), [0 4]) > 0.05*Vpk, 2);
    ii = find(hit, 1);
    if isempty(ii), fprintf('%-28s incepción no detectada\n', fn); continue; end
    t_incep = t(idx_test(ii));

    % 3) clasificación en ventana deslizante desde 1 ciclo antes de la incepción
    i0 = find(t >= t_incep, 1) - cyc;
    i1 = min(find(t >= t_incep, 1) + 12*cyc, numel(t) - W);
    labels = {}; tend = [];
    for s = i0:STRIDE:i1
        try, ft = ga_fault_classifier(V(s:s+W-1,:), fsL, f); catch, ft = 'Error'; end
        labels{end+1} = ft; tend(end+1) = t(s+W-1); %#ok<AGROW>
    end
    % 4) latencias respecto a t_incep (solo ventanas que terminan tras la incepción)
    post = tend > t_incep;
    d1 = NaN; dM = NaN;
    for i = find(post)
        if strcmp(labels{i}, expected) && isnan(d1), d1 = tend(i) - t_incep; end
        if i >= M && all(strcmp(labels(i-M+1:i), expected)) && isnan(dM)
            dM = tend(i) - t_incep; break
        end
    end
    fprintf('%-28s %-6s %8.4f s %8.2f ms %8.2f ms\n', fn, expected, t_incep, 1e3*d1, 1e3*dM);
    lat1(end+1) = d1; latM(end+1) = dM;
end
lat1 = lat1(~isnan(lat1)); latM = latM(~isnan(latM));
fprintf('\nRESUMEN (n=%d): 1a correcta mediana %.2f ms (p90 %.2f) | asegurada M=%d mediana %.2f ms (p90 %.2f)\n', ...
    numel(lat1), 1e3*median(lat1), 1e3*prctile(lat1,90), M, 1e3*median(latM), 1e3*prctile(latM,90));

% distribución de etiquetas pre-falta (calidad ambiente del banco), último registro
fprintf('\nEtiquetas pre-falta (últimos 0.4 s antes de la región, último registro):\n');
segP = t >= t_incep-0.45 & t <= t_incep-0.05; Vp = V(segP,:); labP = {};
for s = 1:STRIDE:(size(Vp,1)-W+1)
    try, ft = ga_fault_classifier(Vp(s:s+W-1,:), fsL, f); catch, ft='Error'; end
    labP{end+1} = ft; %#ok<AGROW>
end
u = unique(labP);
for i=1:numel(u), fprintf('  %-8s: %d (%.0f%%)\n', u{i}, sum(strcmp(labP,u{i})), 100*mean(strcmp(labP,u{i}))); end
