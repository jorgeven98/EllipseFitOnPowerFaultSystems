function E = r2_sector_event(S, truth, t_on, t_off, M)
% R2_SECTOR_EVENT  Metricas por evento sobre las etiquetas ventana a ventana (T21).
% Regla de decision: declaracion asegurada cuando las M ultimas etiquetas son
% identicas y distintas de 'A-B-C' (como prate_strict de
% robustness_sensor_errors_r12.m). Marcas temporales al final de la ventana.
% No incluye la puerta adaptativa de redondez: con senales sinteticas limpias el
% nivel ambiente es redondez 1 y su umbral, min(1-0.01, 0.99) = 0.99, coincide
% con el umbral fijo interno del clasificador.
lab = S.label; n = numel(lab);
fire = false(1,n); cls = cell(1,n); cls(:) = {''};
for i = M:n
    w = lab(i-M+1:i);
    if ~strcmp(w{1},'A-B-C') && all(strcmp(w, w{1})), fire(i) = true; cls{i} = w{1}; end
end
post = S.tend > t_on & S.tend <= t_off;          % ventanas que terminan dentro de la falta
full = S.tstart >= t_on & S.tend <= t_off;       % ventanas completamente en falta
ok = strcmp(lab, truth);
i1 = find(post & ok, 1);
E.t_first_correct_ms = NaN; if ~isempty(i1), E.t_first_correct_ms = 1e3*(S.tend(i1)-t_on); end
id = find(post & fire, 1);
E.declared = ~isempty(id);
E.declared_class = ''; E.t_declared_ms = NaN;
if E.declared, E.declared_class = cls{id}; E.t_declared_ms = 1e3*(S.tend(id)-t_on); end
ic = find(post & fire & strcmp(cls, truth), 1);
E.correct_secured = ~isempty(ic);
E.t_correct_secured_ms = NaN; if ~isempty(ic), E.t_correct_secured_ms = 1e3*(S.tend(ic)-t_on); end
% declaraciones erroneas aseguradas: rachas distintas con clase asegurada ~= truth
wrong = post & fire & ~strcmp(cls, truth);
E.n_wrong_secured = sum(diff([false wrong]) == 1);
wl = unique(cls(wrong)); E.wrong_classes = strjoin(wl, '|');
% racha correcta maxima dentro de la falta
r = 0; rmax = 0;
for i = find(post), if ok(i), r = r+1; rmax = max(rmax,r); else, r = 0; end, end
E.max_correct_run = rmax;
E.pct_full_windows_correct = 100*mean(ok(full));
E.pct_full_windows_healthy = 100*mean(strcmp(lab(full),'A-B-C'));
E.n_prefault_secured = sum(fire & S.tend <= t_on);
th = S.theta(full); E.theta_median_deg = NaN; E.round_median = NaN;
if any(full)
    mu = 0.5*rad2deg(angle(mean(exp(1j*2*deg2rad(th)))));
    E.theta_median_deg = mod(mu + median(mod(th-mu+90,180)-90), 180);
    E.round_median = median(S.round(full));
end
lf = lab(full); [u,~,j] = unique(lf); c = accumarray(j(:),1); [~,im] = max(c);
E.modal_label = u{im};
end
