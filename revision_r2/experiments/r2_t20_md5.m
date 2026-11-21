function s = r2_t20_md5()
% R2_T20_MD5  Huella (MD5) de los archivos que fijan la decision en T20, para dejar
% constancia en los logs de con que version de la cadena se obtuvo cada resultado.
% Otros agentes amplian r2_decision_chain.m durante la ronda 2 (opciones nuevas con los
% mismos valores por defecto); los logs de T20 escriben esta huella al empezar y al
% terminar, y avisan si cambia durante la ejecucion.
here = fileparts(mfilename('fullpath')); src = fullfile(here, '..', 'src', 'fault_classifiers');
files = {fullfile(here, 'r2_decision_chain.m'), fullfile(here, 'r2_ga_classifier_param.m'), ...
         fullfile(src, 'ga_fault_classifier_pure.m'), fullfile(src, 'fit_ellipse_gac_improved.m')};
p = cell(1, numel(files));
for k = 1:numel(files)
    [~, n, e] = fileparts(files{k});
    if ~isfile(files{k}), p{k} = sprintf('%s%s=NO_EXISTE', n, e); continue, end
    fid = fopen(files{k}, 'r'); b = fread(fid, Inf, '*uint8'); fclose(fid);
    h = java.security.MessageDigest.getInstance('MD5'); h.update(b);
    p{k} = sprintf('%s%s=%s', n, e, sprintf('%02x', typecast(h.digest(), 'uint8')));
end
s = strjoin(p, ' ');
end
