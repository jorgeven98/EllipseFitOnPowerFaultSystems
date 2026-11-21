function opts = r2_t20_defaults(varargin)
% R2_T20_DEFAULTS  Opciones comunes de la tarea T20 (dataset EMT de SINTEF y laboratorio).
%
%   opts = r2_t20_defaults()                       % valores por defecto
%   opts = r2_t20_defaults('tag', 'G4', 'chain', struct('kirchhoff_thr', 0.02))
%
% UNICO sitio donde se cambia la cadena de decision. opts.chain es un struct que se
% pasa integro a r2_decision_chain (cualquier opcion de T03_API.md: M, gate, Delta, K,
% kirchhoff_thr, cls_params, classifier, amb_mode...). Con chain vacio se usa la cadena
% de referencia de T03 tal cual. Si el clasificador o la puerta de secuencia cero
% cambian, basta con:
%
%   r2_t20_run_all('tag', 'variante', 'chain', struct('kirchhoff_thr', 0.02))
%
% y todo (CSV por registro, resumenes, tabla y cifras) se regenera en
% resultados\T20\variante\ sin tocar la linea base.
%
% Campos:
%   tag            ''  etiqueta de la variante; si no esta vacia las salidas van a
%                      resultados\T20\<tag>\ y los textos a textos\T20_<tag>_*.tex
%   chain          struct() opciones de r2_decision_chain
%   reuse_scan     false; con true guarda y reutiliza la clasificacion por ventana del
%                      dataset SINTEF (cache de unos 95 MB) si las opciones que obligan a
%                      reclasificar no han cambiado
%   sintef.base    carpeta del dataset (solo lectura)
%   sintef.scales  modos de escalado que se evaluan; el primero es el principal
%                      'pu_pre' pu con el pico de secuencia directa de prefalta de cada
%                               registro y bus; 'volts' sin escalar; 'pu_nom' base fija
%   sintef.v_nom   pico de fase de la base fija (400 V de linea)
%   sintef.limit   [] o numero maximo de registros (prueba de humo)
%   lab.signals    senales del laboratorio; la primera es la principal
%   lab.t_off_ms   ventana de evento tras la incepcion (como r2_lab_dy_check)
here = fileparts(mfilename('fullpath'));
opts.raiz = fileparts(here);
opts.tag = '';
opts.chain = struct();
opts.reuse_scan = false;
opts.sintef.base = getenv('SINTEF_DATA');   % carpeta del dataset de SINTEF (no publico; disponible bajo peticion)
opts.sintef.scales = {'pu_pre', 'volts', 'pu_nom'};
opts.sintef.v_nom = 400*sqrt(2)/sqrt(3);
opts.sintef.t_fault = 1.0;        % incepcion nominal del generador (s), mas el retardo del registro
opts.sintef.dur = 0.1;            % duracion nominal de la falta (s)
opts.sintef.limit = [];
opts.sintef.pick = {};            % lista de nombres de archivo (prueba de humo)
opts.lab.signals = {'conv1_ref_red', 'box2', 'conv1_fase'};
opts.lab.t_off_ms = 120;
opts.lab.pre_s = 1.5;
for k = 1:2:numel(varargin)
    name = varargin{k}; val = varargin{k+1};
    p = strsplit(name, '.');
    if numel(p) == 1, opts.(p{1}) = val; else, opts.(p{1}).(p{2}) = val; end
end
opts.outdir = fullfile(opts.raiz, '_revision_TPWRD', 'ronda2', 'resultados', 'T20');
if ~isempty(opts.tag), opts.outdir = fullfile(opts.outdir, opts.tag); end
if ~exist(opts.outdir, 'dir'), mkdir(opts.outdir); end
opts.textdir = fullfile(opts.raiz, '_revision_TPWRD', 'ronda2', 'textos');
% opciones resueltas de la cadena y clave de las que obligan a reclasificar
addpath(here);
d = r2_decision_chain('defaults'); c = opts.chain; fn = fieldnames(c);
for k = 1:numel(fn), d.(fn{k}) = c.(fn{k}); end
if isa(d.classifier, 'function_handle'), d.classifier = func2str(d.classifier); else, d.classifier = char(d.classifier); end
opts.chain_resolved = d;
% clave = todas las opciones resueltas menos las que solo afectan a la decision posterior
% (puerta de redondez, nivel ambiente, persistencia y evento): cualquier opcion nueva de la
% cadena entra en la clave sin tocar este archivo
dec = {'gate','Delta','K','thr_cap','gate_rc','amb_mode','r_amb','amb_init','amb_feed','amb_holdoff','M','healthy_label','truth','t_on','t_off','t_range','scale','fs'};
opts.scan_key = rmfield(d, intersect(fieldnames(d), dec));
end
