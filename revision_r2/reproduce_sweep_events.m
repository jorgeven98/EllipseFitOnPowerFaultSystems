function reproduce_sweep_events(casedir, work)
% REPRODUCE_SWEEP_EVENTS  Recomputes the 170 events of the fault-resistance and
% inception-angle sweeps of the EMT campaign (Tables S11 to S13 of supplementary.pdf)
% with the files of the decision chain alone, and compares the outcome, the declared
% class and the latency of each event with reference_results/S2_final.csv.
%
%   reproduce_sweep_events(casedir)         casedir: folder with the EMT cases (S2a_*.mat, ...)
%   reproduce_sweep_events(casedir, work)   work: empty folder for the copy of the chain files
%
% The files that the decision chain needs are found with the dependency analysis of
% MATLAB, copied to an empty folder, and the MATLAB path is reset before the run, so
% that nothing else of this repository or of the local installation is used.
% The EMT cases are not in the repository; they are available from the authors on request.
root = fileparts(mfilename('fullpath'));
if nargin < 2 || isempty(work), work = tempname; end
if ~exist(work, 'dir'), mkdir(work); end
s2csv = fullfile(root, 'reference_results', 'S2_final.csv');
addpath(fullfile(root, 'experiments'), fullfile(root, 'src', 'fault_classifiers'), fullfile(root, 'src', 'utils'));
ent = {'r2_final_opts', 'r2_decision_chain', 'extract_converter_data'}; fl = {};
for k = 1:numel(ent), fl = union(fl, matlab.codetools.requiredFilesAndProducts(which(ent{k}))); end
fl = fl(startsWith(fl, root));   % files of this repository only
for k = 1:numel(fl), copyfile(fl{k}, work); fprintf('chain file: %s\n', strrep(fl{k}, [root filesep], '')); end
restoredefaultpath; addpath(work); w = warning('off', 'all');
o = r2_final_opts(); warning(w);
T = readtable(s2csv, 'TextType', 'char'); n = height(T); same = 0; dmax = 0; bad = {};
for k = 1:n
    Sx = load(fullfile(casedir, [T.id{k} '.mat'])); c = Sx.caso;
    [V, t] = extract_converter_data(Sx, T.terminal{k}); V = double(V); t = double(t(:));
    t_on = c.t_falta; t_off = t_on + c.dur_falta;
    out = r2_decision_chain(V, t, o, 'truth', c.tipo, 't_on', t_on, 't_off', t_off, 'f', c.red_f, ...
        't_range', [t_on - 1.0, min(t(end), t_off + 0.12)]);
    e = out.event; ok = strcmp(char(e.outcome), T.outcome{k}) && strcmp(char(e.decl_class), T.decl_class{k});
    d = abs(e.decl_ms - T.decl_ms(k)); if isnan(e.decl_ms) && isnan(T.decl_ms(k)), d = 0; end
    if ok && d < 1e-6
        same = same + 1;
    else
        bad{end+1} = sprintf('%s %s: %s/%s/%.2f, reference %s/%s/%.2f', T.id{k}, T.terminal{k}, ...
            char(e.outcome), char(e.decl_class), e.decl_ms, T.outcome{k}, T.decl_class{k}, T.decl_ms(k)); %#ok<AGROW>
    end
    if ~isnan(d), dmax = max(dmax, d); end
end
fprintf('events: %d | identical outcome, class and latency: %d | largest latency difference: %.3g ms\n', n, same, dmax);
for k = 1:numel(bad), fprintf('  DIFF %s\n', bad{k}); end
fprintf('MATLAB %s on %s\n', version, computer);
end
