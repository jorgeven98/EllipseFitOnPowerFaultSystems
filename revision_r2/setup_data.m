function setup_data()
% SETUP_DATA  Copies the public simulation records of the repository (../data/*.mat)
% to data/simulation_curated, the folder that the scripts of this revision expect.
% Run it once from any folder:  run('revision_r2/setup_data.m')
here = fileparts(mfilename('fullpath'));
src = fullfile(fileparts(here), 'data');
dst = fullfile(here, 'data', 'simulation_curated');
if ~exist(dst, 'dir'), mkdir(dst); end
f = dir(fullfile(src, 'DataGrid*.mat'));
for k = 1:numel(f)
    copyfile(fullfile(src, f(k).name), fullfile(dst, f(k).name));
end
fprintf('setup_data: %d files copied to %s\n', numel(f), dst);
end
