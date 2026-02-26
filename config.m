function cfg = config()

cfg = struct();
cfgFile = mfilename('fullpath');
cfg.repoRoot = fileparts(cfgFile);

% Load machine-specific paths from config_local.m
localCfgPath = fullfile(cfg.repoRoot, 'config_local.m');
if exist(localCfgPath, 'file') == 2
    cfg = config_local(cfg);
else
    error('config:LocalConfigMissing', ...
        ['Missing config_local.m in repo root. Create it with:\n' ...
         'function cfg = config_local(cfg)\n' ...
         'cfg.dataRoot = ''/path/to/data'';\n' ...
         'cfg.extrasRoot = ''/path/to/Extras'';\n' ...
         'cfg.resultsRoot = ''/path/to/Dropbox'';\n' ...
         'end']);
end

if ~isfield(cfg, 'dataRoot') || isempty(cfg.dataRoot)
    error('config:DataRootMissing', ...
        'config_local.m must set cfg.dataRoot to your local data folder.');
end
if ~isfield(cfg, 'extrasRoot') || isempty(cfg.extrasRoot)
    error('config:ExtrasRootMissing', ...
        'config_local.m must set cfg.extrasRoot to your local Extras folder.');
end
if exist(cfg.dataRoot, 'dir') ~= 7
    error('config:DataRootNotFound', ...
        'Data root not found: %s', cfg.dataRoot);
end
if exist(cfg.extrasRoot, 'dir') ~= 7
    error('config:ExtrasRootNotFound', ...
        'Extras root not found: %s', cfg.extrasRoot);
end

cfg.logsDir = fullfile(cfg.extrasRoot, 'monkeyN', '_logs');
cfg.stimDir = fullfile(cfg.extrasRoot, 'monkeyN', 'comp1');

if ~isfield(cfg, 'matRoot') || isempty(cfg.matRoot)
    cfg.matRoot = cfg.repoRoot;
end
cfg.matDir = fullfile(cfg.matRoot, 'data_mat');

if ~isfield(cfg, 'resultsRoot') || isempty(cfg.resultsRoot)
    error('config:ResultsRootMissing', ...
        'config_local.m must set cfg.resultsRoot to your Dropbox root.');
end
cfg.resultsDir = fullfile(cfg.resultsRoot, 'results');

if exist(cfg.logsDir, 'dir') ~= 7
    error('config:LogsDirNotFound', ...
        'Expected logs directory not found: %s', cfg.logsDir);
end
if exist(cfg.stimDir, 'dir') ~= 7
    error('config:StimDirNotFound', ...
        'Expected stimulus directory not found: %s', cfg.stimDir);
end

if ~exist(cfg.matDir, 'dir')
    mkdir(cfg.matDir);
end
if ~exist(cfg.resultsDir, 'dir')
    mkdir(cfg.resultsDir);
end

end
