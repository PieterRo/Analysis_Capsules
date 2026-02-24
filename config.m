function cfg = config()

cfg = struct();
cfgFile = mfilename('fullpath');
cfg.rootDir = fileparts(cfgFile);

% Environment-dependent raw data locations (priority order).
cfg.dataDirCandidates = { ...
    '/Users/pieter/Library/CloudStorage/Dropbox/Pieter/data', ...
    '/Users/pieter/Dropbox/Pieter/data' ...
};
idx = find(cellfun(@(p) exist(p, 'dir') == 7, cfg.dataDirCandidates), 1, 'first');
if isempty(idx)
    msg = sprintf('No valid data directory found. Checked:\\n  - %s', ...
        strjoin(cfg.dataDirCandidates, '\\n  - '));
    error('config:DataDirNotFound', '%s', msg);
end
cfg.dataDir = cfg.dataDirCandidates{idx};

projectRoot = fileparts(cfg.rootDir);
cfg.extrasDir  = fullfile(projectRoot, 'Extras');
cfg.logsDir    = fullfile(cfg.extrasDir, 'monkeyN', '_logs');
cfg.stimDir    = fullfile(cfg.extrasDir, 'monkeyN', 'comp1');
cfg.matDir     = fullfile(cfg.rootDir, 'data_mat');
cfg.resultsDir = fullfile(cfg.rootDir, 'results');

if ~exist(cfg.matDir,'dir')
    mkdir(cfg.matDir);
end

if ~exist(cfg.resultsDir,'dir')
    mkdir(cfg.resultsDir);
end

end
