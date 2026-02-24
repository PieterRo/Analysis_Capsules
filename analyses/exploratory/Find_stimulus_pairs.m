cfg = config();
folder = cfg.stimDir;

% Pass 1: yellow is the anchor
R_yellow = findStimulusPairs_capsules_overlap_fast(folder, ...
    'Anchor','yellow', ...
    'AnchorIoUmin', 0.86, ...
    'OtherIoUmax',  0.60, ...
    'MaxShiftPx',   3);

pairs = R_yellow.pairs;
files = string(R.files);
disp([files(pairs(1:min(10,end),1)) files(pairs(1:min(10,end),2))])


% Pass 2: purple is the anchor
R_purple = findStimulusPairs_capsules_overlap_fast(folder, ...
    'Anchor','purple', ...
    'AnchorIoUmin', 0.86, ...
    'OtherIoUmax',  0.60, ...
    'MaxShiftPx',   3);



pairs = R_purple.pairs;
files = string(R.files);
disp([files(pairs(1:min(10,end),1)) files(pairs(1:min(10,end),2))])


