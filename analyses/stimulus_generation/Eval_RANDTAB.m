% LogDir is defined in Line_Stimuli

RandFile1 = "RANDTAB_ObjAtt_lines_monkeyN_20220131_B1";
RandFile2 = "RANDTAB_ObjAtt_lines_monkeyN_20220131_B2";
RandFile3 = "RANDTAB_ObjAtt_lines_monkeyN_20220131_B3";
RandFile4 = "RANDTAB_ObjAtt_lines_monkeyN_20220131_B4";
RandFile5 = "RANDTAB_ObjAtt_lines_monkeyN_20220201_B1";
RandFile6 = "RANDTAB_ObjAtt_lines_monkeyN_20220201_B2";
RandFile7 = "RANDTAB_ObjAtt_lines_monkeyN_20220201_B3";
RandFile8 = "RANDTAB_ObjAtt_lines_monkeyN_20220201_B4";
%RandFile9 = "RANDTAB_ObjAtt_lines_monkeyN_20220204_B4";  % this one did
%not match the rest

% combine all RANDTABS in one, because single ones do not have all stimuli.
load(LogDir+"/"+RandFile1);  RTAB=RANDTAB;
load(LogDir+"/"+RandFile2);  RTAB=[RTAB;RANDTAB];
load(LogDir+"/"+RandFile3);  RTAB=[RTAB;RANDTAB];
load(LogDir+"/"+RandFile4);  RTAB=[RTAB;RANDTAB];
load(LogDir+"/"+RandFile5);  RTAB=[RTAB;RANDTAB];
load(LogDir+"/"+RandFile6);  RTAB=[RTAB;RANDTAB];
load(LogDir+"/"+RandFile7);  RTAB=[RTAB;RANDTAB];
load(LogDir+"/"+RandFile8);  RTAB=[RTAB;RANDTAB];
%load(LogDir+"/"+RandFile9);  RTAB=[RTAB;RANDTAB];

RTAB384 = nan(384, 8);   % output
badStim  = [];           % stimuli with mismatches
badCount = [];           % how many rows differ from the first

for s = 1:384
    idx = (RTAB(:,1) == s);
    if ~any(idx)
        % stimulus not present in RTAB; leave as NaNs (or error if you prefer)
        continue
    end

    rows = RTAB(idx, :);
    base = rows(1, :);
    RTAB384(s, :) = base;

    % check all rows match the first row
    mismatch = any(rows ~= base, 2);   % row-wise mismatch
    if any(mismatch)
        badStim(end+1)  = s; %#ok<AGROW>
        badCount(end+1) = sum(mismatch);
    end
end

% Report
if isempty(badStim)
    fprintf('All repeats are identical for all stimuli found.\n');
else
    fprintf('Mismatches found in %d stimuli:\n', numel(badStim));
    disp(table(badStim(:), badCount(:), 'VariableNames', {'stim','nMismatchingRows'}));
end

save(LogDir+"/"+"RTAB384.mat",'RTAB384');

