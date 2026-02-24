function [T, widthPx] = add_GC_normalization(T, RTAB384, stimNum)
% ADD_GC_NORMALIZATION
% Adds GC-normalized coordinates to table T using widthPx = RTAB384(stimNum,7).
%
% Requirements:
%   - T.assignment exists and contains "target","distractor","background" (and possibly others)
%   - For target/distractor RFs: T.along and T.perp_signed exist
%   - For background RFs: T.r_s exists
%
% Adds columns:
%   - T.along_GC
%   - T.perp_signed_GC
%   - T.r_s_GC
%
% Returns widthPx for reference.

widthPx = double(RTAB384(stimNum, 7));
if ~(isfinite(widthPx) && widthPx > 0)
    error('Invalid widthPx from RTAB384(%d,7): %g', stimNum, widthPx);
end

N = height(T);

% Ensure assignment is string for comparisons
assign = T.assignment;
if iscategorical(assign), assign = string(assign); end
if ischar(assign), assign = string(assign); end

isObj = (assign == "target") | (assign == "distractor");
isBg  = (assign == "background");

% Preallocate as NaN
along_GC = nan(N,1);
perp_GC  = nan(N,1);
r_s_GC   = nan(N,1);

% Object RFs
if any(isObj)
    if ~ismember("along", string(T.Properties.VariableNames))
        error('T.along is missing. Run add_arm_projection_metrics first.');
    end
    if ~ismember("perp_signed", string(T.Properties.VariableNames))
        error('T.perp_signed is missing. Run add_arm_projection_metrics first.');
    end
    along_GC(isObj) = T.along(isObj) ./ widthPx;
    perp_GC(isObj)  = T.perp_signed(isObj) ./ widthPx;
end

% Background RFs
if any(isBg)
    if ~ismember("r_s", string(T.Properties.VariableNames))
        error('T.r_s is missing. Run add_polar_about_s first.');
    end
    r_s_GC(isBg) = T.r_s(isBg) ./ widthPx;
end

% Attach to table
T.along_GC = along_GC;
T.perp_signed_GC = perp_GC;
T.r_s_GC = r_s_GC;

end
