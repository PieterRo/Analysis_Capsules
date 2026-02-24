function Tall = build_all_stim_tables(ALLCOORDS, RTAB384, x_rf, y_rf)

nStim = size(RTAB384,1);
Tall = struct('stimNum', cell(nStim,1), ...
              'widthPx', cell(nStim,1), ...
              'T', cell(nStim,1));

fprintf('Processing %d stimuli...\n', nStim);

for stimNum = 1:nStim
    
    % ---- pipeline ----
    T = rf_table_target_distractor(ALLCOORDS, RTAB384, stimNum, x_rf, y_rf);
    T = add_arm_projection_metrics(T, ALLCOORDS, stimNum);
    T = add_polar_about_s(T, ALLCOORDS, stimNum);
    T = add_arc_about_s_edges(T, ALLCOORDS, RTAB384, stimNum);
    [T, widthPx] = add_GC_normalization(T, RTAB384, stimNum);
    
    Tall(stimNum).stimNum = stimNum;
    Tall(stimNum).widthPx = widthPx;
    Tall(stimNum).T       = T;
    
    % ---- progress reporting ----
    if mod(stimNum,25) == 0 || stimNum == 1 || stimNum == nStim
        fprintf('  Stimulus %d / %d completed\n', stimNum, nStim);
    end
    
end

fprintf('Done.\n');

end
