%% 28 dec 2025, try to get a good latency estimate, first focusing on Figaro


SaveFilePdf = 'Nilson.pdf';

Monkey = 1; % 1 for Nilson, 2 for Figaro
SNR_th = 1;                 % activity threshold for good enough SNR 

p_val_local = 1;               % if true, looks at the p-value in the actual distance condition
p_th = 0.075;

p_th_arr = [0.2,0.2,0.2];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;         % if set to 1, it will use the response in the object attention task, not the global SNR

Latency_analysis_SVM

% combination 2
SNR_th = 1.4;                 % activity threshold for good enough SNR 
p_val_local = 1;               % if true, looks at the p-value in the actual distance condition
p_th = 0.075;
p_th_arr = [0.2,0.2,0.2];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;         % if set to 1, it will use the response in the object attention task, not the global SNR
SaveFilePdf = 'Nilson_1.4_0.075_Loc.pdf';
Latency_analysis_SVM

% combination 3
SNR_th = 1.0;                 % activity threshold for good enough SNR 
p_val_local = 1;               % if true, looks at the p-value in the actual distance condition
p_th = 0.05;
p_th_arr = [0.2,0.2,0.2];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;         % if set to 1, it will use the response in the object attention task, not the global SNR
SaveFilePdf = sprintf('Nilson_%.2f_%.2f_Loc.pdf', SNR_th, p_th);
Latency_analysis_SVM

% combination 4
SNR_th = 0.5;                 % activity threshold for good enough SNR 
p_val_local = 1;               % if true, looks at the p-value in the actual distance condition
p_th = 0.05;
p_th_arr = [0.2,0.2,0.2];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;         % if set to 1, it will use the response in the object attention task, not the global SNR
SaveFilePdf = sprintf('Nilson_%.2f_%.2f_Loc.pdf', SNR_th, p_th);
Latency_analysis_SVM

% combination 5
SNR_th = 0.85;                 % activity threshold for good enough SNR 
p_val_local = 1;               % if true, looks at the p-value in the actual distance condition
p_th = 0.05;
p_th_arr = [0.2,0.2,0.2];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;         % if set to 1, it will use the response in the object attention task, not the global SNR
SaveFilePdf = sprintf('Nilson_%.2f_%.2f_Loc.pdf', SNR_th, p_th);
Latency_analysis_SVM

% combination 6  Check p-values in all cone conditions, i.e. global
SNR_th = 0.85;                 % activity threshold for good enough SNR 
p_val_local = 0;               % if true, looks at the p-value in the actual distance condition
p_th = 0.05;
p_th_arr = [0.15,0.15,0.15];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;             % if set to 1, it will use the response in the object attention task, not the global SNR
pStr = sprintf('%.2f-', p_th_arr);     % '0.20-0.20-0.20-'
pStr(end) = [];                        % remove trailing '-'
SaveFilePdf = sprintf('Nilson_%.2f_%s_Glob.pdf', SNR_th, pStr);
Latency_analysis_SVM

% combination 7  Check p-values in all cone conditions, i.e. global
SNR_th = 0.85;                 % activity threshold for good enough SNR 
p_val_local = 0;               % if true, looks at the p-value in the actual distance condition
p_th = 0.05;
p_th_arr = [0.2,0.2,0.2];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;             % if set to 1, it will use the response in the object attention task, not the global SNR
pStr = sprintf('%.2f-', p_th_arr);     % '0.20-0.20-0.20-'
pStr(end) = [];                        % remove trailing '-'
SaveFilePdf = sprintf('Nilson_%.2f_%s_Glob.pdf', SNR_th, pStr);
Latency_analysis_SVM

% combination 8  Check p-values in all cone conditions, i.e. global
SNR_th = 0.7;                 % activity threshold for good enough SNR 
p_val_local = 0;               % if true, looks at the p-value in the actual distance condition
p_th = 0.05;
p_th_arr = [0.5,0.5,0.5];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;             % if set to 1, it will use the response in the object attention task, not the global SNR
pStr = sprintf('%.2f-', p_th_arr);     % '0.20-0.20-0.20-'
pStr(end) = [];                        % remove trailing '-'
SaveFilePdf = sprintf('Nilson_%.2f_%s_Glob.pdf', SNR_th, pStr);
Latency_analysis_SVM

% combination 9  Check p-values in all cone conditions, i.e. global
SNR_th = 0.8;                 % activity threshold for good enough SNR 
p_val_local = 0;               % if true, looks at the p-value in the actual distance condition
p_th = 0.05;
p_th_arr = [1,1,1];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;             % if set to 1, it will use the response in the object attention task, not the global SNR
pStr = sprintf('%.2f-', p_th_arr);     % '0.20-0.20-0.20-'
pStr(end) = [];                        % remove trailing '-'
SaveFilePdf = sprintf('Nilson_%.2f_%s_Glob.pdf', SNR_th, pStr);
Latency_analysis_SVM

% combination 10  Check p-values in all cone conditions, i.e. global
SNR_th = 1;                 % activity threshold for good enough SNR 
p_val_local = 0;               % if true, looks at the p-value in the actual distance condition
p_th = 0.05;
p_th_arr = [0.5,0.5,0.5];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;             % if set to 1, it will use the response in the object attention task, not the global SNR
pStr = sprintf('%.2f-', p_th_arr);     % '0.20-0.20-0.20-'
pStr(end) = [];                        % remove trailing '-'
SaveFilePdf = sprintf('Nilson_%.2f_%s_Glob.pdf', SNR_th, pStr);
Latency_analysis_SVM

% combination 11  Check p-values in all cone conditions, i.e. global
SNR_th = 1;                 % activity threshold for good enough SNR 
p_val_local = 0;               % if true, looks at the p-value in the actual distance condition
p_th = 0.05;
p_th_arr = [0.2,0.2,0.2];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;             % if set to 1, it will use the response in the object attention task, not the global SNR
pStr = sprintf('%.2f-', p_th_arr);     % '0.20-0.20-0.20-'
pStr(end) = [];                        % remove trailing '-'
SaveFilePdf = sprintf('Nilson_%.2f_%s_Glob.pdf', SNR_th, pStr);
Latency_analysis_SVM

% combination 12  Check p-values in all cone conditions, i.e. global
SNR_th = 0.6;                 % activity threshold for good enough SNR 
p_val_local = 0;               % if true, looks at the p-value in the actual distance condition
p_th = 0.05;
p_th_arr = [0.2,0.2,0.2];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;             % if set to 1, it will use the response in the object attention task, not the global SNR
pStr = sprintf('%.2f-', p_th_arr);     % '0.20-0.20-0.20-'
pStr(end) = [];                        % remove trailing '-'
SaveFilePdf = sprintf('Nilson_%.2f_%s_Glob.pdf', SNR_th, pStr);
Latency_analysis_SVM

% combination 13  Check p-values in all cone conditions, i.e. global
SNR_th = 0.65;                 % activity threshold for good enough SNR 
p_val_local = 0;               % if true, looks at the p-value in the actual distance condition
p_th = 0.05;
p_th_arr = [0.2,0.2,0.2];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;             % if set to 1, it will use the response in the object attention task, not the global SNR
pStr = sprintf('%.2f-', p_th_arr);     % '0.20-0.20-0.20-'
pStr(end) = [];                        % remove trailing '-'
SaveFilePdf = sprintf('Nilson_%.2f_%s_Glob.pdf', SNR_th, pStr);
Latency_analysis_SVM

% combination 14  Check p-values in all cone conditions, i.e. global
SNR_th = 0.65;                 % activity threshold for good enough SNR 
p_val_local = 0;               % if true, looks at the p-value in the actual distance condition
p_th = 0.05;
p_th_arr = [0.5,0.5,0.5];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;             % if set to 1, it will use the response in the object attention task, not the global SNR
pStr = sprintf('%.2f-', p_th_arr);     % '0.20-0.20-0.20-'
pStr(end) = [];                        % remove trailing '-'
SaveFilePdf = sprintf('Nilson_%.2f_%s_Glob.pdf', SNR_th, pStr);
Latency_analysis_SVM

% combination 15  Check p-values in all cone conditions, i.e. global
SNR_th = 0.85;                 % activity threshold for good enough SNR 
p_val_local = 0;               % if true, looks at the p-value in the actual distance condition
p_th = 0.05;
p_th_arr = [1,1,1];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;             % if set to 1, it will use the response in the object attention task, not the global SNR
pStr = sprintf('%.2f-', p_th_arr);     % '0.20-0.20-0.20-'
pStr(end) = [];                        % remove trailing '-'
SaveFilePdf = sprintf('Nilson_%.2f_%s_Glob.pdf', SNR_th, pStr);
Latency_analysis_SVM


% now Figaro

Monkey = 2; % 1 for Nilson, 2 for Figaro
SNR_th = 0.85;                 % activity threshold for good enough SNR 
p_val_local = 0;               % if true, looks at the p-value in the actual distance condition
p_th = 0.05;
p_th_arr = [0.15,0.15,0.15];   % if p_val_local is false, only includes channesl that reach these significances in the three distance conditions
MinSigma = 15;                  % maximal signal of fitted curve
MaxSigma = 75;                  % maximal signal of fitted curve
SelectLocalSNR = 1;         % if set to 1, it will use the response in the object attention task, not the global SNR
pStr = sprintf('%.2f-', p_th_arr);     % '0.20-0.20-0.20-'
pStr(end) = [];                        % remove trailing '-'
SaveFilePdf = sprintf('bbFigaro_%.2f_%s_Glob.pdf', SNR_th, pStr);
Latency_analysis_SVM




