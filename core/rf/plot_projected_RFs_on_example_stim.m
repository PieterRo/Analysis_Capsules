function h = plot_projected_RFs_on_example_stim(Tall, ALLCOORDS, RTAB384, exampleStimNum, varargin)
% PLOT_PROJECTED_RFS_ON_EXAMPLE_STIM
%
% Projects RFs from all stimuli into the coordinate frame of one example stimulus.
% Uses:
%   - along_GC / perp_signed_GC for object RFs
%   - r_s_GC / arc_frac_edge / arc_isInner_edge for background RFs
%
% Produces:
%   - Expanded axes (no clipping)
%   - Grey background everywhere
%   - Thin light-grey frame marking original 1024x768 screen

%% -------------------- Options --------------------
p = inputParser;
p.addParameter('MarkerSize', 4);
p.addParameter('Alpha', 0.15);
p.addParameter('SiteIdx', [], @(x) isempty(x) || (isnumeric(x) && isvector(x)));
p.parse(varargin{:});
opt = p.Results;

W = 1024;
H = 768;

%% -------------------- Geometry of Example Stimulus --------------------
toPx = @(q) [q(1) + W/2, H/2 - q(2)];

fieldName = sprintf('stim_%d', exampleStimNum);
s     = double(ALLCOORDS.(fieldName).s(:))';
tFig  = double(ALLCOORDS.(fieldName).t_fig(:))';
tBack = double(ALLCOORDS.(fieldName).t_back(:))';

s_px = toPx(s);
tT   = toPx(tFig);
tD   = toPx(tBack);

uT = (tT - s_px) / norm(tT - s_px);
uD = (tD - s_px) / norm(tD - s_px);

nT = perpTowardOther(s_px, tT, tD);
nD = perpTowardOther(s_px, tD, tT);

widthEx = double(RTAB384(exampleStimNum,7));
radEx   = widthEx/2;

alphaEx = signedAngle(uD, uT);
sgn = sign(alphaEx); if sgn==0, sgn=1; end
alphaLongEx = alphaEx - sgn*2*pi;

%% -------------------- Collect All Points --------------------
xT=[]; yT=[]; xD=[]; yD=[]; xB=[]; yB=[];

for k = 1:numel(Tall)
    T = Tall(k).T;

    % Optional: restrict to a subset of sites (rows)
    if ~isempty(opt.SiteIdx)
        T = T(opt.SiteIdx, :);   % works if T is a table (it is in your case)
    end
    assign = string(T.assignment);

    % ---- Object RFs ----
    idxT = assign=="target" & ~isnan(T.along_GC);
    idxD = assign=="distractor" & ~isnan(T.along_GC);

    if any(idxT)
        along = T.along_GC(idxT) * widthEx;
        perp  = T.perp_signed_GC(idxT) * widthEx;
        pT = s_px + along.*uT + perp.*nT;
        xT=[xT;pT(:,1)]; yT=[yT;pT(:,2)];
    end

    if any(idxD)
        along = T.along_GC(idxD) * widthEx;
        perp  = T.perp_signed_GC(idxD) * widthEx;
        pD = s_px + along.*uD + perp.*nD;
        xD=[xD;pD(:,1)]; yD=[yD;pD(:,2)];
    end

    % ---- Background RFs (EDGE-BASED ARC) ----
    idxB = assign=="background" & ~isnan(T.r_s_GC);
    if any(idxB)
        r_px = T.r_s_GC(idxB) * widthEx;
        frac = T.arc_frac_edge(idxB);
        isIn = T.arc_isInner_edge(idxB);

        r_eff = max(r_px, radEx + 1e-6);
        Delta = asin(min(1, radEx ./ r_eff));

        alphaFree     = alphaEx     - sgn*(Delta+Delta);
        alphaLongFree = alphaLongEx + sgn*(Delta+Delta);

        beta = zeros(size(frac));

        % INNER ARC (start at inner edge)
        beta(isIn)  =  sgn*Delta(isIn) + frac(isIn).*alphaFree(isIn);

        % OUTER ARC (start at outer edge)
        beta(~isIn) = -sgn*Delta(~isIn) + frac(~isIn).*alphaLongFree(~isIn);

        c = cos(beta); s = sin(beta);
        vx = c*uD(1) - s*uD(2);
        vy = s*uD(1) + c*uD(2);

        pB = s_px + [r_px.*vx, r_px.*vy];
        xB=[xB;pB(:,1)]; yB=[yB;pB(:,2)];
    end
end

%% -------------------- Plotting --------------------
figure('Color',[0.5 0.5 0.5]);  % grey figure background
ax = axes('Position',[0 0 1 1]);
hold on

% Show stimulus image
img = render_stim_from_ALLCOORDS(ALLCOORDS, RTAB384, exampleStimNum);
imshow(img,'Parent',ax,'InitialMagnification','fit');

% Force axes to fill figure AFTER imshow (important!)
set(ax,'Position',[0 0 1 1]);
set(ax,'Color',[0.5 0.5 0.5]);  % grey outside image
axis ij

% Scatter points
scatter(xT,yT,opt.MarkerSize,[0.5 0.1 0.1],'filled','MarkerFaceAlpha',opt.Alpha);
scatter(xD,yD,opt.MarkerSize,[0.1 0.1 0.8],'filled','MarkerFaceAlpha',opt.Alpha);
scatter(xB,yB,opt.MarkerSize,[0.8 0.8 0.8],'filled','MarkerFaceAlpha',opt.Alpha);

% Expand axes
xAll = [xT; xD; xB; 1; W];
yAll = [yT; yD; yB; 1; H];

margin = 20;
xlim([min(xAll)-margin, max(xAll)+margin]);
ylim([min(yAll)-margin, max(yAll)+margin]);
axis equal
set(gca,'YDir','reverse')

% Thin frame around original screen
hFrame = rectangle('Position',[0.5 0.5 W H], ...
                   'EdgeColor',[0.85 0.85 0.85], ...
                   'LineWidth',1);
uistack(hFrame,'top');
% Thin frame around original screen
hFrame = rectangle('Position',[0.5 0.5 W H], ...
                   'EdgeColor',[0.85 0.85 0.85], ...
                   'LineWidth',1);
uistack(hFrame,'top');
h = struct();
h.fig = gcf;
h.ax  = ax;
h.nTarget = numel(xT);
h.nDistr  = numel(xD);
h.nBack   = numel(xB);
h.nTotal  = h.nTarget + h.nDistr + h.nBack;
% ------------------

hold off
nTarget = numel(xT);
nDistr  = numel(xD);
nBack   = numel(xB);

nTotal  = nTarget + nDistr + nBack;
fprintf('Plotted RFs: Target = %d | Distractor = %d | Background = %d | Total = %d\n', ...
    nTarget, nDistr, nBack, nTotal);

end

%% -------------------- Helpers --------------------
function n = perpTowardOther(s_px, t_arm, t_other)
v = t_arm - s_px;
u = v / norm(v);
w = t_other - s_px;
w_perp = w - dot(w,u)*u;
n = w_perp / norm(w_perp);
end

function ang = signedAngle(a,b)
ang = atan2(a(1)*b(2)-a(2)*b(1), dot(a,b));
end
