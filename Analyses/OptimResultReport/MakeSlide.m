% MakeSlide.m — assemble the four slide figures into a single 16:9 .pptx slide.
%
% Three of the four panels are two cells wide (3.04:1) and only F-V is 3:2, so
% they are laid out by their OWN aspect ratios rather than on a uniform grid:
% the wide time-series panels stack in a left column, F-V and the headline
% numbers take the right. Everything sits inside the slide's content container
% (below the title, inside the margins).
%
% No Report Generator toolbox, no Office automation, no Python: Auxiliary/
% buildPptxSlide.m writes the OOXML parts and zips them.
%
%   run('Analyses/OptimResultReport/MakeSlideFigures.m')   % the four PNGs
%   run('Analyses/OptimResultReport/MakeSlide.m')          % the .pptx

clearvars -except TAG_NEW; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resDir = fullfile(here,'results');
if ~exist('TAG_NEW','var'); TAG_NEW = 'spreadw_opt'; end

L = load(fullfile(resDir, sprintf('fit_report_%s.mat', TAG_NEW)));
c1 = sum(L.c1); c2 = sum(L.c2);
iFV = find(strcmp(L.nm,'FV_normpowerAvg'),1);
fvShare = 100*(L.c1(iFV)-L.c2(iFV))/(c1-c2);

% Staircase step-peak overshoot — computed, not typed, so the headline number
% on the slide cannot go stale behind a refit.
Sst = load(fullfile(root,'data','protocol_03_27_2026_velocitytable_stairs.mat'));
ds = Sst.datatable; vs = Sst.velocitytable;
tst = vs(vs(:,2)~=0,1); tst = tst(diff([0;tst]) > 0.02);
t0s = vs(find(vs(:,2)~=0,1),1);
o = L.R(2).out_stairs;
F0 = interp1(o.t,o.Force,t0s,'linear','extrap'); if F0 == 0; F0 = 1; end
Fn = o.Force/F0; e = nan(1,numel(tst));
for j = 1:numel(tst)
    w = ds(:,1) > tst(j) & ds(:,1) < tst(j)+0.06;  if ~any(w); continue; end
    mm = o.t > tst(j) & o.t < tst(j)+0.06;
    e(j) = max(Fn(mm))/max(movmean(ds(w,3),9)) - 1;
end
stairsOver = 100*mean(e,'omitnan');

%% ---- geometry (inches) --------------------------------------------------
% Three panels are now two cells wide (3.04:1) and only F-V is 3:2, so a
% uniform grid no longer tiles. Layout instead: the three wide time-series
% panels stacked in a left column, F-V and the numbers in a right column. Every
% height is derived from the PNG's OWN pixel aspect, so nothing is distorted
% and the slide follows whatever MakeSlideFigures.m produced.
SW = 13.3333; SH = 7.5;                 % 16:9
MG = 0.5;                               % slide margin
TITLE_H = 0.52; SUB_H = 0.30;
Y0 = MG + TITLE_H + SUB_H + 0.16;       % top of the figure block
GAP = 0.14;
BH = SH - Y0 - MG;                      % block height the container allows

WIDE = {'slide_restretch_all.png','slide_restretch_det.png','slide_stairs.png'};
asp = @(f) local_aspect(fullfile(resDir,f));           % width / height
aW  = mean(cellfun(asp, WIDE));                        % the wide panels' aspect
% three stacked wide panels must fit BH:  3*(W/aW) + 2*GAP = BH
LW  = min((BH - 2*GAP) * aW / 3, 0.58*(SW - 2*MG));
LH  = LW / aW;
X0  = MG;  XR = X0 + LW + 0.34;  WR = SW - XR - MG;    % right column
FVW = min(WR, 4.3);  FVH = FVW / asp('slide_fv.png');  % F-V, top of that column
fprintf('wide panels %.2f x %.2f in (x3 stacked) | F-V %.2f x %.2f | right column %.2f in\n', ...
        LW, LH, FVW, FVH, WR);

s = struct('width_in',SW,'height_in',SH,'bg','FFFFFF');

%% ---- text ---------------------------------------------------------------
ACC = 'B01118'; FG = '1F2933'; MUT = '6B7280';
T = struct('str',{},'x',{},'y',{},'w',{},'h',{},'size',{},'bold',{},'color',{},'align',{});
T(end+1) = struct('str','Cross-bridge model — spread-weighted refit', ...
    'x',MG,'y',MG-0.12,'w',SW-2*MG,'h',TITLE_H,'size',30,'bold',true,'color',FG,'align','l');
T(end+1) = struct('str',sprintf(['8 mM, preparation 03/27 · objective %.4f → %.4f · ' ...
    'staircase and k_tr were never fitted'], c1, c2), ...
    'x',MG,'y',MG+TITLE_H-0.10,'w',SW-2*MG,'h',SUB_H,'size',14,'bold',false,'color',MUT,'align','l');

XS = XR; WS = WR;                            % the right column
ys = Y0 + FVH + 0.30;                        % below the F-V panel
stats = { sprintf('%.2f', c2),   sprintf('objective, from %.2f', c1)
          sprintf('%.0f %%', fvShare), 'of the gain is force–velocity'
          sprintf('%+.1f %%', stairsOver), 'staircase step-peak overshoot — never fitted' };
for i = 1:size(stats,1)
    T(end+1) = struct('str',stats{i,1},'x',XS,'y',ys,'w',WS,'h',0.62, ...
        'size',34,'bold',true,'color',ACC,'align','l');                 %#ok<SAGROW>
    T(end+1) = struct('str',stats{i,2},'x',XS,'y',ys+0.56,'w',WS,'h',0.62, ...
        'size',13,'bold',false,'color',MUT,'align','l');                %#ok<SAGROW>
    ys = ys + 1.42;
end
T(end+1) = struct('str',sprintf('%s · not converged', TAG_NEW), ...
    'x',XS,'y',Y0+BH-0.30,'w',WS,'h',0.30,'size',11,'bold',false,'color',MUT,'align','l');
s.texts = T;

%% ---- figures ------------------------------------------------------------
% file, x, y, width — each height comes from the file's own aspect ratio
F = { WIDE{1}, X0, Y0,                  LW      % slack protocol
      WIDE{2}, X0, Y0 + (LH+GAP),       LW      % slack cycle 4
      WIDE{3}, X0, Y0 + 2*(LH+GAP),     LW      % staircase
      'slide_fv.png', XR, Y0,           FVW };  % force-velocity, right column
I = struct('file',{},'x',{},'y',{},'w',{},'h',{});
for i = 1:size(F,1)
    f = fullfile(resDir, F{i,1});
    if ~isfile(f)
        error('MakeSlide:missingFigure','run MakeSlideFigures.m first — missing %s', F{i,1});
    end
    I(end+1) = struct('file',f,'x',F{i,2},'y',F{i,3},'w',F{i,4}, ...
                      'h',F{i,4}/local_aspect(f));                       %#ok<SAGROW>
end
s.images = I;

out = buildPptxSlide(fullfile(resDir, sprintf('slide_%s.pptx', TAG_NEW)), s);
fprintf('open: %s\n', out);

function a = local_aspect(f)
% width / height of a PNG, so placements never distort a figure.
    if ~isfile(f); f = fullfile('Analyses','OptimResultReport','results', f); end
    info = imfinfo(f);
    a = info(1).Width / info(1).Height;
end
