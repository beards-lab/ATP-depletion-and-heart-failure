% MakeSlide.m — assemble the four slide figures into a single 16:9 .pptx slide.
%
% The figures are 3:2, so a 2x2 block is 3:2 too. Fitted to the INNER content
% container of the slide (below the title, inside the margins) at its natural
% aspect ratio, which leaves a column free on the right — used for the three
% numbers a reader needs to interpret the panels.
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

%% ---- geometry (inches) --------------------------------------------------
SW = 13.3333; SH = 7.5;                 % 16:9
MG = 0.5;                               % slide margin
TITLE_H = 0.52; SUB_H = 0.30;
Y0 = MG + TITLE_H + SUB_H + 0.16;       % top of the figure block
GAP = 0.14;
BH = SH - Y0 - MG;                      % block height the container allows
FH = (BH - GAP)/2;  FW = FH*1.5;        % keep each figure 3:2
BW = 2*FW + GAP;
X0 = MG;

fprintf('figure %.2f x %.2f in | block %.2f x %.2f in | free column %.2f in\n', ...
        FW, FH, BW, BH, SW - X0 - BW - MG);

s = struct('width_in',SW,'height_in',SH,'bg','FFFFFF');

%% ---- text ---------------------------------------------------------------
ACC = 'B01118'; FG = '1F2933'; MUT = '6B7280';
T = struct('str',{},'x',{},'y',{},'w',{},'h',{},'size',{},'bold',{},'color',{},'align',{});
T(end+1) = struct('str','Cross-bridge model — spread-weighted refit', ...
    'x',MG,'y',MG-0.12,'w',SW-2*MG,'h',TITLE_H,'size',30,'bold',true,'color',FG,'align','l');
T(end+1) = struct('str',sprintf(['8 mM, preparation 03/27 · objective %.4f → %.4f · ' ...
    'staircase and k_tr were never fitted'], c1, c2), ...
    'x',MG,'y',MG+TITLE_H-0.10,'w',SW-2*MG,'h',SUB_H,'size',14,'bold',false,'color',MUT,'align','l');

XS = X0 + BW + 0.34; WS = SW - XS - MG;      % the free right column
ys = Y0 + 0.10;
stats = { sprintf('%.2f', c2),   sprintf('objective, from %.2f', c1)
          sprintf('%.0f %%', fvShare), 'of the gain is force–velocity'
          '+37 %',              'staircase / k_tr overshoot — the open defect' };
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
F = { 'slide_fv.png',            X0,        Y0
      'slide_stairs.png',        X0+FW+GAP, Y0
      'slide_restretch_all.png', X0,        Y0+FH+GAP
      'slide_restretch_det.png', X0+FW+GAP, Y0+FH+GAP };
I = struct('file',{},'x',{},'y',{},'w',{},'h',{});
for i = 1:size(F,1)
    f = fullfile(resDir, F{i,1});
    if ~isfile(f)
        error('MakeSlide:missingFigure','run MakeSlideFigures.m first — missing %s', F{i,1});
    end
    I(end+1) = struct('file',f,'x',F{i,2},'y',F{i,3},'w',FW,'h',FH);     %#ok<SAGROW>
end
s.images = I;

out = buildPptxSlide(fullfile(resDir, sprintf('slide_%s.pptx', TAG_NEW)), s);
fprintf('open: %s\n', out);
