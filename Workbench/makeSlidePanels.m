function makeSlidePanels(cacheFile, outdir, tag, ttl)
%MAKESLIDEPANELS  Slide-legible panels from a cached slack run.
%
%   makeSlidePanels(cacheFile, outdir, tag, ttl)
%
%   cacheFile must contain out_slack, fm (model features) and fd (data
%   features), as written by Export8mMPanels / Export2mMPanels.
%
%   Writes  <tag>_features.png   curated scored-feature grid, data o vs model x
%           <tag>_zoom1.png      first slack segment, force + SL, features marked
%
%   Unlike Auxiliary/plotFeatures (which tiles EVERY selector and becomes
%   unreadable when exported), this plots a hand-picked set at a size that
%   survives a projector.

if nargin < 4, ttl = tag; end
S = load(cacheFile);
out = S.out_slack; fm = S.fm; fd = S.fd;
if ~exist(outdir,'dir'); mkdir(outdir); end

% ---- feature grid -------------------------------------------------------
% name, label, unit-class  (used only for the y-label)
F = { 'A',                  'A  (amplitude, kPa)'
      'steady',             'steady force (kPa)'
      'peak1_y',            'peak1  (kPa)'
      'vall_y',             'valley  (kPa)'
      'peak2',              'peak2  (kPa)'
      'vall2_dy',           'vall2\_dy  (kPa)'
      'ktr',                'k_{tr}  (s^{-1})'
      'rsK',                'rsK — post-restretch rate (s^{-1})'
      'restretchSlopeStart','restretch slope (kPa/s)'
      'peak1_dSL',          'peak1\_dSL  (\mum)'
      'vall_t',             'valley time (s)'
      'ovrsht_dy',          'ovrsht\_dy  (kPa)' };

fig = figure('Units','pixels','Position',[20 20 1900 1150],'Color','w');
tl = tiledlayout(3,4,'Padding','compact','TileSpacing','compact');
for i = 1:size(F,1)
    nexttile; hold on; box on; grid on
    nm = F{i,1};
    d = getfielddef(fd, nm); m = getfielddef(fm, nm);
    x = 1:max([numel(d) numel(m) 1]);
    if ~isempty(d), plot(1:numel(d), d, 'ko-','MarkerSize',9,'LineWidth',1.6,'MarkerFaceColor','w'); end
    if ~isempty(m), plot(1:numel(m), m, 'rx--','MarkerSize',11,'LineWidth',1.8); end
    xlim([0.6 max(x)+0.4]); xticks(1:max(x));
    ttlstr = F{i,2};
    if ~isempty(d) && ~isempty(m) && numel(d)==numel(m)
        rel = mean(abs(m(:)-d(:))./max(abs(d(:)),eps),'omitnan');
        ttlstr = sprintf('%s   (%.0f%%)', ttlstr, 100*rel);
    end
    title(ttlstr,'FontSize',12,'FontWeight','bold','Interpreter','tex');
    if i > 8, xlabel('slack segment','FontSize',11); end
    set(gca,'FontSize',11);
end
lg = legend({'data','model'},'Orientation','horizontal','FontSize',13);
lg.Layout.Tile = 'north';
title(tl, sprintf('%s — scored slack features, data (o) vs model (x); %% = mean relative error', ttl), ...
      'FontSize',15,'FontWeight','bold');
exportgraphics(fig, fullfile(outdir,[tag '_features.png']), 'Resolution', 170);
fprintf('  saved %s_features.png\n', tag);

% ---- zoom on the first slack segment -----------------------------------
t = out.t(:); Fm = out.Force(:);
dt = out.datatable;                    % [t, SL/ML, F] as loaded by the experiment
td = dt(:,1); Fd = dt(:,end);
t0 = fd.t_seg(1);
w  = [t0-0.06, t0+0.62];

fig2 = figure('Units','pixels','Position',[20 20 1700 1000],'Color','w');
tl2 = tiledlayout(4,1,'Padding','compact','TileSpacing','none');

ax1 = nexttile([3 1]); hold on; box on; grid on
plot(td, Fd, 'k-', 'LineWidth',1.0, 'DisplayName','data');
plot(t,  Fm, 'r-', 'LineWidth',2.0, 'DisplayName','model');
xlim(w); set(ax1,'XTickLabel',[]); ylabel('force (kPa)','FontSize',13);
set(ax1,'FontSize',12);

% feature markers, read from the DATA features of segment 1
mk = @(tt,yy,txt,col) [plot(tt,yy,'o','MarkerSize',10,'LineWidth',2,'Color',col, ...
                            'MarkerFaceColor','w','HandleVisibility','off'), ...
                       text(tt,yy,['  ' txt],'Color',col,'FontSize',12,'FontWeight','bold', ...
                            'VerticalAlignment','bottom')];
cA = [0.85 0.33 0.10]; cB = [0.00 0.45 0.74]; cC = [0.47 0.67 0.19];
try
    ys = fd.steady(1);
    yline(ys,'--','steady','Color',cC,'LineWidth',1.5,'FontSize',12, ...
          'LabelHorizontalAlignment','left','HandleVisibility','off');
    % Locate peak1 / valley / peak2 ON THE DATA TRACE rather than from feature
    % arithmetic: the re-stretch is the largest positive dF/dt in the segment.
    inw = td >= t0 & td <= t0 + 0.62;
    tw  = td(inw); Fw = Fd(inw); Fs = movmean(Fw,9);   % smooth: the trace is noisy
    [~,iR] = max(diff(Fs));                       % re-stretch instant
    seg = (1:numel(tw))' > iR & tw < tw(iR) + 0.09;
    [y1,i1] = max(Fs(seg));  idx = find(seg);  i1 = idx(i1);
    aft = (1:numel(tw))' > i1 & tw < tw(i1) + 0.05;
    [yv,iv] = min(Fs(aft));  idx = find(aft);  iv = idx(iv);
    aft2 = (1:numel(tw))' > iv & tw < tw(iv) + 0.030;
    [y2,i2] = max(Fs(aft2)); idx = find(aft2); i2 = idx(i2);
    mk(tw(i1), y1, 'peak1', cA);
    mk(tw(iv), yv, 'valley', cA);
    h = mk(tw(i2), y2, 'peak2', cA); set(h(2),'VerticalAlignment','top');
    text(t0+0.055, fd.A(1)*0.42, sprintf('k_{tr} = %.0f s^{-1}\nA = %.0f kPa', fd.ktr(1), fd.A(1)), ...
         'Color',cB,'FontSize',13,'FontWeight','bold');
    text(tw(i2)+0.055, fd.steady(1)*0.50, ...
         sprintf('post-restretch recovery\nrsK: data %.0f  |  model %.0f s^{-1}', fd.rsK(1), fm.rsK(1)), ...
         'Color',cB,'FontSize',13,'FontWeight','bold');
catch ME
    fprintf('  [markers skipped] %s\n', ME.message);
end
legend('Location','southeast','FontSize',13);
title(sprintf('%s — first slack segment (the window every scored feature is read from)', ttl), ...
      'FontSize',15,'FontWeight','bold');

ax2 = nexttile; hold on; box on; grid on
plot(t, out.SL(:), '-', 'Color',[0.6 0.1 0.1],'LineWidth',2);
xlim(w); ylabel('SL (\mum)','FontSize',13); xlabel('time (s)','FontSize',13);
slw = out.SL(t>=w(1) & t<=w(2));
if ~isempty(slw)
    pad = 0.12*max(range(slw), 1e-3);
    ylim([min(slw)-pad, max(slw)+pad]);
end
set(ax2,'FontSize',12);
linkaxes([ax1 ax2],'x');
exportgraphics(fig2, fullfile(outdir,[tag '_zoom1.png']), 'Resolution', 170);
fprintf('  saved %s_zoom1.png\n', tag);
end

function v = getfielddef(s, nm)
if isfield(s,nm); v = double(s.(nm)(:)); v = v(:)'; else; v = []; end
end
