params0 = getParams();
ModelOptimParam_TL5_iter_1;

datastruct8 = load(['data/' 'protocol_03_27_2026_8mM_slack.mat']);
datastruct2 = load(['data/' 'protocol_03_27_2026_2mM_slack.mat']);

datatable8 = datastruct8.datatable;
datatable2 = datastruct2.datatable;

if isfield(datastruct8, 'features_data')
    features_data = datastruct8.features_data;
    % fill in missing parts
    features_data.t0_crossing = features_data.t0;
end

% just plot the data. Keep if needed at some point
figure(1);clf;

set(gcf, 'Position',[488 242 325*3 270])
tiledlayout('flow', 'TileSpacing','tight')
nexttile;
t0 = 74.3190;
plot(datatable8(:, 1)-t0, datatable8(:, 2)/2, 'k', LineWidth=2);
title('Slack-restretch protocol')
xlim([0  3.1]);
xlabel('Time (s)');
ylabel('Muscle length (-)');

nexttile;
l1 = plot(datatable2(:, 1)-t0, datatable2(:, 3), LineWidth=2);
hold on;
l2 = plot(datatable8(:, 1)-t0, datatable8(:, 3), LineWidth=2);
lgd = legend([l2 l1], '8 mM', '2 mM', Location='best');
title(lgd, 'ATP');

title('Measured muscle force')
xlim([0  3.1]);
ylim([0, inf])
xlabel('Time (s)');
ylabel('Force (kPa)');

nexttile;
l1 = plot(datatable2(:, 1)-t0, datatable2(:, 3), LineWidth=2);
hold on;
l2 = plot(datatable8(:, 1)-t0, datatable8(:, 3), LineWidth=2);
title('Zoom in')
xlim([0.6940    1.2279]);
ylim([0, inf])
xlabel('Time (s)');
ylabel('Force (kPa)');


%%
LoadData;

figure(2);clf;
dat = Data_ATP;

set(gcf, 'Position',[488 242 325*2 270]);
tiledlayout('flow', 'TileSpacing','compact')
nexttile;
plot(dat(:, 2), dat(:, 1), 'o-k', dat(:, end), dat(:, 1), 'o--k', LineWidth=2);
lgd = legend('8 mM', '2 mM');
title(lgd, 'ATP');

title('Force-velocity protocol');
xlabel('Force (kPa)');ylabel('Velocity (ML/s)');

nexttile;
plot(dat(:, 1), dat(:, 1).*dat(:, 2), 'o-k',dat(:, 1), dat(:, 1).*dat(:, end), 'o--k', LineWidth=2);
title('Low ATP reduces contraction power');
ylabel('Power (W/ML)');xlabel('Velocity (ML/s)');

%% Figure 3
plotThat = true;
if plotThat
    figure(3);clf;
    set(gcf, 'Position',[488 242 325 270]);
    xlim([t0+0.6940    t0+1.2279]);
end

features_data8 = extractSlackAttributes(datatable8(:, 1), datatable8(:, 3), datatable8(:, 2), datastruct8.velocitytable, struct(), [], plotThat);
hold on;
features_data2 = extractSlackAttributes(datatable2(:, 1), datatable2(:, 3), datatable2(:, 2), datastruct2.velocitytable, struct(), [], plotThat);
%%
figure(3);clf;
% params0 = getParams();
% ModelParams_tuesdayLunch;
% RunBakersExp;

fn = {'ktr|SLslack'        , ...
'A|SLslack'          , ...
't0|SLdiff'          , ...
'restretchSlopeStart', ...
'peak1_y|SLslack', ...
'peak1_dSL|SLslack'          , ...
'vall_y'             , ...
'peak2'              , ...
'steady'             , ...
'vall2_dy|SLslack'           , ...
'ovrsht_dy'          , ...
'ovrsht_t'           }
% plotFeatures(features_data8, features_data2, features_model, fn)
plotFeatures(features_data8, features_data2, [], fn)
delete(gca);

%% Synchronize ylims: apply fig 3 limits to fig 7
fig3 = figure(80085);
axes_list3 = findobj(fig3, 'type', 'axes');
load('feature_panels_ylims.mat')
fig_titles = {
    'Overshoot (dt)'            , ...
'Overshoot (dF)'           , ...
'Valley2 (dF)' , ...
'Isometric (F)'              , ...
'Peak2 (F)'               , ...
'Valley (dF)' , ...
'Peak1 (dSL)', ...
'Peak1 (dF)'  , ...
'Restretch slope (dF/dt)' , ...
'T0 (dt)'        , ...
'A (F)'        , ...
'ktr (1/t)'      };

for i = 1:min(length(axes_list3), length(axes_list3))
    axes_list3(i).YLim = ylims(:, i);
    axes_list3(i).Title.String = fig_titles{i};
end

nexttile(1);lgd = legend('8 mM', '2 mM', fontsize=12, location='best', orientation='horizontal');
% title(lgd, 'ATP', FontSize=12);
set(gcf, 'Position', [264.2000 90.6000 727.2000 663.2000])
%% Figure 4
% see identify passive
%% Figure 5
figure(5);clf;
params0 = getParams();
ModelParams_tuesdayLunch;
params0.justPlotStateTransitionsFlag = true;
try
RunBakersExp;
catch
end
set(gcf, "Position",  [385 154.6000 1.0128e+03 500])
% exportgraphics(gca, 'transitions.png', 'Resolution', 300);

%% Export each panel separately as high-res PNG
output_dir = './Figures/TransitionPanels/';
if ~isfolder(output_dir)
    mkdir(output_dir);
end

fig = figure(5);
tlo = fig.Children(1);

panel_names = {
    'Panel_1_Hydrolysis'
    'Panel_2_Attachment_P1_Detachment'
    'Panel_3_Power_Stroke'
    'Panel_4_P2_Detachment'
    'Panel_5_SRX_ATP'
    'Panel_6_SRX_ADP'
    'Panel_7_Inter_SRX'
    'Panel_8_Summary'
};

axes_list = findobj(tlo, 'type', 'axes');
axes_list = flip(axes_list);

fprintf('Exporting %d panels at 600 DPI...\n', length(axes_list));
for i = 1:length(axes_list)
    ax = axes_list(i);
    output_file = fullfile(output_dir, [panel_names{i} '.png']);
    exportgraphics(ax, output_file, 'Resolution', 600, 'Colorspace', 'rgb');
    fprintf('  [%d/8] %s\n', i, panel_names{i});
end
fprintf('Done! Panels saved to %s\n', output_dir);
%% 
figure(6); clf;
ModelOptimParam_TL5_iter_81
params0.justPlotStateTransitionsFlag = false;
params0.RunForceVelocity = false;
params0.RunSlack = true;

RunBakersExp;
%% manipulate the output
set(gcf, "Position",  [385 154.6000 1.0128e+03 500]);
% tl1 = nexttile(1);
% delete(tl1)
%
h = findobj(gca, 'Type', 'Line');

offset = xlim;
for k = 1:numel(h)
    h(k).XData = h(k).XData - offset(1);
end

xlim([0 3]);
delete(h(6));delete(h(3));delete(h(4));
yyaxis left;ylim([-1 100]);ylabel('Force (kPa)')
yyaxis right; ylim([1.8 2.6]);ylabel('Length (\mum)')
legend('Location', 'southwest');
xlabel('Time');
%% an inset
title('Slack segment zoom-in')
set(gcf, "Position",  [385 154.6000 1.0128e+03 500]);
xlim([0.5889    1.1169]);
%% Figure 7: features
figure(7); clf;
plotFeatures(features_data, features_model, [], fn);
delete(gca);


%% Synchronize ylims: apply fig 3 limits to fig 7
fig3 = figure(80085);
axes_list3 = findobj(fig3, 'type', 'axes');

% plotFeatures(features_data, features_model, features_data2, fn);
% delete(gca);
% fig3 = figure(80085);
% axes_list3 = findobj(fig3, 'type', 'axes');
% ylims = [];
% for i = 1:min(length(axes_list3), length(axes_list3))
%     ylims(:, i) = axes_list3(i).YLim;
% end
% save('feature_panels_ylims.mat', 'ylims')


for i = 1:min(length(axes_list3), length(axes_list3))
    axes_list3(i).YLim = ylims(:, i);
    axes_list3(i).Title.String = fig_titles{i};
end
nexttile(1);lgd = legend('Data 8 mM', 'Model', fontsize=12, location='best', orientation='vertical');
set(gcf, 'Position', [264.2000 90.6000 727.2000 663.2000])
%%
figure(8); clf;
ModelOptimParam_TL5_iter_81
params0.justPlotStateTransitionsFlag = false;
params0.FV_velocities = -[0 0.5, 1, 2, 3, 4, 5, 6];
params0.RunForceVelocity = true;
params0.RunSlack = false;

RunBakersExp;
set(gcf, 'Position', [488 322.6000 377.8000 339.4000]);
legend('Data 8 mM', 'Model')