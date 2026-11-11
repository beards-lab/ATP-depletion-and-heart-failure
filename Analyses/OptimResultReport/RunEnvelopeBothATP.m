% RunEnvelopeBothATP.m -- the three-preparation feature envelope at BOTH [ATP],
% built through one code path so the 8 mM and 2 mM numbers are defined
% identically and can be drawn on the same axes.
%
% This extends RunDataEnvelope.m (8 mM only) in one way that matters:
%
% NORMALISATION. RunDataEnvelope divides every force-dimension feature by the
% mean `steady` of ITS OWN trace, which makes shape comparable across
% preparations but silently removes the ATP effect -- 2 mM force divided by
% 2 mM steady is 1 by construction, and the ~1.18x force elevation at low ATP,
% the very thing these figures exist to show, disappears. Here BOTH conditions
% of a preparation are divided by THAT PREPARATION'S mean 8 mM `steady`. One
% rule, both conditions: cross-preparation scale (77 / 65 / 48 kPa, which is
% cross-section and damage) is still divided out, while the within-preparation
% ATP effect survives.
%
% THE COST OF THAT CHOICE, stated because it is real: the 8 mM and 2 mM blocks
% were recorded at different times in a session, so the 2 mM/8 mM ratio carries
% rundown as well as the ATP effect, and the three days do not share a rundown
% (04/10 deliberately ran 2 mM FIRST -- see DataCuration). That contamination
% lands in the 2 mM band width, not in the 8 mM one.
%
% Output -> results/bothatp/envelope_bothatp.mat
%   E.mM8 / E.mM2, each with .feat.<name>.{cls,raw,mean,sd,se}, .prep, .scale8

clearvars -except IREF; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
outDir = fullfile(here,'results','bothatp'); if ~isfolder(outDir); mkdir(outDir); end

PREP = {'03/27 M','04/03 F','04/10 M'};
F8 = { 'protocol_03_27_2026_8mM_slack.mat'
       'protocol_04_03_2026_8mM_slack.mat'
       'protocol_04_10_2026_8mM_slack.mat' };
F2 = { 'protocol_03_27_2026_2mM_slack.mat'
       'protocol_04_03_2026_2mM_slack.mat'
       'protocol_04_10_2026_2mM_slack.mat' };
IREF = 1;                       % 03/27 -- the preparation the models are fitted to
nP = numel(PREP);

% class F = force-dimension (scaled by the prep's 8 mM steady); everything else
% is already scale-free. `steady` is class F here, NOT an untouched anchor as in
% RunDataEnvelope -- that is the whole point of this script.
FEAT = { 'steady','F'; 'A','F'; 'Am','F'; 'peak1_y','F'; 'vall_y','F'; 'peak2','F'
         'vall2_dy','F'; 'rsA','F'; 'ktr','R'; 'rsK','R'; 'restretchSlopeStart','R'
         't0','T'; 'peak1_t','T'; 'peak1_dSL','T'; 'vall_t','T'; 'vall2_t','T'
         'rsT0','T'; 'rsR2','Q'; 'SLslack','P'; 'v_restretch','P' };
% class P = protocol descriptor. SLslack is read off the actual SL trace and
% v_restretch off the commanded velocity table, so between-preparation spread is
% near zero by design -- these two document the ladder rather than test the model.

%% ---- load and extract, both conditions ---------------------------------
P = struct();
for i = 1:nP
    for c = 1:2
        f = F8{i}; if c == 2; f = F2{i}; end
        S = load(['data/' f]);
        d = struct('t',S.datatable(:,1),'ML',S.datatable(:,2),'F',S.datatable(:,3), ...
                   'vt',S.velocitytable, ...
                   'rel',S.velocitytable(S.velocitytable(:,2) < -1, 1));
        d.ft = extractSlackAttributes(d.t, d.F, d.ML, d.vt, [], [], false, true);
        d.sc = mean(double(d.ft.steady),'omitnan');
        if c == 1; P(i).d8 = d; else; P(i).d2 = d; end
    end
    P(i).scale8 = P(i).d8.sc;                       % the one divisor for this prep
    fprintf('%-9s  steady: 8 mM %5.1f | 2 mM %5.1f kPa  (2/8 = %.3f)\n', ...
            PREP{i}, P(i).d8.sc, P(i).d2.sc, P(i).d2.sc/P(i).d8.sc);
end
nC = min([arrayfun(@(x) numel(x.d8.rel), P), arrayfun(@(x) numel(x.d2.rel), P)]);
fprintf('common cycles: %d\n', nC);

%% ---- assemble ----------------------------------------------------------
E = struct('prep',{PREP},'iRef',IREF,'nC',nC,'scale8',[P.scale8]);
for c = 1:2
    fld = 'mM8'; if c == 2; fld = 'mM2'; end
    Q = struct();
    for j = 1:size(FEAT,1)
        nm = FEAT{j,1}; V = nan(nP,nC); ok = true;
        for i = 1:nP
            dd = P(i).d8; if c == 2; dd = P(i).d2; end
            if ~isfield(dd.ft,nm); ok = false; break; end
            v = double(dd.ft.(nm)(:))';
            if numel(v) < nC; ok = false; break; end
            v = v(1:nC);
            if FEAT{j,2} == 'F'; v = v / P(i).scale8; end
            V(i,:) = v;
        end
        if ~ok; fprintf('  skip %s (%s)\n', nm, fld); continue; end
        n = sum(~isnan(V),1);
        Q.(nm) = struct('cls',FEAT{j,2},'raw',V, ...
                        'mean',mean(V,1,'omitnan'), 'sd',std(V,0,1,'omitnan'), ...
                        'se',std(V,0,1,'omitnan')./sqrt(max(n,1)), 'n',n);
    end
    E.(fld) = Q;
end
E.featNames = fieldnames(E.mM8)';
E.built = datestr(now,'yyyy-mm-dd HH:MM:SS');                          %#ok<TNOW1,DATST>

save(fullfile(outDir,'envelope_bothatp.mat'),'E');
fprintf('\n%d features -> %s\n', numel(E.featNames), fullfile(outDir,'envelope_bothatp.mat'));

%% ---- the ATP effect, per feature, for the record -----------------------
fprintf('\n%-22s %10s %10s %8s %8s\n','feature','8 mM mean','2 mM mean','2/8','z(03/27)');
for j = 1:numel(E.featNames)
    nm = E.featNames{j};
    if ~isfield(E.mM2,nm); continue; end
    a = mean(E.mM8.(nm).mean,'omitnan'); b = mean(E.mM2.(nm).mean,'omitnan');
    sd = mean(E.mM8.(nm).sd,'omitnan');
    z  = (mean(E.mM2.(nm).raw(IREF,:),'omitnan') - a) / max(sd,eps);
    fprintf('%-22s %10.4g %10.4g %8.3f %8.2f\n', nm, a, b, b/a, z);
end
