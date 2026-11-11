% RunNoiseBudget.m -- how much of the "not exponential" in a recovery window
% is measurement noise, how much is a reproducible ripple, and how much is a
% real shape the model should be chasing?
%
% The claim under test is that every force recovery is a single exponential
% plus something the rig contributes. If that is right, then the residual of
% a one-exponential fit should decompose into
%
%     var(res) = sigma_w^2            white transducer noise
%              + sigma_c^2            reproducible, phase-locked to the ramp
%              + sigma_s^2            genuine shape the model is missing
%
% and only the third term is worth a weight. The three are separated without
% any model of the ripple:
%
%   sigma_w   from the flat plateau of the SAME window (last 40 % of the
%             span, linearly detrended, high-pass above 60 Hz) -- this is the
%             transducer, measured in situ, with no cross-bridge content.
%   sigma_c   by COHERENT AVERAGING the five slack cycles of a run on their
%             common post-ramp clock. Anything phase-locked to the ramp
%             survives; white noise falls as 1/sqrt(5). The excess of the
%             average's variance over sigma_w^2/5 IS the locked part.
%   sigma_s   whatever is left.
%
% The decisive control is the PNB+Mava run: the same rig, the same ramps, the
% cross-bridges chemically switched off.
%
% Output -> results/noise_budget.mat / .txt / .png

clearvars; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here,'..','..'); cd(root); addpath(genpath(root));
resDir = fullfile(here,'results'); if ~isfolder(resDir); mkdir(resDir); end

Q = load(fullfile(resDir,'oscillation_id.mat')); W = Q.W;
days = unique({W.day},'stable'); conds = {'8mM','2mM','PNB'};
TY   = {'postRestretch','postSlack'};

B = struct([]);
for d = 1:numel(days)
 for c = 1:numel(conds)
  for ty = 1:numel(TY)
    s = find(strcmp({W.day},days{d}) & strcmp({W.cond},conds{c}) & strcmp({W.type},TY{ty}));
    if numel(s) < 3; continue; end
    fs = W(s(1)).fs;
    n  = min(arrayfun(@(k) numel(W(k).res), s));
    Rm = cell2mat(arrayfun(@(k) W(k).res(1:n), s, 'UniformOutput', false));
    Zm = cell2mat(arrayfun(@(k) W(k).z(1:n),   s, 'UniformOutput', false));
    nc = size(Rm,2);

    % ---- sigma_w : the transducer, from the plateau of each window -------
    tt   = W(s(1)).t(1:n);
    plat = tt >= 0.60*tt(end);
    sw = zeros(1,nc);
    for j = 1:nc
        y = detrend(Rm(plat,j));
        sw(j) = std(highpass60(y, fs));           % above 60 Hz: no biology
    end
    % a first-difference estimator as a cross-check, immune to any trend
    swD = std(diff(Rm(plat,:)))/sqrt(2);
    sigW = mean(sw);

    % ---- sigma_c : phase-locked to the ramp ------------------------------
    mR   = mean(Rm,2);
    varA = var(mR);                                % variance of the average
    varC = max(varA - sigW^2/nc, 0);               % excess over the noise floor
    sigC = sqrt(varC);
    % significance: how many sd of the null (white) is the excess?
    nullSD = sigW^2/nc * sqrt(2/(n-1));
    zC     = (varA - sigW^2/nc) / max(nullSD, eps);

    % ---- the total, and what is left -------------------------------------
    sigTot = mean(std(Rm));
    sigS   = sqrt(max(sigTot^2 - sigW^2 - sigC^2, 0));

    A   = mean([W(s).A], 'omitnan');
    r2  = mean([W(s).r2]);
    % R2 a NOISE-FREE model could reach against this data
    r2ceil = 1 - (sigW^2 + sigC^2)/mean(var(Zm));

    % dominant period of the locked part
    ra = movmean(mR, max(3,round(0.002*fs)));
    [px,fx] = periodogram(ra-mean(ra), hann(numel(ra)), 8192, fs);
    ok = fx > 1/tt(end) & fx < 300; fo = fx(ok); [~,jj] = max(px(ok));
    perC = 1e3/fo(jj);

    b = struct('day',days{d},'cond',conds{c},'type',TY{ty},'nc',nc, ...
        'sigW',sigW,'sigWdiff',mean(swD),'sigC',sigC,'sigS',sigS,'sigTot',sigTot, ...
        'zC',zC,'A',A,'r2',r2,'r2ceil',r2ceil,'perC',perC,'mR',mR,'t',tt, ...
        'fracW',sigW^2/sigTot^2,'fracC',sigC^2/sigTot^2,'fracS',sigS^2/sigTot^2);
    if isempty(B); B = b; else; B(end+1) = b; end %#ok<AGROW>
  end
 end
end

save(fullfile(resDir,'noise_budget.mat'),'B');
reportNoiseBudget(B, resDir);
plotNoiseBudget(B, resDir);

function y = highpass60(x, fs)
    n = max(3, round(fs/60));                 % 60 Hz corner as a moving mean
    y = x - movmean(x, n);
end
