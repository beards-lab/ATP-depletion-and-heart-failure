% WritePooledTarget.m -- freeze buildPooledTarget() into a version-controlled
% target the model eval path can read without the raw data.
%
% WHY A GENERATED FILE. `buildPooledTarget` needs the six curated protocol
% files, and `data/` is not in git. The objective must not depend on files that
% are not tracked, and it must not re-extract features from raw traces once per
% optimiser evaluation. So the pooled numbers are frozen ONCE, here, into
%
%     Model/pooledTarget8mM.m
%
% a plain function returning the struct. Regenerate it whenever a new
% preparation joins the pool -- that is the whole point of the exercise: the
% target moves as data comes in, and the file's header records which
% preparations and which options produced the numbers in it.
%
% UNITS. buildPooledTarget divides every force-dimension feature by its own
% preparation's pre-slack isometric at SL 2.2 (`Fss22`), which removes the
% absolute scale. The model, however, produces kPa. So the force features are
% multiplied back by the ONE pooled level `Fref = mean(Fss22)` -- i.e. the
% target handed to the objective is "the pooled SHAPE at the pooled LEVEL".
% Scale-free features (rates, times, dSL) pass through untouched.
%
% Output -> Model/pooledTarget8mM.m
%           results/pooled_target_report.txt   (old single-prep vs new pooled)

clearvars -except OPTS OUTF; clc;
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..');
cd(root); addpath(genpath(root));
resDir = fullfile(here, 'results'); if ~isfolder(resDir); mkdir(resDir); end

% OUTF may be set by the caller, to write somewhere other than Model/ first
% (the file is written in place, so a campaign reading it mid-write would see
% half a target -- stage, verify, then move it in).
if ~exist('OUTF','var') || isempty(OUTF); OUTF = fullfile(root, 'Model', 'pooledTarget8mM.m'); end
REF_FILE = 'protocol_03_27_2026_8mM_slack.mat';   % the single-prep target it replaces

if ~exist('OPTS','var'); OPTS = struct(); end     % passed through to buildPooledTarget

%% ---- the pooled target -------------------------------------------------
T    = buildPooledTarget(OPTS);
Fref = T.scale.mean;
nF   = numel(T.name);
isF  = strcmp(T.cls, 'F');

val = T.m8;
val(isF, :) = val(isF, :) * Fref;                 % back into kPa (or kPa/s)

%% ---- what it replaces --------------------------------------------------
% The old target is whatever `runSlackExperiment` reads out of the curated
% file, in that preparation's own kPa. Printing the two side by side is the
% only honest way to see how far the target itself moved -- a cost change
% after this switch is NOT a fit change until this table is read.
S   = load(fullfile('data', REF_FILE));
old = nan(nF, size(val,2));
if isfield(S, 'features_data')
    for k = 1:nF
        if ~isfield(S.features_data, T.name{k}); continue; end
        v = double(S.features_data.(T.name{k})(:))';
        m = min(size(val,2), numel(v));
        old(k, 1:m) = v(1:m);
    end
end

lines = {};
lines{end+1} = sprintf('pooled 8 mM target -- %d preparations (%s), anchor %s', ...
    T.nP, strjoin(T.scale.days, ', '), T.opts.anchor);
lines{end+1} = sprintf('Fref = %.2f kPa (CV %.1f %%, per prep %s), moderated dof = %.2f', ...
    Fref, 100*T.scale.cv, strtrim(sprintf('%.1f ', T.scale.perPrep)), T.d0 + T.nP - 1);
lines{end+1} = '';
lines{end+1} = sprintf('%-22s %4s %12s %12s %9s %9s %9s', ...
    'feature', 'cls', 'old (03/27)', 'new (pooled)', 'change', 's (used)', 'w=1/s^2');
lines{end+1} = repmat('-', 1, 86);
for k = 1:nF
    mo = mean(old(k,:), 'omitnan');
    mn = mean(val(k,:), 'omitnan');
    if isfinite(mo) && mo ~= 0
        chg = sprintf('%+8.1f%%', 100*(mn/mo - 1));
    else
        chg = sprintf('%9s', '-');
    end
    lines{end+1} = sprintf('%-22s %4s %12.4g %12.4g %s %8.1f%% %9.1f', ...
        T.name{k}, T.cls{k}, mo, mn, chg, 100*T.s8(k), 1/T.s8(k)^2);  %#ok<SAGROW>
end
lines{end+1} = repmat('-', 1, 86);
lines{end+1} = ['NOTE the force block moves as ONE number: every ''F'' feature is ' ...
                'rescaled by'];
lines{end+1} = sprintf(['     Fref / Fss22(03/27) = %.2f / %.2f = %.3f, so a uniform ' ...
                        '%+.1f %% on force.'], Fref, T.scale.perPrep(1), ...
                        Fref/T.scale.perPrep(1), 100*(Fref/T.scale.perPrep(1)-1));
lines{end+1} = '     Anything beyond that factor is a SHAPE difference between 03/27 and the pool.';

txt = strjoin(lines, newline);
disp(txt);
fid = fopen(fullfile(resDir, 'pooled_target_report.txt'), 'w');
fprintf(fid, '%s\n', txt); fclose(fid);

%% ---- serialise ---------------------------------------------------------
% %.17g throughout: the objective's numerical roughness is ~0.076 (see
% Docs objective_noise_floor), so a target that does not round-trip exactly
% would be indistinguishable from search noise and impossible to audit.
fid = fopen(OUTF, 'w');
fprintf(fid, 'function PT = pooledTarget8mM()\n');
fprintf(fid, '%% POOLEDTARGET8MM  Frozen 8 mM pooled fit target (GENERATED -- do not hand-edit).\n%%\n');
fprintf(fid, '%%   Written by Analyses/DataPooling/WritePooledTarget.m on %s\n', ...
             char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')));
fprintf(fid, '%%   from %d preparations: %s\n', T.nP, strjoin(T.scale.days, ', '));
fprintf(fid, '%%   anchor %s | rundown %d | shrink %d | phi %.4g\n', ...
             T.opts.anchor, T.opts.rundown, T.opts.shrink, T.opts.phi);
fprintf(fid, '%%\n');
fprintf(fid, '%% Fields\n');
fprintf(fid, '%%   name/cls  feature name and class ("F" force-dimension, "0" scale-free)\n');
fprintf(fid, '%%   val       pooled 8 mM value, ABSOLUTE units, one column per slack segment\n');
fprintf(fid, '%%   s/sraw    relative between-preparation SD, moderated / raw (the chi-2 scale)\n');
fprintf(fid, '%%   seg       which slack cycles s was measured on. The objective MUST score the\n');
fprintf(fid, '%%             feature on those cycles and no others -- a chi-2 is only a chi-2\n');
fprintf(fid, '%%             if the residual and the spread come from the same manoeuvre.\n');
fprintf(fid, '%%   r/sr/z    2 mM / 8 mM ratio, its spread, and |r-1|/sqrt(s^2+sr^2)\n');
fprintf(fid, '%%   Fref      the pooled absolute level in kPa, and the range the preps span\n');
fprintf(fid, '%%   dof       degrees of freedom behind s (d0 + n - 1), for confidence intervals\n');
fprintf(fid, '%%\n');
fprintf(fid, '%% Regenerate after every new preparation. Consumed by pooledTargetFeatures.\n\n');
fprintf(fid, 'persistent P\nif ~isempty(P); PT = P; return; end\n\n');
fprintf(fid, 'P.name = {%s};\n', strjoin(cellfun(@(s) ['''' s ''''], T.name, 'uni', 0), ', '));
fprintf(fid, 'P.cls  = {%s};\n', strjoin(cellfun(@(s) ['''' s ''''], T.cls,  'uni', 0), ', '));
fprintf(fid, 'P.val  = [ ...\n');
for k = 1:nF
    fprintf(fid, '    %s; ... %% %s\n', strtrim(sprintf('%.17g ', val(k,:))), T.name{k});
end
fprintf(fid, '    ];\n');
writeRow(fid, 'P.s',    T.s8);
writeRow(fid, 'P.sraw', T.s8raw);
fprintf(fid, 'P.seg      = {%s};\n', strjoin(cellfun(@(v) mat2str(v), T.seg, 'uni', 0), ', '));
writeRow(fid, 'P.r',    T.r);
writeRow(fid, 'P.sr',   T.sr);
writeRow(fid, 'P.z',    T.z);
fprintf(fid, 'P.Fref     = %.17g;\n', Fref);
fprintf(fid, 'P.FrefCV   = %.17g;\n', T.scale.cv);
fprintf(fid, 'P.FrefPrep = [%s];\n', strtrim(sprintf('%.17g ', T.scale.perPrep)));
fprintf(fid, 'P.level    = %.17g;   %% 2 mM / 8 mM force level\n', T.level);
fprintf(fid, 'P.sLevel   = %.17g;\n', T.sLevel);
fprintf(fid, 'P.dof      = %.17g;\n', T.d0 + T.nP - 1);
fprintf(fid, 'P.nP       = %d;\n', T.nP);
fprintf(fid, 'P.days     = {%s};\n', strjoin(cellfun(@(s) ['''' s ''''], T.scale.days, 'uni', 0), ', '));
fprintf(fid, '\nPT = P;\nend\n');
fclose(fid);

fprintf('\nwrote %s  (%d features x %d segments)\n', OUTF, nF, size(val,2));
fprintf('report -> %s\n', fullfile(resDir, 'pooled_target_report.txt'));

function writeRow(fid, nm, v)
    fprintf(fid, '%-10s = [%s];\n', nm, strtrim(sprintf('%.17g ', v)));
end
