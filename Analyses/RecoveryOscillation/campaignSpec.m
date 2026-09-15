function S = campaignSpec(name)
%CAMPAIGNSPEC  What a named campaign starts from, and what it writes.
%
%   S = campaignSpec('floored')
%
% One place holding the mapping, because five things need it and agreeing by
% hand is how a report ends up consolidating one campaign against another
% campaign's seed. The launcher, the SLURM job, the campaign script, the
% consolidation report and campaignStatus all read it from here.
%
% S.stem        tag prefix; the profiles are <stem>_a .. <stem>_d
% S.seedSnap    params/<...>.m the search starts from
% S.seedTag     that snapshot's tag, for the report's baseline row
% S.seedExpect  its cost HERE. RunOptimSpreadW refuses a seed more than 0.25
%               away -- the objective's own peak-to-peak roughness, not a
%               tolerance to tighten (Docs objective_noise_floor). Note that
%               the Great Lakes nodes read this objective ~+0.015 high
%               (different MATLAB release); that offset is inside the noise
%               and cancels within a machine, but do not quote absolute costs
%               across machines without it.
% S.tags        the four profile tags
% S.note        one line for the log, so a log says which campaign it is
%
% THE CAMPAIGNS
%
%   rigaware  the rsK80[1:4] objective (conclusions.md sec. 4). RAN, 14 h,
%             and put 110 % of its gain into force-velocity (sec. 7). Kept so
%             its result is reproducible -- do not relaunch it.
%   floored   the same objective with every data term given a one-
%             preparation-width deadband (sec. 9), seeded from the rigaware
%             campaign's own best point so its F-V gain is kept and not paid
%             for again. NEVER RAN -- the compute went to MechanismAblation,
%             which used floored_seed as its base. Superseded by `lean5`.
%   lean5     `floored` plus the three things that campaign taught us:
%             slack cycle 5 is scored again (it is a DIFFERENT manoeuvre --
%             the only re-stretch to ~2.0 um SL -- and carries a real ATP
%             signal, z = 3.01); the draw pool loses the twelve parameters
%             MechanismAblation measured as inert, deleted or redundant and
%             gains `xrate`; and coolDownLS is retired in favour of rsK80 x2.
%             Seed by BuildLean5Seed.m. THIS IS THE ONE TO RUN.
%
% A NOTE ON COMPARING COSTS ACROSS CAMPAIGNS: don't. Each seed is the previous
% optimum reweighted with lambda chosen for cost-neutrality, so the numbers
% reset at every generation. Only within a campaign is a drop a fit improvement.

if nargin < 1 || isempty(name); name = 'lean5'; end
name = lower(char(name));

switch name
    case 'rigaware'
        S.stem       = 'rigaware';
        S.seedTag    = 'rigaware_seed';
        S.seedExpect = 2.6858;
        S.note       = 'rsK80[1:4], unfloored -- SUPERSEDED, see conclusions.md sec. 7';
    case 'floored'
        S.stem       = 'floored';
        S.seedTag    = 'floored_seed';
        S.seedExpect = 0.9081;
        S.note       = 'rsK80[1:4] + one-width deadband -- NEVER RAN, superseded by lean5';
    case 'lean5'
        S.stem       = 'lean5';
        S.seedTag    = 'lean5_seed';
        % Measured by the pre-flight on THIS machine, not copied from the
        % snapshot header (which reads 2.3167 -- the UNFLOORED cost-neutral
        % total). Reading an expected cost off a header is how a campaign ends
        % up refusing a good seed; see Docs objective_noise_floor.
        S.seedExpect = 0.9800;
        S.note       = ['cycle 5 scored again + lean pool (12 out, xrate in) ' ...
                        '+ coolDownLS retired, rsK80 x2, FV pin 0.25'];
    otherwise
        error('campaignSpec:unknown', ...
              ['unknown campaign "%s" (expected ''rigaware'', ''floored'' ' ...
               'or ''lean5'')'], name);
end

S.name     = name;
S.seedSnap = fullfile('params', [S.seedTag '.m']);
S.tags     = strcat(S.stem, {'_a','_b','_c','_d'});
end
