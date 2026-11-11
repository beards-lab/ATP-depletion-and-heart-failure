function [vt, meta] = buildConditioningTestProtocol(L0, fracAmp, v, tWarm, tRec, tPost, ML)
%BUILDCONDITIONINGTESTPROTOCOL Two identical sawtooths separated by a recovery hold.
%
%   [vt, meta] = buildConditioningTestProtocol(L0, fracAmp, v, tWarm, tRec, tPost, ML)
%
%   The full structure of van der Zee et al. 2026 (PLOS Comput Biol,
%   10.1371/journal.pcbi.1014748): isometric phase, isokinetic CONDITIONING
%   stretch, isokinetic shortening back to the pre-conditioning length,
%   isometric recovery phase, isokinetic TEST stretch.
%
%   Both events are made IDENTICAL sawtooths (stretch then release at the same
%   speed) rather than conditioning-sawtooth-then-bare-test-stretch. That costs
%   nothing and buys the thing the comparison needs: the two ramps see exactly
%   the same length input, so any difference in their onset stiffness is history
%   and nothing else.
%
%   Inputs are as buildSawtoothProtocol, plus
%     tRec   isometric recovery hold between the two events (s).
%
%   Output
%     vt     [ t(s), v(ML/s), v_um(bookkeeping), L(ML) ]; the model reads
%            COLUMN 2 ONLY. t = 0 is the CONDITIONING stretch onset.
%     meta   timings, with meta.onsets = [t0 t3] the two stretch onsets.
%
%   See also buildSawtoothProtocol (the single-event version).

    dL = fracAmp * L0;      % ML
    dt = dL / v;            % s - derived, never typed as a third number

    t0 = 0;                 % conditioning stretch onset
    t1 = t0 + dt;           % apex 1
    t2 = t1 + dt;           % back at L0
    t3 = t2 + tRec;         % TEST stretch onset
    t4 = t3 + dt;           % apex 2
    t5 = t4 + dt;           % back at L0
    t6 = t5 + tPost;

    vt = [ -tWarm,  0,  0,      L0
            t0,    +v, +v*ML,   L0        % CONDITIONING stretch
            t1,    -v, -v*ML,   L0 + dL   % release
            t2,     0,  0,      L0        % recovery hold
            t3,    +v, +v*ML,   L0        % TEST stretch
            t4,    -v, -v*ML,   L0 + dL   % release
            t5,     0,  0,      L0        % final hold
            t6,     0,  0,      L0 ];     % terminal row: velocity never used

    meta = struct('L0', L0, 'fracAmp', fracAmp, 'v', v, 'dL', dL, 'dt', dt, ...
                  't0', t0, 't1', t1, 't2', t2, 't3', t3, 't4', t4, ...
                  't5', t5, 't6', t6, 'onsets', [t0 t3], ...
                  'tWarm', tWarm, 'tRec', tRec, 'tPost', tPost, 'ML', ML, ...
                  'dSL', dL * ML, 'nm_per_HS', dL * ML * 1000 / 2);
end
