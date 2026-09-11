function [vt, meta] = buildSawtoothProtocol(L0, fracAmp, v, tWarm, tPost, ML)
%BUILDSAWTOOTHPROTOCOL Symmetric stretch-release triangle as a velocitytable.
%
%   [vt, meta] = buildSawtoothProtocol(L0, fracAmp, v, tWarm, tPost, ML)
%
%   Ports the length protocol of van der Zee et al. 2026 (PLOS Comput Biol,
%   10.1371/journal.pcbi.1014748): an isokinetic stretch of fracAmp * L0 at
%   speed v, immediately reversed at the SAME speed back to L0 - a sawtooth.
%
%   Inputs
%     L0      starting length (ML). With params.ML = 2 um, L0 = 1.0 <-> SL 2.0 um.
%     fracAmp stretch amplitude as a fraction of L0 (0.04 = 4%).
%     v       ramp speed (ML/s, positive). Ramp DURATION is derived from
%             amplitude and speed - never typed as a third, redundant number.
%     tWarm   isometric hold before the stretch (s). Lets the cross-bridge
%             population reach steady state so the response is to the stretch
%             and not to the initial condition.
%     tPost   isometric hold after the release (s), for force recovery.
%     ML      um per muscle length, for the bookkeeping column 3.
%
%   Output
%     vt      [ t(s), v(ML/s), v_um(bookkeeping), L(ML) ]. The model reads
%             COLUMN 2 ONLY, held constant over [t_i, t_i+1]; columns 3-4 are
%             bookkeeping. t = 0 is stretch onset by convention.
%     meta    struct of the derived timings and amplitudes.
%
%   See also Example/BuildExampleProtocol.m (same pattern, one-way stretch).

    dL = fracAmp * L0;      % ML
    dt = dL / v;            % s - derived, so column 4 cannot drift from column 2

    t0 = 0;                 % stretch onset
    t1 = t0 + dt;           % apex: stretch ends, release begins
    t2 = t1 + dt;           % back at L0: release ends
    t3 = t2 + tPost;        % end of recovery hold

    vt = [ -tWarm,  0,  0,      L0
            t0,    +v, +v*ML,   L0        % STRETCH
            t1,    -v, -v*ML,   L0 + dL   % RELEASE (same speed, reversed)
            t2,     0,  0,      L0        % recovery hold
            t3,     0,  0,      L0 ];     % terminal row: velocity never used

    meta = struct('L0', L0, 'fracAmp', fracAmp, 'v', v, 'dL', dL, 'dt', dt, ...
                  't0', t0, 't1', t1, 't2', t2, 't3', t3, ...
                  'tWarm', tWarm, 'tPost', tPost, 'ML', ML, ...
                  'dSL', dL * ML, 'nm_per_HS', dL * ML * 1000 / 2);
end
