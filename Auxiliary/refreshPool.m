function refreshPool(nThreads)
%REFRESHPOOL Restart the parallel pool so workers pick up current code.
%
%   Thread-based parpool workers ('Threads') cache function code in memory.
%   After editing a Model function mid-session, existing workers keep running
%   the OLD code -- there is no automatic invalidation, and functions
%   dispatched via parfeval (e.g. Model/experiments/runSlackExperiment.m) will
%   silently use stale logic until the pool is restarted.
%
%   REFRESHPOOL deletes any existing pool and starts a fresh 'Threads' pool,
%   and calls rehash so the client session also sees any new/changed files.
%
%   Usage:
%       refreshPool();      % default 5 threads
%       refreshPool(8);
%
%   Call this at the top of run scripts, after addpath/genpath, any time you
%   may have edited Model code since the pool was last started.
%
%   See also: parpool, gcp, rehash

    if nargin < 1 || isempty(nThreads)
        nThreads = 5;
    end

    try
        p = gcp('nocreate');
        if ~isempty(p)
            delete(p);
        end
        parpool('Threads', nThreads);
    catch ME
        warning('refreshPool:noParallel', ...
            'Could not (re)start parallel pool (%s). Continuing without it.', ME.message);
    end

    rehash;
end
