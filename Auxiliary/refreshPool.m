function refreshPool(n)
%REFRESHPOOL  Clear stale compiled code and (re)start a fresh parallel pool.
%
%   REFRESHPOOL() clears the client-side compiled-code cache so that edits to
%   model files (e.g. dPUdT_CombinedTransitions.m or a param snapshot) take
%   effect, then deletes any existing parallel pool and starts a fresh one with
%   the default local profile size. Parpool workers keep their own copy of
%   compiled code, so without restarting them, edits made during a session are
%   silently ignored — a recurring foot-gun in this project.
%
%   REFRESHPOOL(N) starts the fresh pool with N workers.
%   REFRESHPOOL(0) clears the code cache and shuts the pool down without
%   starting a new one.
%
%   Requires the Parallel Computing Toolbox for the pool steps; without it, the
%   code cache is still cleared and the function returns quietly.
%
%   See also: loadParams, parpool, gcp, pctRunOnAll

if nargin < 1 || isempty(n)
    n = -1;   % sentinel: use the default local pool size
end

% (1) Clear client-side compiled code so edits to .m files take effect.
clear functions %#ok<CLFUNC>
rehash

% No Parallel Computing Toolbox -> nothing further to do.
if exist('parpool', 'file') ~= 2
    return;
end

% (2) Delete any existing pool; its workers cache their own stale code.
existing = gcp('nocreate');
if ~isempty(existing)
    delete(existing);
end

if isequal(n, 0)
    return;   % caller asked for no pool
end

% (3) Start a fresh *process*-based pool. Process pools (not thread pools) are
%     required for UseParallel optimizers, which silently no-op on a thread pool.
try
    if n < 0
        parpool('local');
    else
        parpool('local', n);
    end
catch err
    warning('refreshPool:parpoolFailed', ...
        'Could not start parpool: %s', err.message);
end
end
