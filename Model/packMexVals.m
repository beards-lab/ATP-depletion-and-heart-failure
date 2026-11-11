function mv = packMexVals(params)
% PACKMEXVALS  Pack params struct into flat double array for MEX RHS calls.
%
%   MV = PACKMEXVALS(PARAMS) returns a 187-element column vector containing
%   all numerical parameters needed by dPUdT_core_mex.
%
%   Layout (1-based indices):
%     1-43  : scalar rate/geometry parameters (see IDX_* constants below)
%     44-103: strain grid params.s  (ss=60 elements)
%     104-107: pchip R1D breaks (4 values, n=3 segments)
%     108-119: pchip R1D coefs  (3×4 col-major = 12 values)
%     120-124: pchip pwsd breaks (5 values, n=4 segments)
%     125-140: pchip pwsd coefs  (4×4 col-major = 16 values)
%     141-144: pchip R21 breaks (4 values, n=3 segments)
%     145-156: pchip R21 coefs  (3×4 col-major = 12 values)
%     157-163: pchip pwsd2 breaks (7 values, n=6 segments)
%     164-187: pchip pwsd2 coefs  (6×4 col-major = 24 values)
%
%   IMPORTANT: This file is config-specific (hard-codes the active switch
%   set). Regenerate when params.mods or active feature flags change.
%   The index layout here must match the #define constants in dPUdT_core_mex.c.
%
%   Prerequisites:
%     params.cached.pwsdR1D_brk/cof, pwsd_brk/cof, pwsdR21_brk/cof,
%     pwsd2_brk/cof must be pre-populated (set UseFastPPval=true in getParams).
%
%   See also: getParams, dPUdT_CombinedTransitions_mex, compileMex

% --- verify pchip cache exists ---
assert(isfield(params,'cached') && isfield(params.cached,'pwsd_brk'), ...
    'packMexVals: params.cached pchip data missing. Set UseFastPPval=true.');

% --- scalars (43 total) ---
scalars = [
    params.kd;              %  1
    params.k1;              %  2
    params.k_1;             %  3
    params.k2;              %  4
    params.kstiff1;         %  5
    params.kstiff2;         %  6
    params.kstiff1_n;       %  7
    params.kstiff2_n;       %  8
    params.ka;              %  9
    params.kah;             % 10
    params.kamh;            % 11
    params.dS;              % 12
    params.ss;              % 13
    params.dr;              % 14
    params.dr2;             % 15
    params.L_thick;         % 16
    params.L_hbare;         % 17
    params.L_thin;          % 18
    params.k_pas;           % 19
    params.Lsc0;            % 20
    params.gamma;           % 21
    params.kSE;             % 22
    params.ekSE;            % 23
    params.MaxSlackNegativeForce;  % 24
    params.mu;              % 25
    params.mu_neg;          % 26
    params.mu2;             % 27
    params.ksr0;            % 28
    params.sigma1;          % 29
    params.sigma2;          % 30
    params.kmsr;            % 31
    params.kmsrd;           % 32
    params.sigma_srd1;      % 33
    params.ksrd;            % 34
    params.sigma_srd2;      % 35
    params.ksrd2sr;         % 36
    params.ksr2srd;         % 37
    params.slope;           % 38
    params.d_actin;         % 39
    params.s_threshold_R;   % 40
    params.s_threshold_L;   % 41
    params.A1AttachmentWidth; % 42
    params.estiff;          % 43
];

% --- strain grid (padded to N_SS_MAX=60 elements, must be column) ---
% params.s may have fewer entries when Slim_l/Slim_r produce a narrow range
% (e.g. FV conditions with ss=11). Pad with zeros so pchip offsets stay fixed.
% The C MEX reads the actual count from mv(13) = params.ss and loops only that many times.
N_SS_MAX = 60;
s_grid_raw = params.s(:);
assert(numel(s_grid_raw) <= N_SS_MAX, ...
    'packMexVals: ss=%d exceeds N_SS_MAX=%d. Increase N_SS_MAX in packMexVals and recompile MEX.', ...
    numel(s_grid_raw), N_SS_MAX);
s_grid = zeros(N_SS_MAX, 1);
s_grid(1:numel(s_grid_raw)) = s_grid_raw;  % 44..103 (padded)

% --- pchip data ---
% R1D: n=3 segments, breaks=4, coefs=3×4
brk_R1D = params.cached.pwsdR1D_brk(:);   % 4 values
cof_R1D = params.cached.pwsdR1D_cof(:);   % 12 values (col-major 3×4)

% pwsd (R12): n=4 segments, breaks=5, coefs=4×4
brk_pwsd = params.cached.pwsd_brk(:);     % 5 values
cof_pwsd = params.cached.pwsd_cof(:);     % 16 values (col-major 4×4)

% R21: n=3 segments, breaks=4, coefs=3×4
brk_R21 = params.cached.pwsdR21_brk(:);   % 4 values
cof_R21 = params.cached.pwsdR21_cof(:);   % 12 values (col-major 3×4)

% pwsd2 (R2): n=6 segments, breaks=7, coefs=6×4
brk_pwsd2 = params.cached.pwsd2_brk(:);   % 7 values
cof_pwsd2 = params.cached.pwsd2_cof(:);   % 24 values (col-major 6×4)

% --- assemble ---
mv = [scalars; s_grid;
      brk_R1D; cof_R1D;
      brk_pwsd; cof_pwsd;
      brk_R21; cof_R21;
      brk_pwsd2; cof_pwsd2];

% Sanity check: expected total length
% 43 + 60 (padded s_grid) + (4+12) + (5+16) + (4+12) + (7+24) = 187
% s_grid is always padded to 60 regardless of actual params.ss (which may be < 60
% for narrow-range conditions like FV). Actual ss is in mv(13) = params.ss.
expected_len = 187;
if numel(mv) ~= expected_len
    error('packMexVals: expected %d elements but got %d. Check pchip segment counts.', ...
        expected_len, numel(mv));
end
end
