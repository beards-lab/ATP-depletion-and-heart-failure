function target = buildAtpRatioTarget(f2d, f8d, mode)
%BUILDATPRATIOTARGET  Build the low-ATP relative-ratio target from 2 mM & 8 mM data.
%   target = buildAtpRatioTarget(f2d, f8d)            % 'full' (default)
%   target = buildAtpRatioTarget(f2d, f8d, 'sparse')  % legacy 8-feature set
%
%   'full'   : EVERY well-behaved slack feature scored PER SEGMENT (5-vec), mirroring the
%              high-ATP feature selector list. This constrains the whole restretch timecourse
%              shape, not just segment-averaged amplitudes.
%   'sparse' : the legacy set (steady/A/Am/peak2/peak1_y/vall_y as scalar means + ktr/t0
%              per-segment) — kept for continuity with the 0.1188 baseline.
%
%   All comparisons are ratios (2 mM / 8 mM), so the 8 mM baseline's own fit error cancels and
%   only the [ATP]-driven CHANGE is scored. Cost is evaluated by atpRatioCost.
    if nargin < 3 || isempty(mode); mode = 'full'; end
    target = struct('scalar',struct(),'vector',struct(),'weights',struct());

    switch lower(mode)
      case 'sparse'
        sf = {'steady','A','Am','peak2','peak1_y','vall_y'};   % scalar means
        vf = {'ktr','t0'};                                     % per-segment
        w  = struct('steady',50,'A',50,'Am',10,'peak2',5,'peak1_y',10,'vall_y',10,'ktr',2,'t0',1);
        for i=1:numel(sf); f=sf{i}; target.scalar.(f)=mean(double(f2d.(f))./double(f8d.(f)),'omitnan'); end
        for i=1:numel(vf); f=vf{i}; target.vector.(f)=double(f2d.(f))./double(f8d.(f)); end

      case 'full'
        % feature -> weight (mirrors the high-ATP fn list; FV / bound-regularisers dropped,
        % as there is no 2 mM force-velocity data to ratio against).
        vw = { 'steady',50; 'A',50; 'Am',10; 'peak1_y',10; 'vall_y',10; 'peak2',5; ...
               'ktr',2; 't0',2; 'peak1_dSL',1; 'restretchSlopeStart',0.5; ...
               'vall2_dy',0.5; 'ovrsht_dy',0.3 };
        for i=1:size(vw,1)
            f = vw{i,1};
            if ~isfield(f2d,f) || ~isfield(f8d,f); continue; end
            r = double(f2d.(f)) ./ double(f8d.(f));
            target.vector.(f) = r;
            w.(f) = vw{i,2};
        end
      otherwise
        error('buildAtpRatioTarget:mode','unknown mode "%s"', mode);
    end
    target.weights = w;
end
