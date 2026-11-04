% AttachmentPointerArtefact.m
% The attachment source in dPUdT_CombinedTransitions (UseA1AttachmentKernel = 0,
% A1AttachmentWidth = Ax > 0 -> the "else" branch) is a triangular kernel of
% half-width Ax sampled on the strain grid around the moving zero. As SL - LSE
% changes, the zero slides across the bins and the SAMPLED kernel's moments
% change. A moments model would use the continuous triangle exactly:
%   mean 0, variance Ax^2/6, 3rd moment 0, 4th Ax^4/15.
% This script quantifies how far the grid version departs from that.

params0 = momentsSetup('fastgate_best');
dS = params0.dS; Ax = params0.A1AttachmentWidth;
offs = linspace(0, 1, 201) * dS;                 % pointer position within one bin
M = zeros(numel(offs), 4);
for k = 1:numel(offs)
    s = ((-30:30)' * dS) - offs(k);              % zero sits offs(k) past a bin
    % --- verbatim logic of the dxdt else-branch ---
    s_0 = 1 + round(-s(1)/dS, 6);
    halfspan = ceil(Ax/dS);
    att = zeros(size(s));
    for jj = floor(max(1, s_0-halfspan)):ceil(min(numel(s), s_0+halfspan))
        att(jj) = max(0, 1 - abs(s(jj))/Ax);
    end
    att_norm = att/max(1, sum(att))/dS;
    % ---
    m = att_norm * dS;                            % masses per unit RD1
    M(k,:) = [sum(m), sum(s.*m), sum(s.^2.*m), sum(s.^3.*m)];
end
fprintf('Attachment kernel moments as the pointer slides through one bin (Ax = %.1f nm, dS = %.1f nm):\n', Ax*1e3, dS*1e3);
fprintf('  mass        %.4f .. %.4f   (exact 1)\n', min(M(:,1)), max(M(:,1)));
fprintf('  mean (nm)   %+.3f .. %+.3f (exact 0)\n', 1e3*min(M(:,2)), 1e3*max(M(:,2)));
fprintf('  sd (nm)     %.3f .. %.3f   (exact %.3f)\n', 1e3*sqrt(min(M(:,3)-M(:,2).^2)), 1e3*sqrt(max(M(:,3)-M(:,2).^2)), 1e3*Ax/sqrt(6));
K = xbKernels(((-30:30)'*dS), params0, 0);
k12ofs = zeros(numel(offs),1);
for k = 1:numel(offs)
    s = ((-30:30)' * dS) - offs(k); Kk = xbKernels(s, params0, 0);
    s_0 = 1 + round(-s(1)/dS, 6); att = zeros(size(s));
    for jj = floor(max(1, s_0-ceil(Ax/dS))):ceil(min(numel(s), s_0+ceil(Ax/dS))), att(jj) = max(0, 1 - abs(s(jj))/Ax); end
    m = att/max(1,sum(att));
    k12ofs(k) = sum(m .* Kk.k12);                 % initial stroke rate of a fresh cohort
end
xs = linspace(-Ax, Ax, 4001)'; tri = max(0, 1-abs(xs)/Ax); tri = tri/trapz(xs, tri);
k12exact = trapz(xs, tri .* xbKernels(xs, params0, 0).k12);
fprintf('  mean k12 seen by a fresh cohort: %.0f .. %.0f 1/s (continuous triangle %.0f 1/s)\n', min(k12ofs), max(k12ofs), k12exact);
save(fullfile(fileparts(mfilename('fullpath')), 'results', 'attachPointer.mat'), 'offs', 'M', 'k12ofs', 'k12exact');
