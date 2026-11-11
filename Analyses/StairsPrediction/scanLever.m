function R = scanLever(base, fD, vt, spec, doFull)
%SCANLEVER  Sweep one parameter, scoring the staircase AND the paid-for objective.
%
%   R = scanLever(BASE, FD, VT, SPEC, DOFULL)
%
%   SPEC.name    parameter field to sweep
%   SPEC.vals    values to set it to
%   SPEC.comp    'none' | 'kstiff2'  — hold the isometric ACTIVE force constant
%                by rescaling kstiff2/kstiff2_n on the frozen isometric p2
%                distribution (only meaningful for 'estiff').
%   SPEC.ref     struct from stairsIsoRef() — required when comp = 'kstiff2'
%
%   Returns a struct array with the staircase LSQ cost, the block feature cost,
%   the isometric force actually reached, and the full objective + breakdown.

n = numel(spec.vals);
R = struct('val',num2cell(spec.vals(:)'), 'scale',[], 'Elsq',[], 'Cfeat',[], ...
           'parts',[], 'Fiso',[], 'Et',[], 'featCost',[], 'fM',[]);
tMov = vt(find(vt(:,2)~=0,1),1);

for i = 1:n
    p = base;
    p.(spec.name) = spec.vals(i);
    sc = 1;
    if isfield(spec,'comp') && strcmp(spec.comp,'kstiff2')
        r  = spec.ref;
        sc = (r.Fiso - r.Fp1) / (base.kstiff2 * r.pw(spec.vals(i)));
        p.kstiff2   = base.kstiff2   * sc;
        p.kstiff2_n = base.kstiff2_n * sc;
    end
    R(i).scale = sc;

    ps = p; ps.RunStairs = 1; ps.PlotEachSeparately = false;
    ps.ghostSave = ''; ps.ghostLoad = '';
    try
        [Elsq, os] = runStairsExperiment(ps);
        Fm = interp1(os.t, os.Force, tMov, 'linear','extrap');
        [Cf, parts, fM] = stairsCost(os.t, os.Force/Fm, vt, fD);
        R(i).Elsq=Elsq; R(i).Cfeat=Cf; R(i).parts=parts; R(i).fM=fM; R(i).Fiso=Fm;
    catch ME
        fprintf('   stairs FAILED %s=%g: %s\n', spec.name, spec.vals(i), ME.message);
        R(i).Elsq=NaN; R(i).Cfeat=NaN; R(i).Fiso=NaN;
    end

    if doFull
        pf = p; pf.RunStairs = 0;
        try
            [Et,~,fc] = evaluateBakersExp(pf.g, pf, true, false);
            R(i).Et = Et; R(i).featCost = fc;
        catch ME
            fprintf('   objective FAILED %s=%g: %s\n', spec.name, spec.vals(i), ME.message);
            R(i).Et = NaN;
        end
    end
    fprintf('  %s=%-9.4g scale=%8.3f Fiso=%6.1f | LSQ=%.5f feat=%7.3f | Et=%9.4f\n', ...
        spec.name, spec.vals(i), R(i).scale, R(i).Fiso, R(i).Elsq, R(i).Cfeat, R(i).Et);
end
end
