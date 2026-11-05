function L = drugLevers()
%DRUGLEVERS  Single-mechanism levers of the human-EHT twitch model, each with a dose grid.
%   L(i).name, L(i).label, L(i).grid (factor values; 1 = untreated), L(i).apply(base, th, f) ->
%   [base, th], L(i).drugs (which compounds of the ARPA panel the lever stands in for).
%   Ca transient levers act on th.CaOpt (risedecay: Cadia, Camax, tau_r, t50, hdec); thin-filament
%   levers on th.tf_*; cross-bridge levers on BASE fields. ATP is set in mM (grid = mM, not a factor).
k = 0;
k = k+1; L(k) = mk('CaAmp',   'Ca amplitude x f (Camax-Cadia)',            [0.5 0.6 0.7 0.8 0.9 1.15 1.3 1.5], @(b,t,f) caAmp(b,t,f), 'diltiazem, clozapine (ICaL block); isoproterenol (+); doxorubicin/carfilzomib (SR load)');
k = k+1; L(k) = mk('CaDecay', 'Ca decay t50 x f (SERCA)',                  [0.7 0.85 1.2 1.5 2],   @(b,t,f) caOpt(b,t,'t50',f), 'isoproterenol (PLB, f<1); SERCA impairment (f>1): doxorubicin, carfilzomib (PP2A)');
k = k+1; L(k) = mk('CaRise',  'Ca rise tau_r x f',                          [0.7 1.5 2],           @(b,t,f) caOpt(b,t,'tau_r',f), 'ICaL kinetics; minor');
k = k+1; L(k) = mk('ka',      'attachment ka x f',                         [0.1 0.2 0.35 0.5 0.7], @(b,t,f) sc(b,t,{'ka'},f), 'mavacamten (weak-binding trap), blebbistatin');
k = k+1; L(k) = mk('k1',      'Pi release / power stroke k1,k_1 x f',      [0.1 0.2 0.35 0.5 0.7], @(b,t,f) sc(b,t,{'k1','k_1'},f), 'mavacamten, aficamten (canonical target: Pi release)');
k = k+1; L(k) = mk('SRXstab', 'SRX stabilisation: ksr0,ksrd x f, kmsrd / f', [1.5 2 3 5 10],       @(b,t,f) srx(b,t,f), 'mavacamten (IHM/SRX); PKA on MyBP-C is the inverse (f<1)');
k = k+1; L(k) = mk('SRXrel',  'SRX release: kmsrd x f, ksr0,ksrd / f',    [1.5 2 3],             @(b,t,f) srx(b,t,1/f), 'isoproterenol (MyBP-C phosphorylation), danicamtiv');
k = k+1; L(k) = mk('k2',      'ADP release k2,k_2 x f',                    [0.35 0.5 0.7 1.5],    @(b,t,f) sc(b,t,{'k2','k_2'},f), 'omecamtiv-like (f<1); proxy for ATP-limited detachment');
k = k+1; L(k) = mk('MgATP',   'MgATP (mM), both ATP gates on',             [4 2 1 0.5 0.25],      @(b,t,f) setf(b,t,'MgATP',f), 'doxorubicin, sunitinib, carfilzomib (mitochondrial ATP depletion)');
k = k+1; L(k) = mk('tf_koff', 'thin filament k_off x f',                   [0.7 1.5 2 3],         @(b,t,f) tf(b,t,'tf_k_off',f), 'isoproterenol via cTnI-PKA (absent with ssTnI in hiPSC-CM)');
k = k+1; L(k) = mk('tf_kon',  'thin filament k_on x f',                    [0.7 0.85 1.2 1.5],    @(b,t,f) tf(b,t,'tf_k_on',f), 'Ca sensitisers; ssTnI; pH/Pi desensitisation (f<1)');
k = k+1; L(k) = mk('kd',      'weak-binding detachment kd x f',            [2 4],                 @(b,t,f) sc(b,t,{'kd'},f), 'aficamten-like');
k = k+1; L(k) = mk('kstiff',  'cross-bridge stiffness kstiff1,2 x f',      [0.5 0.7],             @(b,t,f) sc(b,t,{'kstiff1','kstiff2'},f), 'sarcomere disarray (doxorubicin, carfilzomib): per-head force');
end

function s = mk(name, label, grid, apply, drugs)
s = struct('name', name, 'label', label, 'grid', grid, 'apply', apply, 'drugs', drugs);
end
function [b, t] = caAmp(b, t, f)
o = t.CaOpt;  t.CaOpt.Camax = o.Cadia + f*(o.Camax - o.Cadia);
end
function [b, t] = caOpt(b, t, fld, f)
t.CaOpt.(fld) = f*t.CaOpt.(fld);
end
function [b, t] = sc(b, t, flds, f)
for i = 1:numel(flds), b.(flds{i}) = f*b.(flds{i}); end
end
function [b, t] = srx(b, t, f)
b.ksr0 = f*b.ksr0;  b.ksrd = f*b.ksrd;  b.kmsrd = b.kmsrd/f;
end
function [b, t] = setf(b, t, fld, v)
b.(fld) = v;
end
function [b, t] = tf(b, t, fld, f)
t.(fld) = f*t.(fld);
end
