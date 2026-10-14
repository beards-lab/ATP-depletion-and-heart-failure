function [base, th] = buildModel(cfg, x)
%BUILDMODEL  Params struct and thin-filament/Ca options for parameter vector X under tier CFG.
%   CFG.base is the resolved seed (numeric fields). CFG.fixedScale (struct field->factor) is
%   applied first, then every free parameter CFG.par(i):
%     type 'scale'  multiply CFG.par(i).fields by 10^x (log) or x
%     type 'abs'    set field CFG.par(i).fields{1} to 10^x or x
%     type 'tf'     thin filament: th.(name) = 10^x or x
%     type 'shift'  th.shift = x (s)
%   th also carries CFG.ca (CaKind, CaOpt), CFG.pre, CFG.nb, CFG.CaModel.
base = cfg.base;
if isfield(cfg, 'fixedScale')
    for f = string(fieldnames(cfg.fixedScale))', base.(f) = base.(f)*cfg.fixedScale.(f); end
end
th = struct('CaKind', cfg.ca.kind, 'CaOpt', cfg.ca.opt, 'pre', cfg.pre, 'nb', cfg.nb);
if isfield(cfg, 'CaModel'), th.CaModel = cfg.CaModel; end
if isfield(cfg, 'CaInput'),  th.CaInput = cfg.CaInput; end      % 'tworld': T-World as the Ca source
if isfield(cfg, 'CaTworld'), th.CaTworld = cfg.CaTworld; end    % its options (coupling, prebeats, param, ...)
if isfield(cfg, 'tfFixed')
    for f = string(fieldnames(cfg.tfFixed))', th.(f) = cfg.tfFixed.(f); end
end
for i = 1:numel(cfg.par)
    p = cfg.par(i);  v = x(i);  if p.log, v = 10^v; end
    switch p.type
        case 'scale', for f = string(p.fields), base.(f) = base.(f)*v; end
        case 'abs',   base.(p.fields{1}) = v;
        case 'tf',    th.(p.name) = v;
        case 'tw',    th.tw.(p.name) = v;     % T-World parameter (getDefaultParameters.m name), absolute
        case 'shift', th.shift = v;
    end
end
end
