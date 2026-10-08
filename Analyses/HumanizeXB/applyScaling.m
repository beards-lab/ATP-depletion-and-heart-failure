function p = applyScaling(p, S)
% APPLYSCALING  Multiply named parameter fields of a loaded (mods-free) params
% struct by the factors in S (struct: field -> multiplier). Fields that are
% '=...' expressions are resolved first so the product is numeric. Call BEFORE
% getParams. Unknown fields raise an error so typos cannot pass silently.
fn = fieldnames(S);
for i = 1:numel(fn)
    f = fn{i};
    if ~isfield(p, f), error('applyScaling:unknown', 'no field %s', f); end
    if ischar(p.(f)) || isstring(p.(f))
        error('applyScaling:expr', '%s is an expression (%s); resolve it before scaling', f, p.(f));
    end
    p.(f) = p.(f) * S.(f);
end
end
