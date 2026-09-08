function showBlocks(S, lbl, name, fn)
%SHOWBLOCKS  Print a lever sweep as block-cost + objective-breakdown tables.
ok = ~cellfun(@isempty, {S.parts});
P  = cat(1, S(ok).parts); v = [S(ok).val];
bl = fieldnames(P);
fprintf('\n== %s ==\n%-22s', lbl, name); fprintf('%9.4g', v); fprintf('\n');
for j = 1:numel(bl)
    fprintf('%-22s', bl{j}); fprintf('%9.3f', [P.(bl{j})]); fprintf('\n');
end
fprintf('%-22s','STAIRS feat total'); fprintf('%9.3f',[S(ok).Cfeat]); fprintf('\n');
fprintf('%-22s','stairs LSQ');        fprintf('%9.4f',[S(ok).Elsq]);  fprintf('\n');
if nargin >= 4 && ~isempty(fn) && ~all(cellfun(@isempty,{S.featCost}))
    FC = cat(1, S(ok).featCost);
    big = max(FC,[],1) > 0.05;
    fprintf('   -- objective terms above 0.05 --\n');
    for j = find(big)
        fprintf('%-22s', fn{j}); fprintf('%9.3f', FC(:,j)'); fprintf('\n');
    end
    fprintf('%-22s','OBJECTIVE Et'); fprintf('%9.3f',[S(ok).Et]); fprintf('\n');
end
end
