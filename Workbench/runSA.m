function [E, E0] = runSA(params0, fieldNames, delta, drawPlots)
E0 = evaluateBakersExp(params0.g, params0);

if nargin < 2
    fieldNames = fieldnames(params0);
end

if nargin < 3 
    delta = 0.1; % 10% difference
end

if nargin < 4 
    drawPlots = true;
end

% compare with defaults
defaults = getParams();

for k = 1:98 %length(fieldNames)
    
    params = params0;
    if isfield(defaults, fieldNames{k}) && islogical(defaults.(fieldNames{k}))
        
        % just switch
        params.(fieldNames{k}) = ~params.(fieldNames{k});
        fprintf("%d: flipping %s...", k,  fieldNames{k});

        tic
        errval = evaluateBakersExp(params.g, params);
        t = toc;
        E(k, 4) = t;
        
        if errval > 1e6
            E(k, 1) = NaN;
        else
            E(k, 1) = errval;
        end

        if errval == E0
            fprintf(">>> no effect \n", errval/E0*1e2 - 100);
        else
            fprintf(">>> %+0.3f%% \n", errval/E0*1e2 - 100);
        end

    elseif iscell(params.(fieldNames{k}))
        continue;
    else %if length(params.(fieldNames{k})) == 1
        fprintf("%d: %s %+0.3f%%...", k,  fieldNames{k}, delta(1));
        params.(fieldNames{k}) = params.(fieldNames{k}).*(1 + delta(1));

        tic
        errval = evaluateBakersExp(params.g, params);
        t = toc;
        E(k, 4) = t;
        
        if errval > 1e6
            E(k, 1) = NaN;
        else
            E(k, 1) = errval;
        end

        if errval == E0
            fprintf(">>> no effect \n", errval/E0*1e2 - 100);
        else
            fprintf(">>> %+0.3f%% \n", errval/E0*1e2 - 100);
        end
    
        if length(delta) > 1
            fprintf("%d: %s %+0.3f%% ...", k,  fieldNames{k}, delta(2));
            params.(fieldNames{k}) = params.(fieldNames{k})*(1 + delta(2));
            tic
            errval = evaluateBakersExp(params.g, params);
            t = toc;
            E(k, 4) = (E(k, 4) + t)/2;
            
            if errval > 1e6
                E(k, 2) = NaN;
            else
                E(k, 2) = errval;
            end
    
            if errval == E0
                fprintf(">>> no effect \n", errval/E0*1e2 - 100);
            else
                fprintf(">>> %+0.3f%% \n", errval/E0*1e2 - 100);
            end
        end
    end



    % does this even matter? 
    % i.e. does using zero produces error lower than double?
%     if tryZeros
%         g_s(k) = 0;
%         try
%             E_g0(k) = evaluateProblem(fcn, g_s, false);
%         catch
%             E_g0(k) = NaN;
%         end
%     end

end
disp("Sensitivities Done!")
% else
%     E = sens;
%     disp('sensitivities provided, plotting..')
% end

% if ~drawPlots
%     return;
% end
%%
save('sens_02.mat', 'E');
figure(52);clf; hold on;
% subplot(211);cla;hold on;

E(E == 0) = NaN;
plotInds = ~isnan(E(:, 1));
plotE = E(plotInds, :);
bar(plotE(:, 1:2)/E0*100 - 100);
% plot([1, length(E)], E0*[1 1], '--r')
for i = 1:length(fieldNames)
    g_lab{i} = fieldNames{i} + " (" + num2str(i) + ")";
end
xticklabels(g_lab(plotInds));
xtickangle(45);
%% 

% subplot(212);cla;hold on;
% title('Relative sensitivity to one-at-a-time perturbation of G by \delta');
% bar(E/E0*100 - 100);
% plot([0, 16], [1 1], '--r')
% 
% ylabel('E / E0');
% ylim([0, 10])
% xticks(1:20)
% xticklabels(g_names);
% xtickangle(45);
% 
% is the param zeroable?
% if tryZeros
% for k = 1:length(g)
%     if ~isnan(g(k)) && E_g0(k) < 2*E0
%         text(k, E_g0(k)/E0,'*', 'FontSize', 18);
%     end
% end
% end
% plot(nan, nan, '*k');
% legend('g(x) - \delta', 'baseline', 'g(x) + \delta', 'baseline', 'g(x) = 0');
% ylim([0.95 1.05])