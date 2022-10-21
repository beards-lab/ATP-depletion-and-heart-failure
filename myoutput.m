function stop = myoutput(gr,optimvalues,state)
    stop = false;    
    if ~isequal(state,'iter') || mod(optimvalues.iteration, 10) > 0
        return;
    end
%     g_all = [gr(1) 3 gr(2) 0.8 gr(3:end)];
    evaluateProblem(@dPUdTCa, gr, true, [1 1 1 1 0 1]);
%     evaluateProblem(@dPUdT, [gr(1) 0.0555 gr(2) 3.6484 gr(3:end)], true);
end
