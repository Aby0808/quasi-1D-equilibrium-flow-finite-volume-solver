function stop = FUN_check_stop_criteria()
% This function returns true if solver should stop (hardcoded for now)

global SOL COUNTER

stop = false;

maxIter = 5000;
tolLinf = 1e-6;

% Blow-up guard
if any(~isfinite(SOL.U(:)))
    warning('NaN/Inf detected in SOL.U at iter=%d. Stopping.', COUNTER.iter);
    stop = true;
    return;
end

% Convergence (needs FUN_solution_residual() to store SOL.res.Linf)
if isfield(SOL,'res') && isfield(SOL.res,'Linf') && numel(SOL.res.Linf) >= 3
    resMax = max(SOL.res.Linf(1:3));  % cont, mom, energy
    if resMax < tolLinf
        fprintf('Converged at iter=%d (max Linf = %.3e < %.3e)\n', COUNTER.iter, resMax, tolLinf);
        stop = true;
        return;
    end
end

% Max iterations
if COUNTER.iter >= maxIter
    fprintf('Reached maxIter=%d. Stopping.\n', maxIter);
    stop = true;
    return;
end

end
