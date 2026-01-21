function FUN_solution_residual()
% This function prints residual norms (explicit solver): continuity, momentum, energy

global SOL COUNTER

rhs = SOL.RHS;
if isempty(rhs)
    error('FUN_solution_residual: SOL.RHS is empty.');
end

% Assume Nv = 3 -> [rho, rho*u, rho*E]
if size(rhs,2) < 3
    error('FUN_solution_residual: expected SOL.RHS to have at least 3 columns.');
end

it = 0;
if ~isempty(whos('global','COUNTER')) && isstruct(COUNTER) && isfield(COUNTER,'iter')
    it = COUNTER.iter;
    ft = COUNTER.flowtime;
end

N = size(rhs,1);

% Linf norms
Linf = max(abs(rhs), [], 1);              % 1x3
% L2 norms (RMS)
L2   = sqrt(sum(rhs.^2, 1) / max(N,1));   % 1x3


% SOL.res.Linf = Linf(1:3);
% SOL.res.L2   = L2(1:3);

fprintf('\n %d\t %d\t %d\t %d\t %d',it, ft, L2(1), L2(2), L2(3));
% fprintf('it=%6d | Linf: cont=%.3e mom=%.3e ener=%.3e | L2: cont=%.3e mom=%.3e ener=%.3e\n', ...
%     it, Linf(1), Linf(2), Linf(3), L2(1), L2(2), L2(3));

end
