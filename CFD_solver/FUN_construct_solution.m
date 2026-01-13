function sol = FUN_construct_solution()
% Creates flowfield/solution structure with consistent sizes.
% Default problem: quasi-1D Euler with U = [rho, m, Et].

global MESH
N = MESH.N;

sol.Nv = 3;
sol.U  = zeros(N, sol.Nv);     % physical cells only (N x 3)

% Optional: names for readability
sol.varNames = {'rho','m','Et'};

% Optional: placeholder for time, iteration count, etc.
sol.it = 0;
sol.t  = 0.0;

end
