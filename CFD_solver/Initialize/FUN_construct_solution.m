function sol = FUN_construct_solution()
% Creates flowfield/solution structure with consistent sizes.
% Default problem: quasi-1D Euler with U = [rho, m, Et].

global MESH
N = MESH.N;

sol.Nv = 3;
sol.U  = zeros(sol.Nv, N);     % physical cells only (N x 3)
sol.RHS  = zeros(sol.Nv, N);
sol.res  = zeros(sol.Nv, N);

% Names for readability
sol.varNames = {'rho','m','Et'};

% Placeholder for time, iteration count, etc.
sol.it = 0;
sol.t  = zeros(1,N);
sol.dt = zeros(1,N);
sol.cfl = 0.0;

end
