function Uext = FUN_construct_Uext(sol)
% Builds extended state array with ghost cells from sol.U (physical).
% Default BC: transmissive (zero-gradient) on both ends.

global MESH
U  = sol.U;          % N x Nv
N  = MESH.N;
ng = MESH.ng;
Nv = sol.Nv;

Uext = zeros(MESH.Ntot, Nv);

% fill physical region
Uext(MESH.iC,:) = U;

% left ghosts
for g = 1:ng
    Uext(ng+1-g,:) = U(1,:);
end

% right ghosts
for g = 1:ng
    Uext(ng+N+g,:) = U(N,:);
end

end