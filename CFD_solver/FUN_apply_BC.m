function Uext = FUN_apply_BC()
% Applies boundary conditions and returns extended state Uext (with ghosts).
% This is the ONLY function allowed to fill ghost cells.

global SOL MESH

U  = SOL.U;               % N x Nv (physical)
N  = MESH.N;
ng = MESH.ng;
Nv = SOL.Nv;

% Allocate extended
Uext = zeros(MESH.Ntot, Nv);
Uext(MESH.iC,:) = U;

% ---- Left boundary procedure
Uext = FUN_bc_left(Uext);

% ---- Right boundary procedure
Uext = FUN_bc_right(Uext);

end