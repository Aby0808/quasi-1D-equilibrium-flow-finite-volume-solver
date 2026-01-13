function [A, RHS] = FUN_construct_lin_sys()
%FUN_construct_lin_sys
% Iterates over the grid and constructs the linear system matrix A and RHS.
%
% Philosophy:
%   - This function only orchestrates assembly.
%   - Fluxes, source terms, BCs, indexing, etc. live in their own functions.

% ----------------------------
% 1) Access grid and solution
% ----------------------------
% These can be globals, persistent structs, or returned by getters.
% Keep it abstract: just "get what you need".
[grid, sol, opts] = FUN_get_context();   % you will implement this later

N  = grid.N;          % number of cells
Nv = sol.Nv;          % number of variables per cell
Ndof = N * Nv;

% ----------------------------
% 2) Allocate system
% ----------------------------
% Start simple: dense RHS, sparse A.
% (Later you can switch to triplets prealloc for speed.)
A   = sparse(Ndof, Ndof);
RHS = zeros(Ndof, 1);

% Local DOF mapper (cell i, var k) -> global index
gid = @(i,k) (i-1)*Nv + k;

% ----------------------------
% 3) Iterate over the grid
% ----------------------------
for i = 1:N

    % ---------------------------------------------------------
    % 3a) Gather stencil state needed for flux reconstruction
    % ---------------------------------------------------------
    stencil = FUN_get_stencil(sol, i);   % returns Ui, Uim1, Uip1, etc.

    % ---------------------------------------------------------
    % 3b) Inviscid flux contributions (residual + jacobian block)
    % ---------------------------------------------------------
    % Contract: return residual contribution for cell i (Nv x 1)
    % and Jacobian blocks wrt neighbor states (Nv x Nv each).
    [Rinv, Jm, J0, Jp] = FUN_get_inv_fluxes(grid, stencil, i, opts);

    % ---------------------------------------------------------
    % 3c) Source term contributions (residual + jacobian block)
    % ---------------------------------------------------------
    [Rsrc, Jsrc] = FUN_get_source_terms(grid, sol, i, opts);

    % ---------------------------------------------------------
    % 3d) Assemble RHS
    % ---------------------------------------------------------
    % Convention: A*dU = -R  (Newton / implicit update form)
    Ri = Rinv + Rsrc;                 % Nv x 1
    for k = 1:Nv
        RHS(gid(i,k)) = RHS(gid(i,k)) - Ri(k);
    end

    % ---------------------------------------------------------
    % 3e) Assemble matrix blocks (only if implicit)
    % ---------------------------------------------------------
    if opts.implicit
        % Center block
        if ~isempty(J0)
            A = FUN_add_block(A, gid, i, i, J0);
        end

        % Left neighbor block
        if i > 1 && ~isempty(Jm)
            A = FUN_add_block(A, gid, i, i-1, Jm);
        end

        % Right neighbor block
        if i < N && ~isempty(Jp)
            A = FUN_add_block(A, gid, i, i+1, Jp);
        end

        % Source Jacobian (diagonal)
        if ~isempty(Jsrc)
            A = FUN_add_block(A, gid, i, i, Jsrc);
        end
    end

end

% ----------------------------
% 4) Apply boundary conditions
% ----------------------------
% Keep BC handling out of the loop for cleanliness.
[A, RHS] = FUN_apply_BC(A, RHS, grid, sol, opts);

end
