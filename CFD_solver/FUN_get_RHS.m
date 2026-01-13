function [A, RHS] = FUN_get_RHS(Uext)
% FUN_GET_RHS  Assemble the semi-discrete residual (RHS) for quasi-1D FV Euler.
%
% Governing semi-discrete form (cell-centered conservative update):
%   dU_i/dt = RHS_i = -(1/dx_i) * ( F_{i+1/2} - F_{i-1/2} ) + S_i
%
% Inputs:
%   Uext : Extended conservative state array including ghost cells.
%          Size = (N + 2*ng) x Nv
%          - Ghost cells must already be filled by the BC layer.
%          - Physical cells are embedded inside Uext at indices MESH.iC.
%
% Outputs:
%   A    : Placeholder Jacobian (empty for explicit time stepping; fill later for implicit).
%   RHS  : Cell-centered RHS for physical cells only. Size = N x Nv.
%
% Notes / Design rules:
%   - SOL.U is the authoritative solution (physical cells only) and is NOT modified here.
%   - Uext is a derived workspace used for reconstruction/flux evaluation near boundaries.
%   - Boundary handling is done ONLY by the BC procedure that created Uext.
%   - Reconstruction MUST NOT special-case i=1 or i=N; ghost cells cover the stencil.
%
% Typical usage in main loop (later):
%   Uext = FUN_apply_BC();         % build/fill ghost cells from SOL.U and BC
%   [~, RHS] = FUN_get_RHS(Uext);  % compute residual
%   SOL.U = FUN_time_integrator(SOL.U, RHS, dt);  % update elsewhere (explicit/implicit)

global MESH SOL

% ----------------------------
% 1) Basic dimensions and geometry
% ----------------------------
N  = MESH.N;      % number of physical cells
Nv = SOL.Nv;      % number of conservative variables per cell (e.g. 3)

dx = MESH.dx;     % cell widths (N x 1 recommended)
if isscalar(dx)
    dx = dx * ones(N,1);
end

% ----------------------------
% 2) Allocate outputs
% ----------------------------
RHS = zeros(N, Nv);   % physical cells only
A   = [];             % empty for explicit (implicit later)

% ----------------------------
% 3) Loop over physical cells and assemble RHS
% ----------------------------
for i = 1:N

    % ---------------------------------------------------------
    % 3a) Reconstruct interface states (piecewise-constant for now)
    %     Returns Nv x 1 vectors:
    %       (UL,UR) at i-1/2 and i+1/2
    % ---------------------------------------------------------
    [UL_imh, UR_imh, UL_iph, UR_iph] = FUN_reconstruct_interfaces(Uext, i);

    % ---------------------------------------------------------
    % 3b) Numerical fluxes at interfaces (Roe LTE)
    %     Flux vectors are Nv x 1
    % ---------------------------------------------------------
    [F_imh, ~] = FUN_flux_num_roe_lte(UL_imh, UR_imh);  % F_{i-1/2}
    [F_iph, ~] = FUN_flux_num_roe_lte(UL_iph, UR_iph);  % F_{i+1/2}

    % ---------------------------------------------------------
    % 3c) Source term (cell-centered)
    %     Must return Nv x 1 vector. For pure Euler w/o geometry,
    %     this can be zeros(Nv,1).
    % ---------------------------------------------------------
    S_i = 0;%FUN_get_source_terms(i);  % Nv x 1

    % ---------------------------------------------------------
    % 3d) Assemble RHS for cell i
    % ---------------------------------------------------------
    rhs_i = -(F_iph - F_imh) / dx(i) + S_i;  % Nv x 1
    RHS(i,:) = rhs_i.';                     % store row-wise (1 x Nv)
end

end
