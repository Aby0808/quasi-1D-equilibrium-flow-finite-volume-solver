function FUN_compute_RHS()
% FUN_GET_RHS  Assemble the semi-discrete residual (RHS) for quasi-1D FV Euler.
%
% Governing semi-discrete form (cell-centered conservative update):
%   dU_i/dt = RHS_i = -(1/dx_i) * ( F_{i+1/2} - F_{i-1/2} ) + S_i

% For explicit schemes only

global MESH SOL COUNTER

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
RHS = zeros(Nv, N);   % physical cells only
% A   = [];             % empty for explicit (implicit later)

% ----------------------------
% 3) Loop over physical cells and assemble RHS
% ----------------------------
for i = 1:N

    % Update global counter
    COUNTER.pos = i;

    % Numerical fluxes at interfaces
    [F_imh, F_iph] = FUN_get_inv_flux();
    % [F_imh, ~] = FUN_flux_num_roe_lte(UL_imh, UR_imh);  % F_{i-1/2}
    % [F_iph, ~] = FUN_flux_num_roe_lte(UL_iph, UR_iph);  % F_{i+1/2}

    % Source term (cell-centered)
    S_i = FUN_get_source_terms();  % Nv x 1

    % Assemble RHS for cell i
    rhs_i = -(F_iph - F_imh) / dx(i) + S_i;  % Nv x 1
    RHS(:,i) = rhs_i;                     % store row-wise (1 x Nv)
end

SOL.RHS = RHS;

end
