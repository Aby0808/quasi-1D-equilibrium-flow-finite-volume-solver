function [] = FUN_set_initial_solution()
% Initializes SOL.U from primitive variables (P, T, u)
% Converts to conservative variables.
% Physical cells only. No ghost cells.

global SOL MESH CONST BC

N = MESH.N;

% --- primitive IC (use inlet values)
P = BC.left.P;     % Pa
T0 = BC.left.T0;     % K
u = 1.0;           % m/s (small value to avoid surprises)
T = T0;             % assuming reservoir conditions for now

% --- ideal gas relations
rho = P / (CONST.Rair * T);

% total specific energy
E = CONST.Cv * T + 0.5 * u^2;

% --- fill conservative variables
SOL.U(1,:) = rho;        % density
SOL.U(2,:) = rho * u;    % momentum
SOL.U(3,:) = rho * E;    % total energy

end