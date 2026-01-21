function U = FUN_prim2cons(V)
% Convert primitive -> conservative for 1D perfect-gas Euler.
% V: [N x 3] or [1 x 3] = [rho, u, T]
% U: [N x 3]            = [rho, rho*u, rho*E]

global CONST
gamma = CONST.gamma;
R     = CONST.Rair;
Cv    = R/(gamma-1);

rho = V(:,1);
u   = V(:,2);
T   = V(:,3);

E = Cv.*T + 0.5.*u.^2;

U = zeros(size(V));
U(:,1) = rho;
U(:,2) = rho .* u;
U(:,3) = rho .* E;
end