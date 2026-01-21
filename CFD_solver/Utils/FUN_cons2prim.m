function V = FUN_cons2prim(U)
% Convert conservative -> primitive for 1D perfect-gas Euler.
% U: [N x 3] or [1 x 3] = [rho, rho*u, rho*E]
% V: [N x 3]            = [rho, u, T]

global CONST
gamma = CONST.gamma;
R     = CONST.Rair;
Cv    = R/(gamma-1);

rho = U(1,:);
m   = U(2,:);
Et  = U(3,:);

u = m ./ rho;

E = Et ./ rho;
e = E - 0.5.*u.^2;

T = e ./ Cv;

V = zeros(size(U));
V(1,:) = rho;
V(2,:) = u;
V(3,:) = T;
end