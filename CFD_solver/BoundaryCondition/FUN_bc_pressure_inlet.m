function Ug = FUN_bc_pressure_inlet(Us, Ug, bc)

% Boundary condition procedure to impose pressure inlet boundary condition

% extract boundary values
P0 = bc.P0;
P = bc.P;
T0 = bc.T0;
dir = bc.dir;

% find primitive variables
% assuming perfect gas
global CONST

gama = CONST.gamma;
R = CONST.Rair;

% find quantities at boundary
M = sqrt((((P0/P)^((gama-1)/gama)) - 1) * (2/(gama-1)));
T = T0/(1 + 0.5*(gama-1)*M^2);
a = sqrt(gama*R*T);
u = M*a;

% check for Mach number

if M<1  % if subsonic case: use interior P
    Us_prim = FUN_cons2prim(Us(:,1));
    P = Us_prim(1)*R*Us_prim(3);  % extrapolate pressure from interior
    
    % compute properties according to extrapolated P
    M = sqrt((((P0/P)^((gama-1)/gama)) - 1) * (2/(gama-1)));
    T = T0/(1 + 0.5*(gama-1)*M^2);
    a = sqrt(gama*R*T);
    u = M*a;
end

u = u*dir(1); % set velocity direction
rho = P/(R*T);

Up = [rho, u, T]; % primitive variables to be imposed
Ub = FUN_prim2cons(Up);  % convert to conservative variables

% find and set ghost state

ng = size(Ug,2);
Ug(:,:) = repmat(Ub, ng, 1);

end