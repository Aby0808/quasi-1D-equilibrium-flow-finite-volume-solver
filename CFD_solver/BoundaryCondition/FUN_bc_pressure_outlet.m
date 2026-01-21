function Ug = FUN_bc_pressure_outlet(Us, Ug, bc)

% Boundary condition procedure to impose pressure outlet boundary condition

% find primitive variables
% assuming perfect gas
global CONST

gama = CONST.gamma;
R = CONST.Rair;

% get interior states
Us_prim = FUN_cons2prim(Us(:,1));
rho = Us_prim(1);
u = Us_prim(2);
T = Us_prim(3);
if T<0
    error('Negative Temperature detected at outlet boundary: T = %f',T);
end

M = sqrt(u^2)/sqrt(gama*R*T);

% check for Mach number

if M<1  % if subsonic case: use P at boundary
    if ~isfield(bc,'P') || isempty(bc.P)
        error('pressure_outlet requires bc.P for subsonic outflow.');
    end
    P = bc.P; % get boundary value
    rho = P/(R*T);
end

Up = [rho, u, T]; % primitive variables to be imposed
Ub = FUN_prim2cons(Up);  % convert to conservative variables

% find and set ghost state

ng = size(Ug,2);
Ug(:,:) = repmat(Ub, ng, 1);

end