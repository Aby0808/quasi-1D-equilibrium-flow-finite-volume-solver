function [] = FUN_initialize_time_integrator_FE()
% This function sets a global timestep size according to the user or local time steps according to CFL

% Access globals
global SOL MESH

CFL = 0.5;
dt = 1e-7;  % placeholder

SOL.dt = ones(1,MESH.N)*dt;

end