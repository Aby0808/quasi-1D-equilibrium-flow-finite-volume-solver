function [] = FUN_initialize()
% This routine performs all required constructors to bring the solver into a
% consistent “ready-to-run” state. After this function returns, the mesh,
% thermochemistry interpolants, boundary-condition parameters, and the
% initial flowfield (in conservative variables) are available globally.
%
% Sequence:
%   1) Initialize constants (gas constants, calorically-perfect fallback, etc.)
%   2) Load equilibrium tables from disk
%   3) Build equilibrium-property interpolants (persistent cache inside FUN_get_eq_prop)
%   4) Generate the 1D finite-volume mesh (including ghost-indexing contract)
%   5) Allocate the solution/flowfield storage (physical cells only)
%   6) Define boundary-condition parameters (BC values; no enforcement here)
%   7) Initialize the flowfield using primitive ICs (P, T, u) and convert to conservative U
%
% Notes:
%   - SOL.U stores only physical cells (size N x Nv). Ghost states are NOT
%     stored here and will be created on demand by BC procedures later.
%   - Boundary-condition parameters are specified here (BC struct), but
%     ghost-filling / enforcement is handled separately during RHS evaluation.
%
% Globals created/updated:
%   CONST  : physical constants and fallback thermodynamic parameters
%   TABLES : equilibrium tables loaded from file
%   MESH   : FV mesh geometry and ghost indexing (ng, iC, phys2ext, etc.)
%   SOL    : solution struct (SOL.U, SOL.Nv, etc.)
%   BC     : boundary-condition parameters (types + values)

% --- clean slate (optional but recommended during development)
clear global CONST MESH TABLES SOL BC

fprintf('\n Initializing constants... ')

global CONST
CONST = struct();

% Universal gas constant (use correct value!)
CONST.Ru = 8.314462618;          % J/(mol·K)

% Air properties (for calorically perfect fallback)
CONST.Mair = 28.965e-3;          % kg/mol (mean molar mass of air)
CONST.Rair = CONST.Ru / CONST.Mair; % J/(kg·K)

CONST.gamma = 1.4;
CONST.Cp    = CONST.gamma*CONST.Rair/(CONST.gamma-1); % J/(kg·K)
CONST.Cv    = CONST.Cp - CONST.Rair;

fprintf('done! \n')

fprintf('\n Initializing equilibrium tables... ')
% Load the equilibrium tables
S = load('tables.mat');
global TABLES
TABLES = S.tables;
fprintf('done! \n')
fprintf('\n Creating interpolants for equilibrium tables... ')
% Create interpolants
FUN_get_eq_prop('initialize', TABLES);
fprintf('done! \n')

fprintf('\n Generating grid... ')
% Generate mesh
global MESH
MESH = FUN_generate_mesh(0, 1, 0.35, 100, 1); % 1 ghost cell for now

global GEOM
GEOM = FUN_load_area_profile('Area.dat');

% Construct solution matrix
global SOL
SOL = FUN_construct_solution();

fprintf('done! \n')

fprintf('\n Initializing solution... ')
% Set boundary values
global BC
BC = FUN_set_boundary_values();

% Initialize ghost states
global GHOST
GHOST = FUN_initialize_ghost_states(SOL.Nv, 1, BC); %1 ghost cell for now

% Initialize solution
FUN_set_initial_solution();

% set counters
global COUNTER
COUNTER = FUN_set_global_count();

fprintf('done! \n')

fprintf('\n Initializing time integrator... ')
% Set boundary values
FUN_initialize_time_integrator_FE();
fprintf('done! \n')

end
