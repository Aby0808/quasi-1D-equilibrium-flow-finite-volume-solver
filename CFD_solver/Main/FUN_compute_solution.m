function [] = FUN_compute_solution()

% This is the driver function that computes and advances solution in time
% Updates ghost states
% Calculates RHS
% Advances solution by a time step
% Computes residuals
% Checks for convergence/stopping criteria

% Access globals
global COUNTER SOL

fprintf('\n Calculating solution....\n');
fprintf('\n iter\t flowtime\t continuity(L2)\t momentum(L2)\t energy(L2)');

% Advance solution

while true  % Main loop to advance solution in time
    COUNTER.iter = COUNTER.iter + 1;
    COUNTER.flowtime = COUNTER.flowtime + mean(SOL.dt);

    % Apply boundary conditions: Update ghost states
    FUN_apply_BC();

    % Compute RHS: Fluxes + Source terms
    FUN_compute_RHS();

    % Advance solution in time: Explicit for now
    FUN_time_intergrator_FE();

    % Compute and print residuals
    FUN_solution_residual();

    % Save solution if needed
    FUN_save_solution();

    % Check stopping criteria
    if FUN_check_stop_criteria()
        break;
    end

end

end