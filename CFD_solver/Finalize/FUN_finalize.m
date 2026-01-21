function [] = FUN_finalize()
% FUN_FINALIZE  Clean up globals, persistent data, and figures.
% Safe to call after successful runs OR after errors/partial initialization.

fprintf('\n Finalizing solver and cleaning up... ')

% ---- finalize equilibrium interpolants (persistent storage)
% Make this safe even if initialize was not called or was interrupted.
try
    FUN_get_eq_prop('finalize');
catch
    % ok: interpolants may not exist yet
end

% ---- clear all globals created by FUN_initialize
clear global CONST MESH TABLES SOL BC GHOST GEOM COUNTER

% ---- close figures (optional)
close all

fprintf('done! \n')

end