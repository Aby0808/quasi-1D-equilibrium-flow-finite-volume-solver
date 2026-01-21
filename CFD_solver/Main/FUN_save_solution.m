function [] = FUN_save_solution()
% This function saves solution periodically (hardcoded for now)

% Access globals
global SOL MESH COUNTER

saveEvery = 200;
outFolder = fullfile(pwd,'Results');

if COUNTER.iter == 1
    if ~exist(outFolder,'dir')
        mkdir(outFolder);
    end
end

if mod(COUNTER.iter, saveEvery) ~= 0
    return;
end

tag = sprintf('iter_%06d', COUNTER.iter);
fname = fullfile(outFolder, ['sol_' tag '.mat']);

save(fname, 'SOL', 'MESH', 'COUNTER');

% Rolling latest file
save(fullfile(outFolder, 'sol_latest.mat'), 'SOL', 'MESH', 'COUNTER');

end
