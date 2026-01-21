function FUN_setup_path()
%SETUP_PATH 
% Add all project subfolders to MATLAB path.

    root = fileparts(mfilename('fullpath'));
    addpath(genpath(root));
end
