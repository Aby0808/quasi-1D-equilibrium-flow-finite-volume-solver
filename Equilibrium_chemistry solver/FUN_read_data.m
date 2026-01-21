function [Species_Data, Species_Mapping, Atoms_Mapping] = FUN_read_data(file)
% FUN_read_data Reads thermo data and caches results to avoid re-reading.
%
% Cache behavior:
% - If called again with the same file and unchanged timestamp, returns cached data.
% - If file was modified, re-reads and refreshes cache.

persistent cache

if nargin < 1 || isempty(file)
    error('FUN_read_data:MissingFile', 'You must provide a file path.');
end

% Get file info for change detection
info = dir(file);
if isempty(info)
    error('FUN_read_data:FileNotFound', 'File not found: %s', file);
end

% Build cache key
key = string(which(file));
if key == ""   % if not on path, keep as given
    key = string(file);
end

% Return cached if valid
if ~isempty(cache) && isfield(cache, matlab.lang.makeValidName(key))
    entry = cache.(matlab.lang.makeValidName(key));
    if isequal(entry.datenum, info.datenum)
        Species_Data    = entry.Species_Data;
        Species_Mapping = entry.Species_Mapping;
        Atoms_Mapping   = entry.Atoms_Mapping;
        return;
    end
end

%% --------- (Your original parsing logic) ----------
fileID = fopen(file, 'r');
if fileID == -1
    error('FUN_read_data:OpenFail', 'Cannot open file: %s', file);
end

Species_Data = struct();
Species_Mapping = containers.Map();
Atoms_Mapping = containers.Map();
currentSpecies = '';

while ~feof(fileID)
    line = strtrim(fgetl(fileID));
    if isempty(line); continue; end

    if startsWith(line, 'Species')
        tokens = strsplit(line);
        currentSpecies = tokens{2};

        validFieldName = matlab.lang.makeValidName(currentSpecies);
        Species_Mapping(currentSpecies) = validFieldName;

        Species_Data.(validFieldName) = {};
        Atoms_Mapping(currentSpecies) = zeros(1, 3); % [H O N] example

    elseif contains(line, '=') && ~isempty(currentSpecies)
        tokens = regexp(line, '(\w+)=(\d+)', 'tokens');
        atom_counts = zeros(1, 3);

        for k = 1:length(tokens)
            element = tokens{k}{1};
            count = str2double(tokens{k}{2});
            switch element
                case 'H', atom_counts(1) = count;
                case 'O', atom_counts(2) = count;
                case 'N', atom_counts(3) = count;
            end
        end
        Atoms_Mapping(currentSpecies) = atom_counts;

    elseif ~isempty(currentSpecies)
        numValues = sscanf(line, '%e');
        if ~isempty(numValues)
            validFieldName = Species_Mapping(currentSpecies);
            Species_Data.(validFieldName){end + 1} = numValues';
        end
    end
end

fclose(fileID);

%% --------- Store into cache ----------
cacheField = matlab.lang.makeValidName(key);
cache.(cacheField) = struct( ...
    'datenum', info.datenum, ...
    'Species_Data', Species_Data, ...
    'Species_Mapping', Species_Mapping, ...
    'Atoms_Mapping', Atoms_Mapping ...
);
end
