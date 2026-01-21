function filtered_reactions = FUN_get_reactions(file, elements)
% FUN_get_reactions  Filter reactions to those involving only the given elements
% (and optionally e- if included in elements).
%
% - Caches the file content to avoid re-reading.
% - Robust tokenization: removes + signs, stoich coefficients, arrows, etc.
% - Ignores third-body 'M' tokens.
%
% INPUTS:
%   file     : path to Reactions.dat
%   elements : cell array of allowed elements, e.g. {'N','O'} or {'H','O','N'}
%              include 'e-' if you want to allow electron in reactions
%
% OUTPUT:
%   filtered_reactions : cell array of reaction lines (strings)

    persistent cache
    if nargin < 2
        error('FUN_get_reactions:MissingInputs', 'Need file and elements.');
    end

    % Electron allowed?
    include_electron = ismember('e-', elements);

    % Read + cache raw lines
    info = dir(file);
    if isempty(info)
        error('FUN_get_reactions:FileNotFound', 'File not found: %s', file);
    end

    key = string(which(file));
    if key == "", key = string(file); end
    cacheField = matlab.lang.makeValidName(key);

    if ~isempty(cache) && isfield(cache, cacheField) && isequal(cache.(cacheField).datenum, info.datenum)
        lines = cache.(cacheField).lines;
    else
        fid = fopen(file, 'r');
        if fid == -1
            error('FUN_get_reactions:OpenFail', 'Cannot open file: %s', file);
        end

        lines = {};
        line = fgetl(fid);
        while ischar(line)
            line = strtrim(line);

            % skip empty/comment lines
            if isempty(line) || startsWith(line, '#') || startsWith(line, '!')
                line = fgetl(fid);
                continue;
            end

            lines{end+1} = line; %#ok<AGROW>
            line = fgetl(fid);
        end
        fclose(fid);

        cache.(cacheField) = struct('datenum', info.datenum, 'lines', {lines});
    end

    % Filter reactions
    filtered_reactions = {};
    for k = 1:numel(lines)
        line = lines{k};

        % Must contain reversible arrow "<->" (as in your file)
        parts = strsplit(line, '<->');
        if numel(parts) ~= 2
            continue; % skip malformed
        end

        % Extract species tokens from both sides
        speciesTokens = [extract_species_tokens(parts{1}), extract_species_tokens(parts{2})];

        % Validate each token
        valid_reaction = true;
        for i = 1:numel(speciesTokens)
            sp = speciesTokens{i};

            % Ignore third-body "M"
            if strcmp(sp, 'M')
                continue;
            end

            % Electron handling
            if strcmp(sp, 'e-')
                if ~include_electron
                    valid_reaction = false;
                    break;
                else
                    continue;
                end
            end

            % Element membership check
            if ~is_member_of_elements(sp, elements)
                valid_reaction = false;
                break;
            end
        end

        if valid_reaction
            filtered_reactions{end+1} = line; %#ok<AGROW>
        end
    end
end

% ---------- helpers ----------

function tokens = extract_species_tokens(sideStr)
% Extract species tokens from one side of reaction
% Removes '+', strips leading stoich coefficients like "2", "3", etc.
% Keeps things like "N2", "NO", "O2+", "N2p" as a single token.

    s = strtrim(sideStr);

    % Split by whitespace
    raw = regexp(s, '\s+', 'split');

    tokens = {};
    for j = 1:numel(raw)
        t = strtrim(raw{j});
        if isempty(t) || strcmp(t, '+')
            continue;
        end

        % Remove leading stoichiometric coefficients: "2N2" or "3 NO"
        % Handle case where coeff and species are stuck together: "2N2"
        t = regexprep(t, '^\d+(\.\d+)?', '');

        % If after stripping it's empty or just '+', skip
        t = strtrim(t);
        if isempty(t) || strcmp(t, '+')
            continue;
        end

        % Remove trailing commas/semicolons (if any)
        t = regexprep(t, '[,;]+$', '');

        % Normalize
        tokens{end+1} = t; %#ok<AGROW>
    end
end

function isValid = is_member_of_elements(species, elements)
% Return true if the species contains only allowed element symbols.
% Note: `elements` should NOT include 'e-' here (handled earlier).
%
% Examples:
% - species = 'N2' with elements {'N','O'} -> ok
% - species = 'NO' with elements {'N','O'} -> ok
% - species = 'CO' with elements {'N','O'} -> reject
%
% This parses element symbols via regex like: [A-Z][a-z]*

    % Special cases that aren't chemical elements
    if strcmp(species, 'M')
        isValid = true;
        return;
    end

    % Remove charge markers or naming conventions if your file uses them.
    % Keep this conservative—adjust if your mechanism uses 'N2p' vs 'N2+'.
    % Example: convert trailing 'p' to '+' semantics: N2p -> N2
    species_clean = species;

    % Common mechanism convention: 'p' for positive ions (N2p, NOp, Op, etc.)
    species_clean = regexprep(species_clean, 'p$', '');

    % Remove explicit '+' if present (O2+, NO+)
    species_clean = strrep(species_clean, '+', '');

    % Now extract element symbols
    individual_elements = regexp(species_clean, '[A-Z][a-z]*', 'match');

    if isempty(individual_elements)
        % If we can't parse elements, treat as invalid (safer)
        isValid = false;
        return;
    end

    isValid = all(ismember(individual_elements, elements));
end
