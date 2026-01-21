function [species_matrix, indicator_matrix, coeff_matrix] = FUN_get_Stoichiometric_Coeff(react)



%% Get the reactions based on the input elements
% react = FUN_get_reactions(elements);

% Initialize a cell array to hold the species
species_matrix = {};

% Initialize the indicator matrix as a character array
num_reactions = length(react);
indicator_matrix = repmat(' ', num_reactions, 10); % Preallocate with spaces

% Initialize the coefficient matrix
coeff_matrix = zeros(num_reactions, 10); % Preallocate for coefficients

% Loop through each reaction
for i = 1:length(react)
    % Split the reaction into reactants and products
    parts = split(react{i}, ' <-> ');

    % Extract reactants and products
    reactants = strtrim(parts{1}); % Trim any whitespace
    products = strtrim(parts{2});   % Trim any whitespace

    % Split reactants and products into individual species
    reactant_species = strsplit(reactants, ' + '); % Split by space and plus
    product_species = strsplit(products, ' + ');   % Split by space and plus

    % Combine reactants and products into a single array for this reaction
    species_matrix{i} = [reactant_species, product_species]; % Combine both species

    % Create indicators for reactants ('r') and products ('p')
    indicators = [repmat('r', 1, length(reactant_species)), ...
        repmat('p', 1, length(product_species))]; % 'r' for reactants, 'p' for products

    % Fill the indicator matrix with indicators
    indicator_matrix(i, 1:length(indicators)) = indicators; % Store in the indicator matrix

    % Get coefficients for reactants
    for j = 1:length(reactant_species)
        % Extract the coefficient as a number
        coeff = str2double(reactant_species{j}(1)); % Get the first character as coefficient
        if isnan(coeff)
            coeff = 1; % Default to 1 if not specified
        end
        coeff_matrix(i, j) = coeff; % Store coefficient in the matrix
    end

    % Get coefficients for products
    for j = 1:length(product_species)
        % Extract the coefficient as a number
        coeff = str2double(product_species{j}(1)); % Get the first character as coefficient
        if isnan(coeff)
            coeff = 1; % Default to 1 if not specified
        end
        coeff_matrix(i, length(reactant_species) + j) = coeff; % Store coefficient in the matrix
    end
end

% Remove coefficients from species matrix
for i = 1:num_reactions
    for j = 1:length(species_matrix{i})
        species_matrix{i}{j} = regexprep(species_matrix{i}{j}, '^\d+', ''); % Remove leading digits
    end
end
indicator_matrix=strtrim(indicator_matrix);
coeff_matrix=coeff_matrix(:,1:size(indicator_matrix,2));

end