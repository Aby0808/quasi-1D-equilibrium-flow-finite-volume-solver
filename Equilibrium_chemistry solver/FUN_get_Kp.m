function [Kp, Species, type, Coeffs, unique_species_list] = FUN_get_Kp(Species_Data, Species_Mapping, reactions, T)

% Get the stoichiometric coefficients, species, and their types
[Species, type, Coeffs] = FUN_get_Stoichiometric_Coeff(reactions);

num_react = size(Species, 2);  % Number of reactions
Kp = zeros(num_react, 1);      % Preallocate Kp array
G0_reactants = 0;
G0_products = 0;

unique_species_set = {};  % To store the unique species

% Loop through each reaction
for i = 1:num_react
    nb_sp_react = length(Species{i});  % Number of species in the current reaction
    
    % Loop through each species in the current reaction
    for j = 1:nb_sp_react
        species_current = char(Species{i}(j));  % Current species name
        species_type = type(i, j);              % Reactant ('r') or Product ('p')
        species_coeff = Coeffs(i, j);           % Stoichiometric coefficient

        % Extract the Gibbs free energy for the current species
        [~, ~, ~, G0] = FUN_extract_species_thermo_data(Species_Data, Species_Mapping, T, species_current);

        % Accumulate Gibbs free energy for reactants or products
        if species_type == 'r'
            G0_reactants = G0_reactants + species_coeff * G0;
        elseif species_type == 'p'
            G0_products = G0_products + species_coeff * G0;
        end

        % Add the current species to the unique species set if not already present
        if ~ismember(species_current, unique_species_set)
            unique_species_set = [unique_species_set, species_current];
        end
    end
    
    % Calculate ΔG0 for the current reaction
    delG0 = (G0_products - G0_reactants);
    
    % Calculate Kp for the current reaction
    Kp(i) = exp(-delG0 / (8.314 * T));

    % Reset the Gibbs free energy accumulators for the next reaction
    G0_reactants = 0;
    G0_products = 0;
end

% Return the unique species list
unique_species_list = unique_species_set;

end
