function [gamma_mix, Cp_mix, Cv_mix, Rmix, Mmix, X, h_mix, e_mix] = FUN_equilibrium_prop(T, P, x0)

if nargin < 3
    x0 = [0.77, 0.01, 0.19, 0.01, 0.02];
end

table_file_name    = 'Thermodynamic_tables.dat';
reaction_file_name = 'Reactions.dat';
element_fraction   = 3.727;
elements = {'N','O'};

Ru = 8.314462618; % J/mol-K

% Molecular weights (kg/mol)
MW.N2 = 28.0134e-3;
MW.N  = 14.0067e-3;
MW.O2 = 31.9988e-3;
MW.O  = 15.9994e-3;
MW.NO = 30.0061e-3;

Sp      = {'N2','N','O2','O','NO'};
MW_list = [MW.N2, MW.N, MW.O2, MW.O, MW.NO];

% Initial guess for equilibrium solver
% x0 = [0.77, 0.01, 0.19, 0.01, 0.02];

% Read thermo tables once
[Species_Data, Species_Mapping, ~] = FUN_read_data(table_file_name);

% solve for equilibrium mole fractions
X = solve_mole_fractions(x0, P, T, element_fraction, ...
                         table_file_name, reaction_file_name, elements);

X = X(:); % ensure column

% ---- Mixture Cp (molar) and mixture h (molar) ----
Cp_molar_mix = 0.0;   % J/mol-K
h_molar_mix  = 0.0;   % J/mol

for i = 1:numel(Sp)
    % Modify FUN_extract_species_thermo_data to ALSO return Hbar = H/(Ru*T)
    % (see the small change below)
    [Cp_i_molar, ~, Hbar_i] = FUN_extract_species_thermo_data( ...
        Species_Data, Species_Mapping, T, Sp{i});

    Cp_molar_mix = Cp_molar_mix + X(i) * Cp_i_molar;

    h_i_molar = Hbar_i * Ru * T;     % J/mol
    h_molar_mix = h_molar_mix + X(i) * h_i_molar;
end

% ---- Mixture MW and gas constant ----
Mmix = sum(X.' .* MW_list);  % kg/mol
Rmix = Ru / Mmix;            % J/kg-K

% ---- Convert molar -> mass basis ----
Cp_mass_mix = Cp_molar_mix / Mmix;  % J/kg-K
Cv_mass_mix = Cp_mass_mix - Rmix;   % J/kg-K

h_mix = h_molar_mix / Mmix;         % J/kg
e_mix = h_mix - Rmix * T;           % J/kg  (ideal gas: e = h - RT)

gamma_mix = Cp_mass_mix / Cv_mass_mix;

% Optional: return Cp/Cv in kJ/kg-K if you like
Cp_mix = Cp_mass_mix / 1000;  % kJ/kg-K
Cv_mix = Cv_mass_mix / 1000;  % kJ/kg-K
h_mix  = h_mix / 1000;        % kJ/kg
e_mix  = e_mix / 1000;        % kJ/kg

end
