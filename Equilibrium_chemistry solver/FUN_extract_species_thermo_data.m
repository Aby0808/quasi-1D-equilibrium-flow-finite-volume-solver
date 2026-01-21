function [Cp_molar, h_molar, s_molar, g_molar] = FUN_extract_species_thermo_data...
    (Species_Data, Species_Mapping, T, species)
% Returns molar properties:
% Cp_molar [J/mol-K], h_molar [J/mol], s_molar [J/mol-K], g_molar [J/mol]

Ru = 8.314462618; % J/mol-K

Data = Species_Data.(Species_Mapping(species));

% --- Select coefficient range (your table layout) ---
if T < Data{1,1}(2) && T >= Data{1,1}(1)
    a = Data{1,2};  b = Data{1,3};
elseif T < Data{1,4}(2) && T >= Data{1,4}(1)
    a = Data{1,5};  b = Data{1,6};
elseif T <= Data{1,7}(2) && T >= Data{1,7}(1)
    a = Data{1,8};  b = Data{1,9};
else
    error('Thermodynamics data for T = %g K does not exist', T);
end

% In this form, a gives Cp/R, and h/(RT), s/R (NASA-9 style)
Cp_over_R = a(1)*T^-2 + a(2)*T^-1 + a(3) + a(4)*T + a(5)*T^2 + a(6)*T^3 + a(7)*T^4;

h_over_RT = -a(1)*T^-2 + a(2)*log(T)/T + a(3) + a(4)*T/2 + a(5)*T^2/3 + a(6)*T^3/4 + a(7)*T^4/5 + b(1)/T;

s_over_R  = -0.5*a(1)*T^-2 - a(2)*T^-1 + a(3)*log(T) + a(4)*T + a(5)*T^2/2 + a(6)*T^3/3 + a(7)*T^4/4 + b(2);

% Convert to molar SI units
Cp_molar = Cp_over_R * Ru;          % J/mol-K
h_molar  = h_over_RT * Ru * T;      % J/mol
s_molar  = s_over_R  * Ru;          % J/mol-K
g_molar  = h_molar - T*s_molar;     % J/mol
end