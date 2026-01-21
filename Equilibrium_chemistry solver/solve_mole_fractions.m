function [X] = solve_mole_fractions(x0,P,T,element_fract,table_file_name, reaction_file_name, elements)
    % solve_mole_fractions: Solve for the mole fractions given initial guesses
    % x0 - initial guesses for the mole fractions (array)
    % P - total pressure (atm)
    % Kp - equilibrium constants for the reactions (array)

    % Lower bounds for the mole fractions (mole fractions must be >= 0)
    lb = zeros(length(x0), 1); % Set lower bounds to 0

    % Get equilibrium constant and properties
    [Species_Data, Species_Mapping, ~] = FUN_read_data(table_file_name);
    reactions = FUN_get_reactions(reaction_file_name, elements);
    [Kp,~,~,~,~] = FUN_get_Kp(Species_Data, Species_Mapping, reactions, T);

    % Call lsqnonlin to solve the system of nonlinear equations
    % Create an anonymous function to pass P and Kp to mole_fraction_equations
    options = optimoptions('lsqnonlin', 'Display', 'none');
    P = P/101325; %normalized pressure
    [X] = lsqnonlin(@(x) mole_fraction_equations(x, P, Kp,element_fract), x0, lb, [], options);
    
    % Display results
    % disp('Solution for the mole fractions:');
    % disp(['XN2   = ', num2str(X(1))]);
    % disp(['XN    = ', num2str(X(2))]);
    % disp(['XO2   = ', num2str(X(3))]);
    % disp(['XO    = ', num2str(X(4))]);
    % disp(['XNO   = ', num2str(X(5))]);
end

% Function defining the system of nonlinear equations
% function F = mole_fraction_equations(x, P, Kp,ele_frac)
%     % x: Array containing the mole fractions of each species
%     % P: Total pressure (atm)
%     % Kp: Equilibrium constants for each reaction
% 
%     % Mole fractions of species
%     XN2 = x(1);
%     XN = x(2);
%     XO2 = x(3);
%     XO = x(4);
%     XNO = x(5);
% 
%     % Define the system of equations
%     % Kp(i) corresponds to the equilibrium constant for the i-th reaction
%     F(1) = 1 - log((XNO^2) / (XN2 * XO2)) / log(Kp(3));   % Reaction 1 equilibrium relation
%     F(2) = 1 - log((XN2) / (XN^2 * P)) / log(Kp(2));      % Reaction 2 equilibrium relation
%     F(3) = 1 - log((XO^2 * P) / (XO2)) / log(Kp(1));      % Reaction 3 equilibrium relation
% 
%     % Mole fraction balance (sum of mole fractions = 1)
%     F(4) = 1 - (XN2 + XN + XO2 + XO + XNO);
% 
%     % Elemental balance (N and O balance example)
%     % Modify this equation based on the elemental fraction you want to enforce
%     F(5) = ele_frac - (2 * XN2 + XN + XNO) / (2 * XO2 + XO + XNO);
% end


function F = mole_fraction_equations(x, P, Kp, ele_frac)
    % x: Array containing the mole fractions of each species
    % P: Total pressure (atm)
    % Kp: Equilibrium constants for each reaction

    % Mole fractions of species
    XN2 = x(1);
    XN  = x(2);
    XO2 = x(3);
    XO  = x(4);
    XNO = x(5);

    % Define the system of equations
    % Kp(i) corresponds to the equilibrium constant for the i-th reaction
    F(1) = 1 - log((XNO^2) / (XN2 * XO2)) / log(Kp(3));   % Reaction 1 equilibrium relation
    F(2) = 1 - log((XN2) / (XN^2 * P)) / log(Kp(2));      % Reaction 2 equilibrium relation
    F(3) = 1 - log((XO^2 * P) / (XO2)) / log(Kp(1));      % Reaction 3 equilibrium relation

    % Mole fraction balance (sum of mole fractions = 1)
    w = 1e3;  % weight to enforce sum(X)=1 strongly (try 1e2 to 1e6 if needed)
    F(4) = w * (1 - (XN2 + XN + XO2 + XO + XNO));

    % Elemental balance (N and O balance example)
    F(5) = ele_frac - (2 * XN2 + XN + XNO) / (2 * XO2 + XO + XNO);
end
