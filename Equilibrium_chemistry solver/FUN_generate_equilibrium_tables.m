function tables = FUN_generate_equilibrium_tables(Pmin, Pmax, Tmin, Tmax, NP, NT, outFile)
%FUN_generate_equilibrium_tables  Generate equilibrium thermo tables on (T,P) grid.
%
% Inputs:
%   Pmin,Pmax : pressure bounds [Pa]  (use Pa consistently)
%   Tmin,Tmax : temperature bounds [K]
%   NP,NT     : number of grid points
%   outFile   : optional, e.g. 'eqTables.mat'
%
% Output:
%   tbl : struct with fields:
%       .T, .P (vectors)
%       .gamma, .Cp, .Cv, .R, .M, .h, .e (NT x NP arrays)
%       .X (nsp x NT x NP)

    if nargin < 7
        outFile = '';
    end

    % Grids
    T_vec = linspace(Tmin, Tmax, NT);                 % K
    P_vec = logspace(log10(Pmin), log10(Pmax), NP);   % Pa

    % Call once to know nsp and inital guess
    [~,~,~,~,~,X0,~,~] = FUN_equilibrium_prop(T_vec(1), P_vec(1));
    nsp = 5;%numel(X0);

    % Preallocate
    gamma = zeros(NT, NP);
    Cp    = zeros(NT, NP);   % J/kg-K
    Cv    = zeros(NT, NP);   % J/kg-K
    Rmix  = zeros(NT, NP);   % J/kg-K
    Mmix  = zeros(NT, NP);   % kg/mol
    hmix  = zeros(NT, NP);   % J/kg
    emix  = zeros(NT, NP);   % J/kg
    X     = zeros(nsp, NT, NP);

    % Fill table
    % for i = 1:NT
    %     Ti = T_vec(i);
    %     for j = 1:NP
    %         Pj = P_vec(j);
    % 
    %         [gamma(i,j), Cp(i,j), Cv(i,j), Rmix(i,j), Mmix(i,j), X(:,i,j), hmix(i,j), emix(i,j)] = ...
    %             FUN_equilibrium_prop(Ti, Pj, X0);
    %     end
    %     fprintf('Table row %d/%d done (T = %.1f K)\n', i, NT, Ti);
    % end

for j = 1:NP
    Pj = P_vec(j);

    % Initial guess at the first temperature for this pressure
    Xprev = X0(:);

    for i = 1:NT
        Ti = T_vec(i);

        [gamma(i,j), Cp(i,j), Cv(i,j), Rmix(i,j), Mmix(i,j), X(:,i,j), hmix(i,j), emix(i,j)] = ...
            FUN_equilibrium_prop(Ti, Pj, Xprev);

        % Warm start for next temperature
        Xprev = X(:,i,j);
    end

    fprintf('Pressure column %d/%d done (P = %.2e Pa)\n', j, NP, Pj);
end



    % Pack in struct
    tables.T = T_vec;
    tables.P = P_vec;
    tables.gamma = gamma;
    tables.Cp = Cp;
    tables.Cv = Cv;
    tables.R  = Rmix;
    tables.M  = Mmix;
    tables.h  = hmix;
    tables.e  = emix;
    tables.X  = X;

    % Save if requested
    if ~isempty(outFile)
        save(outFile, 'tables', '-v7.3');
        fprintf('Saved equilibrium table to %s\n', outFile);
    end
end