function out = FUN_get_eq_prop(cmd, varargin)

% EQ  Equilibrium table accessor (global via persistent interpolants)
%
% Current testing signature:
%   out = EQ(cmd, table, var1, var2)
%
% Commands (each takes two inputs):
%   EQ('init',      {T_vec,P_vec,tab}, [],  [])
%   EQ('from_pt',   table, P,   T)
%   EQ('from_rhot', table, rho, T)
%   EQ('from_rhop', table, rho, P)
%   EQ('from_ph',   table, P,   h)
%   EQ('from_rhoh', table, rho, h)
%   EQ('from_pe',   table, P,   e)
%   EQ('from_rhoe', table, rho, e)
%
% Output fields (consistent across commands):
%   P, T, rho, Rmix, h, e, gamma, Mmix, (optional) X
%
% Notes:
% - Interpolants are built in (T, log10(P)) with linear + nearest extrap.
% - Reverse lookups use bracketed fzero with fminbnd fallback if not bracketed.

    persistent interp meta  % interp: griddedInterpolant objects, meta: bounds

    cmd = lower(cmd);

    % ------------------------------------------------------------
    % INIT
    % ------------------------------------------------------------
    if strcmp(cmd,'initialize')
        if numel(varargin) ~= 1
            error('INIT requires exactly one argument: table struct.');
        end
        table = varargin{1};

        T_vec = table.T;         % 1xNT
        P_vec = table.P;         % 1xNP

        logP = log10(P_vec(:).');  % 1 x NP
        T    = T_vec(:);           % NT x 1
        method = 'linear';
        extrap = 'nearest';

        meta.Tmin = min(T_vec); meta.Tmax = max(T_vec);
        meta.Pmin = min(P_vec); meta.Pmax = max(P_vec);

        % Required fields in tab: Rmix, h, e, gamma, Mmix (NT x NP)
        interp.Rmix  = griddedInterpolant({T, logP}, table.R,     method, extrap);
        interp.h     = griddedInterpolant({T, logP}, table.h,     method, extrap);
        interp.e     = griddedInterpolant({T, logP}, table.e,     method, extrap);
        interp.gamma = griddedInterpolant({T, logP}, table.gamma, method, extrap);
        interp.Mmix  = griddedInterpolant({T, logP}, table.M,     method, extrap);

        % Optional: species mole fractions, tab.X is NT x NP x Ns
        if isfield(table,'X')
            Ns = size(table.X,1);           % 5
            interp.X = cell(1,Ns);
            for s = 1:Ns
                % squeeze gives [NT x NP]
                interp.X{s} = griddedInterpolant({T, logP}, squeeze(table.X(s,:,:)), method, extrap);
            end
        end

        out = [];
        return;
    end


    % ------------------------------------------------------------
    % CLEAR
    % ------------------------------------------------------------
    if any(strcmp(cmd, {'clear','finalize'}))
        interp = [];
        meta   = [];
        out    = [];
        return;
    end


    % ------------------------------------------------------------
    % Check number of arguments and require initialized for non-initialize/finalize calls
    % ------------------------------------------------------------
    if numel(varargin) ~= 2
        error('Command "%s" requires exactly two inputs.', cmd);
    end

    var1 = varargin{1};
    var2 = varargin{2};
    if isempty(interp)
        error('EQ is not initialized. Call EQ(''init'', {T_vec,P_vec,tab}, [], []) first.');
    end

    % ------------------------------------------------------------
    % DISPATCH: determine (P,T) then compute full state
    % ------------------------------------------------------------
    switch cmd

        case 'from_pt'
            P = clampP(var1, meta);
            T = var2;

            out = fillFromPT(P, T, interp, meta);

            % Enforce exact inputs (optional but useful for consistency)
            out.P = P;
            out.T = T;

        case 'from_rhot'
            rho = var1;
            T   = var2;

            P = invertP_from_rhoT(rho, T, interp, meta);
            out = fillFromPT(P, T, interp, meta);

            out.rho = rho;

        case 'from_rhop'
            rho = var1;
            P   = clampP(var2, meta);

            T = invertT_from_rhoP(rho, P, interp, meta);
            out = fillFromPT(P, T, interp, meta);

            out.rho = rho;

        case 'from_ph'
            P = clampP(var1, meta);
            h = var2;

            T = invertT_from_hP(h, P, interp, meta);
            out = fillFromPT(P, T, interp, meta);

            out.h = h;

        case 'from_rhoh'
            rho = var1;
            h   = var2;

            T = invertT_from_rhoh(rho, h, interp, meta);
            P = invertP_from_rhoT(rho, T, interp, meta);

            out = fillFromPT(P, T, interp, meta);

            out.rho = rho;
            out.h   = h;

        case 'from_pe'
            P = clampP(var1, meta);
            e = var2;

            T = invertT_from_eP(e, P, interp, meta);
            out = fillFromPT(P, T, interp, meta);

            out.e = e;

        case 'from_rhoe'
            rho = var1;
            e   = var2;

            T = invertT_from_rhoe(rho, e, interp, meta);
            P = invertP_from_rhoT(rho, T, interp, meta);

            out = fillFromPT(P, T, interp, meta);

            out.rho = rho;
            out.e   = e;
        otherwise
            error('Unknown cmd: %s', cmd);
    end

end

% ======================================================================
% Helpers (local functions)
% ======================================================================

function out = fillFromPT(P, T, interp, meta)
    % Clamp P here as a last line of defense
    P  = clampP(P, meta);
    lP = log10(P);

    out.P = P;
    out.T = T;

    out.Rmix  = interp.Rmix(T, lP);
    out.h     = interp.h(T,    lP);
    out.e     = interp.e(T,    lP);
    out.gamma = interp.gamma(T,lP);
    out.Mmix  = interp.Mmix(T, lP);

    % EOS
    out.rho = P ./ (out.Rmix .* T);

    % Optional species
    if isfield(interp,'X')
        Ns = numel(interp.X);
        Xout = zeros([size(T), Ns]);
        for s = 1:Ns
            Xout(:,:,s) = interp.X{s}(T, lP);
        end
        out.X = Xout;
    end
end

function P = clampP(P, meta)
    % Prevent log10(0) and keep within table domain
    P = min(max(P, meta.Pmin), meta.Pmax);
end

function T = clampT(T, meta)
    T = min(max(T, meta.Tmin), meta.Tmax);
end

% ------------------- Inversions: T from (rho,P) ------------------------

function T = invertT_from_rhoP(rho, P, interp, meta)
    P  = clampP(P, meta);
    lP = log10(P);

    f = @(T) (P ./ (interp.Rmix(T, lP).*T)) - rho;
    T = robustSolveInT(f, meta);
end

% ------------------- Inversions: T from (h,P) --------------------------

function T = invertT_from_hP(h, P, interp, meta)
    P  = clampP(P, meta);
    lP = log10(P);

    f = @(T) interp.h(T, lP) - h;
    T = robustSolveInT(f, meta);
end

% ------------------- Inversions: T from (e,P) --------------------------

function T = invertT_from_eP(e, P, interp, meta)
    P  = clampP(P, meta);
    lP = log10(P);

    f = @(T) interp.e(T, lP) - e;
    T = robustSolveInT(f, meta);
end

% ------------------- Inversions: P from (rho,T) ------------------------
% Solve for P in [Pmin,Pmax] at fixed T:
%   g(P) = P/(R(T,P)T) - rho = 0
% Solve in log10(P) for scaling

function P = invertP_from_rhoT(rho, T, interp, meta)
    T = clampT(T, meta);

    g = @(lp) (10.^lp) ./ (interp.Rmix(T, lp).*T) - rho;

    lpmin = log10(meta.Pmin);
    lpmax = log10(meta.Pmax);

    lp = robustSolveInLogP(g, lpmin, lpmax);
    P  = 10.^lp;

    P = clampP(P, meta);
end

% ------------------- Inversions: (rho,h) -> T --------------------------

function T = invertT_from_rhoh(rho, htar, interp, meta)
    f = @(T) h_at_rhoT(T, rho, interp, meta) - htar;
    T = robustSolveInT(f, meta);
end

function hval = h_at_rhoT(T, rho, interp, meta)
    T = clampT(T, meta);
    P = invertP_from_rhoT(rho, T, interp, meta);
    hval = interp.h(T, log10(P));
end

% ------------------- Inversions: (rho,e) -> T --------------------------

function T = invertT_from_rhoe(rho, etar, interp, meta)
    f = @(T) e_at_rhoT(T, rho, interp, meta) - etar;
    T = robustSolveInT(f, meta);
end

function evalv = e_at_rhoT(T, rho, interp, meta)
    T = clampT(T, meta);
    P = invertP_from_rhoT(rho, T, interp, meta);
    evalv = interp.e(T, log10(P));
end

% ------------------- Robust bracketed solvers --------------------------

function T = robustSolveInT(f, meta)
    Tmin = meta.Tmin;
    Tmax = meta.Tmax;

    f1 = f(Tmin);
    f2 = f(Tmax);

    if ~isfinite(f1) || ~isfinite(f2)
        T = fminbnd(@(T) abs(f(T)), Tmin, Tmax);
        return;
    end

    if sign(f1) == sign(f2)
        % Not bracketed (target outside range or mild non-monotonicity)
        T = fminbnd(@(T) abs(f(T)), Tmin, Tmax);
    else
        T = fzero(f, [Tmin, Tmax]);
    end
end

function lp = robustSolveInLogP(g, lpmin, lpmax)
    g1 = g(lpmin);
    g2 = g(lpmax);

    if ~isfinite(g1) || ~isfinite(g2)
        lp = fminbnd(@(lp) abs(g(lp)), lpmin, lpmax);
        return;
    end

    if sign(g1) == sign(g2)
        lp = fminbnd(@(lp) abs(g(lp)), lpmin, lpmax);
    else
        lp = fzero(g, [lpmin, lpmax]);
    end
end
