function MESH = FUN_generate_mesh(x0, xe, xt, N, ng)
% Generates 1D finite-volume mesh (physical) + ghost indexing contract
% Faces at x0, xe; cell centers inside domain.
%
% NOTE: Ghost CELLS are an indexing/storage concept, not geometry.
% We store ng + index ranges here, but NOT ghost solution values.

    arguments
        x0 (1,1) double
        xe (1,1) double
        xt (1,1) double
        N  (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(N,3)}
        ng (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(ng,1)} = 1
    end

    if xe <= x0
        error('xe must be greater than x0');
    end

    % ---- clustering strength
    beta = 2.5;

    % normalized throat position in [0,1]
    s_th = (xt - x0)/(xe - x0);
    s_th = min(max(s_th,0),1);

    % parametric coordinate for faces
    sf = linspace(0,1,N+1);

    % clustering towards throat
    sf_mapped = zeros(size(sf));
    for ii = 1:numel(sf)
        if sf(ii) <= s_th
            if s_th == 0
                sf_mapped(ii) = 0;
            else
                s = sf(ii)/s_th;
                sf_mapped(ii) = s_th * tanh(beta*s)/tanh(beta);
            end
        else
            if s_th == 1
                sf_mapped(ii) = 1;
            else
                s = (sf(ii)-s_th)/(1-s_th);
                sf_mapped(ii) = s_th + (1-s_th)*(1 - tanh(beta*(1-s))/tanh(beta));
            end
        end
    end

    % physical faces
    xf = x0 + (xe - x0)*sf_mapped(:);

    % cell centers
    xc = 0.5*(xf(1:end-1) + xf(2:end));

    % cell widths
    dx = diff(xf);

    % throat index (closest center)
    [~, i_throat] = min(abs(xc - xt));

    % ---- store globally
    % global MESH
    MESH.xf = xf;
    MESH.xc = xc;
    MESH.dx = dx;
    MESH.N  = N;
    MESH.x0 = x0;
    MESH.xe = xe;
    MESH.xt = xt;
    MESH.i_throat = i_throat;

    % % ---- ghost indexing contract 
    MESH.ng   = ng;
    % MESH.Ntot = N + 2*ng;
    % 
    % % Indices in an extended array Uext of size (Ntot x Nv)
    % MESH.iL = 1:ng;                 % left ghost cells
    % MESH.iC = (ng+1):(ng+N);        % physical cells inside extended storage
    % MESH.iR = (ng+N+1):(ng+N+ng);   % right ghost cells
    % 
    % % Map physical cell i=1..N to extended index ic
    % MESH.phys2ext = @(i) i + ng;

end
