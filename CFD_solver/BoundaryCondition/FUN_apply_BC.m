function FUN_apply_BC()
% Applies boundary conditions and updates ghost states in global GHOST.
% This is the ONLY place where ghost states should be updated.

global SOL MESH GHOST BC

N  = MESH.N;
ng = GHOST.ng;

% Looping over boundaries
boundaries = fieldnames(BC);

for b = 1:numel(boundaries)

    side = boundaries{b};
    type = BC.(side).type;

    switch type

        case 'pressure_inlet'

            switch side
                case 'left'
                    Uin = SOL.U(:,1:ng);
                    GHOST.left.U = FUN_bc_pressure_inlet(Uin, GHOST.left.U, BC.left);
                case 'right'
                    Uin = SOL.U(:,N:-1:N-ng+1);
                    GHOST.right.U = FUN_bc_pressure_inlet(Uin, GHOST.right.U, BC.right);
            end

        case 'pressure_outlet'
            switch side
                case 'left'
                    Uin = SOL.U(:, 1:ng);
                    GHOST.left.U = FUN_bc_pressure_outlet(Uin, GHOST.left.U, BC.left);
                case 'right'
                    Uin = SOL.U(:, N:-1:N-ng+1);
                    GHOST.right.U = FUN_bc_pressure_outlet(Uin, GHOST.right.U, BC.right);
            end

        otherwise
            error('BC type "%s" not implemented for %s boundary.', type, side);
    end

end

end
