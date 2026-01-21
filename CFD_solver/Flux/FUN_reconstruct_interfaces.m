function [UL_ls,UL_rs,UR_ls,UR_rs] = FUN_reconstruct_interfaces()

% This function reconstructs teh solution at the left and righ faces

% Access globals
global COUNTER SOL GHOST MESH

N = MESH.N;
Nv = SOL.Nv;
i = COUNTER.pos;
Ng = GHOST.ng;

% Find the conservative variables from the required cells according to the stencil

if i<Ng+1  % ghost cells required at left boundary

    % left face
    % Uleft(Ng-i-1,Nv) = [GHOST.left.U(1:Ng+i-1,Nv);SOL.U(1:i-1,Nv)];
    % Uright(Ng-i-1,Nv) = SOL.U(i:Ng,Nv);
    % [UL_ls,UL_rs] = FUN_constant_reconstruction(Uleft,Uright);
    % Need to make it general for any stencil

    Ulleft = GHOST.left.U;
    Ulright = SOL.U(:,i);

    % right face
    Urleft = SOL.U(:,i);
    Urright = SOL.U(:,i+1);

elseif i>N-Ng  % ghost cells required for right boundary

    % left face
    Ulleft = SOL.U(:,i-1);
    Ulright = SOL.U(:,i);

    % right face
    Urleft = SOL.U(:,i);
    Urright = GHOST.right.U;

else  % only interior cells required

    % left face
    Ulleft = SOL.U(:,i-1);
    Ulright = SOL.U(:,i);

    % right face
    Urleft = SOL.U(:,i);
    Urright = SOL.U(:,i+1);

end

% Apply the reconstruction procedure

% constant reconstruction
[UL_ls,UL_rs] = FUN_constant_reconstruction(Ulleft,Ulright);
[UR_ls,UR_rs] = FUN_constant_reconstruction(Urleft,Urright);


% Placeholder for higher order schemes

end
