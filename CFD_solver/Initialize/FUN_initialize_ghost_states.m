function GHOST = FUN_initialize_ghost_states(nVar, ng, bc)
% Initializes ghost-cell containers and boundary metadata.
%
% Convention: GHOST.<side>.U(k,:) is the k-th ghost cell away from boundary,
% with k=1 adjacent to the boundary.

GHOST = struct();
GHOST.ng   = ng;
GHOST.nVar = nVar;

% Left boundary
GHOST.left = struct();
GHOST.left.U    = zeros(nVar, ng);
GHOST.left.type = bc.left.type;

% Right boundary
GHOST.right = struct();
GHOST.right.U    = zeros(nVar, ng);
GHOST.right.type = bc.right.type;

% (Optional) a stamp to track if it’s updated for current stage/iter
GHOST.stamp = -inf;

end
