function BC = FUN_set_boundary_values()

BC = struct();

% Left boundary (inlet) — simplest first version: fix primitive state
BC.left.type = 'pressure_inlet';
BC.left.P0 = 101325;
BC.left.P = 101325;
BC.left.T0 = 1000;
BC.left.dir = [1,0,0];

% Right boundary (outlet) — simplest first version: fixed pressure
BC.right.type = 'pressure_outlet';
BC.right.P = 2000;
BC.right.dir = [1,0,0];

end