function geom = FUN_load_area_profile(fname)
% Loads Area.dat containing [x, r]

% access mesh
global MESH

% Read and load data
D = readmatrix(fname);

x = D(:,1);
r = D(:,2);

geom.x_raw = x;
geom.r_raw = r;
geom.A_raw = pi .* r.^2;

% compute gradient
geom.A = interp1(geom.x_raw, geom.A_raw, MESH.xc, 'pchip', 'extrap');
geom.dlnAdx = gradient(log(geom.A), MESH.xc);

end
