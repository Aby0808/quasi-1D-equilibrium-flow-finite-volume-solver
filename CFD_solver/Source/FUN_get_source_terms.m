function S = FUN_get_source_terms()

%This functions computes the source terms for teh given fluid model

% For now it just computes the source terms for area change

global GEOM SOL COUNTER

S = zeros(SOL.Nv,1);

i = COUNTER.pos;

rho = SOL.U(1,i);
E = SOL.U(3,i);

out = FUN_get_eq_prop('from_rhoe',rho, E);

% calculate slope
slope = GEOM.dlnAdx(i);
S(2) = out.P*slope;

end