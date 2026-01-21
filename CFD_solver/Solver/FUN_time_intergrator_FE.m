function FUN_time_intergrator_FE()

% This function advances the solution in time using forward Euler scheme

% Access solution

global SOL MESH
Unew = zeros(SOL.Nv,MESH.N);

% Advance solution

for  i = 1:SOL.Nv
    Unew(i,:) = SOL.U(i,:) + SOL.dt(1,:).*SOL.RHS(i,:);
end

SOL.res = SOL.RHS;

SOL.U = Unew;

end