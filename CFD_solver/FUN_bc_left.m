function Uext = FUN_bc_left(Uext)
% Default placeholder: transmissive (copy first interior cell)

global MESH
ng = MESH.ng;
i1 = MESH.iC(1);     % first physical cell index in extended storage

for g = 1:ng
    Uext(i1-g,:) = Uext(i1,:);   % copy interior to ghost
end
end
