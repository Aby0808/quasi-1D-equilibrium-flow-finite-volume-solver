function Uext = FUN_bc_right(Uext)
% Default placeholder: transmissive (copy last interior cell)

global MESH
ng = MESH.ng;
iN = MESH.iC(end);   % last physical cell index in extended storage

for g = 1:ng
    Uext(iN+g,:) = Uext(iN,:);   % copy interior to ghost
end
end
