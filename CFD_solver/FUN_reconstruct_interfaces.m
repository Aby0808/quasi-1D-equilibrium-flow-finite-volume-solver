function [UL_imh, UR_imh, UL_iph, UR_iph] = FUN_reconstruct_interfaces(Uext, i)
% Piecewise-constant reconstruction using extended array with ghost cells.
% No boundary branching; BC procedure must have already filled ghost states.
%
% Inputs:
%   Uext : (Ntot x Nv) array with ghosts (rows), conservative vars
%   i    : physical cell index, 1..N
%
% Outputs (all Nv x 1):
%   UL_imh, UR_imh : left/right states at i-1/2
%   UL_iph, UR_iph : left/right states at i+1/2

global MESH
ic = MESH.phys2ext(i);   % extended index for physical cell i

% i-1/2 interface (between ic-1 and ic)
UL_imh = Uext(ic-1,:).';
UR_imh = Uext(ic,:).';

% i+1/2 interface (between ic and ic+1)
UL_iph = Uext(ic,:).';
UR_iph = Uext(ic+1,:).';

end
