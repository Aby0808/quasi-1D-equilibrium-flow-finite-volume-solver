function [U_L, U_R] = FUN_constant_reconstruction(Uleft, Uright)
% Piecewise-constant (1st-order) reconstruction at a single interface

% Safety checks
if isempty(Uleft) || isempty(Uright)
    error('FUN_constant_reconstruction: empty stencil provided.');
end

% Constant reconstruction:
% take nearest cell values to the interface
U_L = Uleft(:,end);
U_R = Uright(:,1);

end
