function [FL,FR] = FUN_get_inv_flux()

% This function computes the left and right side inviscid fluxes

% Reconstruct solution at the faces
[Ul_lf, Ur_lf, Ul_rf, Ur_rf] = FUN_reconstruct_interfaces();

% Compute fluxes

% Roe reconstruction for LTE
[FL,~] = FUN_flux_num_roe_lte(Ul_lf, Ur_lf);
[FR,~] = FUN_flux_num_roe_lte(Ul_rf, Ur_rf);


% Need to add more flux schemes and sanity checks for fluid model used

end