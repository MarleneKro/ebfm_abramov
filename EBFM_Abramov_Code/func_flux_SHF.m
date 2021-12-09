%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute sensible heat flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SHF] = func_flux_SHF(C,Tsurf,clim,cond)

% Sensible heat flux (bulk equations)
C_kat = max(C.k_turb .* (clim.T(cond)-Tsurf) .* sqrt(C.g./(C.T0 .* clim.Theta_lapse .* C.Pr)),0);
C_turb = 0.5*(C.turb + C_kat);
SHF = clim.Dair(cond).*C.Cp.*C_turb.*(clim.T(cond)-Tsurf);

end