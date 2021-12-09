%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute latent heat flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LHF] = func_flux_LHF(C,Tsurf,clim,cond)

% Turbulent exchange coefficient [Klok & Oerlemans, 2002]
C_kat = max(C.k_turb .* (clim.T(cond)-Tsurf) .* sqrt(C.g./(C.T0 .* clim.Theta_lapse .* C.Pr)),0);
C_turb = 0.5*(C.turb + C_kat);

% Saturation vapor pressure (Clausius-Clapeyron)
VPsurf =  C.VP0*exp(C.Ls/C.Rv/273.15*(1.0-273.15./Tsurf)) .* (Tsurf<273.15) ...
    + C.VP0*exp(C.Lv/C.Rv/273.15*(1.0-273.15./Tsurf)) .* (Tsurf>=273.15);

% Latent heat flux (bulk equations)
LHF = C.eps.*clim.Dair(cond).*C.Ls.*C_turb.*(clim.VP(cond)-VPsurf)./clim.Pres(cond) .* (Tsurf<273.15) ...
    + C.eps.*clim.Dair(cond).*C.Lv.*C_turb.*(clim.VP(cond)-VPsurf)./clim.Pres(cond) .* (Tsurf>=273.15);

end