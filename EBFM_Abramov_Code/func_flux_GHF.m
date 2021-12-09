%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute subsurface heat flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [GHF] = func_flux_GHF(Tsurf,A,GHF_C,cond)

GHF = GHF_C(cond).*(A.subT(cond,3)-Tsurf(:));

end