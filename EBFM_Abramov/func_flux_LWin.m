%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute incoming longwave radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LWin] = func_flux_LWin(C,clim)

% Clear-sky emissivity
ecs = 0.23 + C.b.*(clim.VP./clim.T).^0.125;

% Sky emissivity
e = ecs.*(1.0-clim.C.^C.p) + C.ecl.*clim.C.^C.p;
		
% Incoming longwave radiation
LWin = e .* C.boltz .* clim.T.^4;	

end