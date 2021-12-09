%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute reflected shortwave radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SWout,A] = func_flux_SWout(C,time,A,SWin)
%% Albedo
soil_cond = A.subSOIL(:,1)==1;
ice_cond = A.subD(:,1)==C.Dice | A.mbal_snow==0;
snow_cond = A.subD(:,1)<C.Dice  ;  

% For a snow surface
A.tstar(A.Tsurf==C.T0 & snow_cond) = C.tstar_wet;
A.tstar(A.Tsurf<C.T0 & snow_cond) = C.tstar_dry + min(C.T0-A.Tsurf(A.Tsurf<C.T0 & snow_cond),10)*C.tstar_K;
A.alb_snow(snow_cond & A.timelastsnow<time.TCUR) = A.alb_snow(snow_cond & A.timelastsnow<time.TCUR) - max(A.alb_snow(snow_cond & A.timelastsnow<time.TCUR)-C.alb_firn,0.0)./A.tstar(snow_cond & A.timelastsnow<time.TCUR)*time.dt;
A.alb_snow(A.timelastsnow==time.TCUR | ice_cond | soil_cond) = C.alb_fresh; 
A.alb(snow_cond,1) = A.alb_snow(snow_cond,1);

% For an ice surface
A.alb(ice_cond,1) = C.alb_ice;

% For a soil surface
A.alb(soil_cond,1) = C.soil_albedo;

%% Reflected SW radiation
SWout = SWin .* A.alb;
end






% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Compute reflected shortwave radiation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% function [SWout,A] = func_flux_SWout(C,time,A,SWin)
% 
% %% Albedo
% soil_cond = A.subSOIL(:,1)==1;
% ice_cond = A.subD(:,1)>=C.Dice-10; %-10 density correction to avoid rimming effect on albedo
% snow_cond = A.subD(:,1)<C.Dice-10; %-10 density correction to avoid rimming effect on albedo
% 
% % For a snow surface
% A.tstar(A.Tsurf==C.T0 & snow_cond) = C.tstar_wet;
% A.tstar(A.Tsurf<C.T0 & snow_cond) = C.tstar_dry + min(C.T0-A.Tsurf(A.Tsurf<C.T0 & snow_cond),10)*C.tstar_K;
% A.alb_snow(snow_cond & A.timelastsnow<time.TCUR) = A.alb_snow(snow_cond & A.timelastsnow<time.TCUR) - (A.alb_snow(snow_cond & A.timelastsnow<time.TCUR)-C.alb_firn)./A.tstar(snow_cond & A.timelastsnow<time.TCUR)*time.dt;
% A.alb_snow(snow_cond & A.timelastsnow==time.TCUR) = C.alb_fresh;
% A.alb_snow(ice_cond | soil_cond) = C.alb_fresh;
% A.alb(snow_cond,1) = A.alb_snow(snow_cond,1);
% 
% % For an ice surface
% A.alb(ice_cond,1) = C.alb_ice;
% 
% % For a soil surface
% A.alb(soil_cond,1) = C.soil_albedo;
% 
% %% Reflected SW radiation
% SWout = SWin .* A.alb;
% 
% end

