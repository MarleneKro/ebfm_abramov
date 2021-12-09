%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update the mass balance and snow mass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A] = func_massbalance(A,clim,C)

%% Update the mass balance
A.smb = clim.snow + clim.rain - A.runoff ...
    + A.moist_deposition + A.moist_condensation ...
    - A.moist_sublimation - A.moist_evaporation;

A.mbal_snow = max(A.mbal_snow+A.smb,0);
A.mbal_snow(all(A.subD==C.Dice,2),1) = 0;

%set mbal_snow below threshold to 0 to allow for alb_ice
A.mbal_snow(A.mbal_snow<=C.mbal_snow_thresh) = 0;

A.mbal = A.mbal + A.smb;
A.smb_stake = clim.snow - A.melt + A.refr_seasnl_snow ...
    + A.moist_deposition - A.moist_sublimation;
A.mbal_stake = A.mbal_stake + A.smb_stake; 

end


