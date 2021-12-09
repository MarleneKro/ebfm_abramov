%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute energy balance and solve for surface temperature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,insol] = func_energybalance(C,A,clim,insol,t,time,grid)

%% Compute fluxes SWin, SWout, LWin, Qrain and GHF components
[insol]     = func_flux_insol(grid,time,t,insol);
[SWin,insol,A]= func_flux_SWin(C,insol,clim,A); 

[SWout,A]   = func_flux_SWout(C,time,A,SWin);
LWin        = func_flux_LWin(C,clim);

GHF_k       = (0.138-1.01d-3.*A.subD+3.233d-6.*A.subD.^2) .* (A.subSOIL==0) ...
            + C.soil_Kfrozen .* (A.subSOIL==1 & A.subT<273.15) ...
            + C.soil_Kthawed .* (A.subSOIL==1 & A.subT>=273.15);
GHF_C       = (GHF_k(:,1).*A.subZ(:,1)+GHF_k(:,2).*A.subZ(:,2)+0.5.*GHF_k(:,3).*A.subZ(:,3)) ...
                ./(A.subZ(:,1)+A.subZ(:,2)+0.5.*A.subZ(:,3)).^2;

%% Iterative procedure to solve the energy balance for surface temperature:
% - start by setting temperature range
% - compute surface temperature dependent energy fluxes
% - sum all fluxes to obtain the energy budget
% - use bisection method to iteratively update surface temperature
% - stop when surface temperature changes less than predefined limit

tt1 = A.Tsurf-40d0;
tt2 = A.Tsurf+40d0;
dT = tt2-tt1;
for c=1:20																
	dT=0.5*dT;
	ttmid=(tt1+tt2)/2d0;
    
    
    LWout  = func_flux_LWout(C,tt1);
    LHF    = func_flux_LHF(C,tt1,clim,true(grid.gpsum,1));
    SHF    = func_flux_SHF(C,tt1,clim,true(grid.gpsum,1));
    GHF    = func_flux_GHF(tt1,A,GHF_C,true(grid.gpsum,1));

	ebal_tt1 = SWin - SWout + LWin - LWout + SHF + LHF + GHF;
    
    LWout  = func_flux_LWout(C,ttmid);
    LHF    = func_flux_LHF(C,ttmid,clim,true(grid.gpsum,1));
    SHF    = func_flux_SHF(C,ttmid,clim,true(grid.gpsum,1));
    GHF    = func_flux_GHF(ttmid,A,GHF_C,true(grid.gpsum,1));
    
	ebal_ttmid = SWin - SWout + LWin - LWout + SHF + LHF + GHF;
    
    cond_EB = ebal_ttmid.*ebal_tt1<0;
    tt2(cond_EB) = ttmid(cond_EB);
    tt1(~cond_EB) = ttmid(~cond_EB);
  
    if (max(dT)<C.dTacc), break; end
end

%% In case the computed surface temperature is higher than zero:
% - set surface temperature to melting point
% - excess energy is used for melting (Emelt)

ttmid(abs(ttmid-273.15)<C.dTacc) = ttmid(abs(ttmid-273.15)<C.dTacc)-C.dTacc;
ttmid(A.subSOIL(:,1)==0 | clim.snow(:,1)>0) = min(ttmid(A.subSOIL(:,1)==0 | clim.snow(:,1)>0),273.15);

LWout   = func_flux_LWout(C,ttmid);
LHF     = func_flux_LHF(C,ttmid,clim,true(grid.gpsum,1));
SHF     = func_flux_SHF(C,ttmid,clim,true(grid.gpsum,1));
GHF     = func_flux_GHF(ttmid,A,GHF_C,true(grid.gpsum,1));

Emelt   = SWin - SWout + LWin - LWout + SHF + LHF + GHF;% + Qrain;
Emelt(A.subSOIL(:,1)==1 & clim.snow(:,1)==0) = 0.0;
Emelt(ttmid<273.15) = 0.0;

A.melt = 24.*3600.*time.dt.*Emelt./C.Lm./1d3;

%% Avoid melt of soil in case of thin snow cover
max_melt = sum(A.subZ.*A.subD./C.Dwater.*(A.subSOIL==0),2)+clim.snow;
A.melt = min(A.melt,max_melt);
cond = A.melt==max_melt & A.melt>0;
A.Emelt_eff = A.melt(cond) ./ (24.*3600.*time.dt) .* (C.Lm.*1d3);
tt1(cond) = A.Tsurf(cond)-40d0;
tt2(cond) = A.Tsurf(cond)+40d0;
dT = tt2-tt1;

for c=1:20																
	dT=0.5*dT;
	ttmid(cond) = (tt1(cond)+tt2(cond))/2d0;
    
    LWout(cond)  = func_flux_LWout(C,tt1(cond));
    LHF(cond)    = func_flux_LHF(C,tt1(cond),clim,cond);
    SHF(cond)    = func_flux_SHF(C,tt1(cond),clim,cond);
    GHF(cond)    = func_flux_GHF(tt1(cond),A,GHF_C,cond);
    
    ebal_tt1(cond) = SWin(cond) - SWout(cond) + LWin(cond) - LWout(cond) + SHF(cond) + LHF(cond) + GHF(cond) - A.Emelt_eff;
    
    LWout(cond)  = func_flux_LWout(C,ttmid(cond));
    LHF(cond)    = func_flux_LHF(C,ttmid(cond),clim,cond);
    SHF(cond)    = func_flux_SHF(C,ttmid(cond),clim,cond);
    GHF(cond)    = func_flux_GHF(ttmid(cond),A,GHF_C,cond);
	
    ebal_ttmid(cond) = SWin(cond) - SWout(cond) + LWin(cond) - LWout(cond) + SHF(cond) + LHF(cond) + GHF(cond) - A.Emelt_eff;

    cond_EB = ebal_ttmid.*ebal_tt1<0;
    tt2(cond_EB & cond) = ttmid(cond_EB & cond);
    tt1(~cond_EB & cond) = ttmid(~cond_EB & cond);
    
    if (max(dT)<C.dTacc), break; end
end

%% Moisture exchange & melt
A.moist_deposition = 24.*3600.*time.dt.*LHF./C.Ls./1d3 .* (ttmid<273.15 & LHF>0);
A.moist_condensation = 24.*3600.*time.dt.*LHF./C.Lv./1d3 .* (ttmid>=273.15 & LHF>0);
A.moist_sublimation = -24.*3600.*time.dt.*LHF./C.Ls./1d3 .* (ttmid<273.15 & LHF<0);
A.moist_evaporation = -24.*3600.*time.dt.*LHF./C.Lv./1d3 .* (ttmid>=273.15 & LHF<0);

%% Avoid evaporation of absent melt
max_evap = A.melt;
A.moist_evaporation = min(A.moist_evaporation,max_evap);

%% Avoid sublimation of soil
max_subl = sum(A.subZ.*A.subD./C.Dwater.*(A.subSOIL==0),2)+clim.snow;
A.moist_sublimation = min(A.moist_sublimation,max_subl);

%% Store local variables in output structure A 
A.Tsurf = ttmid;
A.Emelt = Emelt;
A.LHF = LHF;
A.TOAshade = insol.TOAshade;
A.SWin = SWin;
A.SWout = SWout;
A.LWin = LWin;
A.LWout = LWout;
A.SHF = SHF;
A.LHF = LHF;
A.GHF = GHF;
A.transm = A.SWin ./ A.TOAshade;


end