%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load climate fields from file or heuristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [clim,A] = func_loadclimate(C,grid,clim,io,t,time,A)

%% Read climate input from file(s)
if (t==1)
    load(io.climfile);
    A.data_clim = climate;
end

tt = (time.TCUR - time.TI)/time.dt + 1;
tt = round(tt); 
%% Temperature [K]
clim.T(:) = A.data_clim(tt,7) + A.data_clim(tt,8)*(grid.z_mask-C.metstat_elev);            
clim.T_lapse = A.data_clim(tt,8);                             


%% Precipitation [m w.e. ts-1]
C.summer_months_cor = [7,8,9];
A.corrP_summer = ismember(A.data_clim(tt,2),C.summer_months_cor);

if A.corrP_summer == 1
   C.prec_corr_f = C.prec_corr_summer;
else
   C.prec_corr_f = C.prec_corr_winter;
end


if A.data_clim(tt,14) > 0                                                                   
    clim.P(:) = (A.data_clim(tt,14)/1d3*24*3600*time.dt)*C.prec_corr_f +...
        (A.data_clim(tt,14)/1d3*24*3600*time.dt)*C.prec_corr_f*C.prec_grad*(grid.z_mask-C.metstat_elev);
else
    clim.P(:) = A.data_clim(tt,14)/1d3 * 24 * 3600 * time.dt;
end

%% Cloud cover [fraction]
clim.C(:) = max(min( A.data_clim(tt,13),1.0),0.0);                                        

%% Relative humidity [fraction]
clim.RH(:) = max(min(A.data_clim(tt,9)/1d2...                                             
    + (grid.z_mask-C.metstat_elev)*A.data_clim(tt,10)/1d2,1.0),0.0);

%% Air pressure [Pa]
clim.Pres(:) = A.data_clim(tt,11) * exp(A.data_clim(tt,12)*(grid.z_mask-C.metstat_elev));   
clim.Pres_lapse = A.data_clim(tt,12);

%% Wind speed [m s-1] (only used for snow densification and only if C.Dfreshsnow < 350.0)
clim.WS(:) = A.data_clim(tt,16);

%% Potential temperature lapse rate [K m-1]
clim.Theta_lapse(:) = A.data_clim(tt,15);                                               

%% Derived climate fields

% Snowfall / Rainfall
clim.snow = clim.P .* (clim.T < C.rainsnowT-1);
clim.rain = clim.P .* (clim.T > C.rainsnowT+1);
clim.snow = clim.snow + clim.P .* (C.rainsnowT-clim.T+1)./2 .* (clim.T < C.rainsnowT+1 & clim.T > C.rainsnowT-1);
clim.rain = clim.rain + clim.P .* (1+clim.T-C.rainsnowT)./2 .* (clim.T < C.rainsnowT+1 & clim.T > C.rainsnowT-1);

% Annual snow accumulation (yearsnow)
A.ys = (1.0-1.0/(365.0/time.dt)).*A.ys + clim.P.*1d3;
logys = log(A.ys);
clim.yearsnow = repmat(A.ys,[1 grid.nl]);
clim.logyearsnow = repmat(logys,[1 grid.nl]);

% Vapor pressure (VP) / specific humidity (q)
VPsat = C.VP0.*exp(C.Lv/C.Rv.*(1.0./273.15-1.0./clim.T)) .* (clim.T>=273.15) + ...
        C.VP0.*exp(C.Ls/C.Rv.*(1.0./273.15-1.0./clim.T)) .* (clim.T<273.15);
clim.VP = clim.RH .* VPsat;
clim.q = clim.RH .* (VPsat .* C.eps ./ clim.Pres);

  
clim.VP = clim.RH .* VPsat;
clim.q = clim.RH .* (VPsat .* C.eps ./ clim.Pres);
% Air density (Dair)
clim.Dair = clim.Pres./C.Rd./clim.T;

% Time since last snow fall event (timelastsnow)
A.timelastsnow(clim.snow/(time.dt*24*3600)>C.Pthres) = time.TCUR;
if t==1
    A.timelastsnow(:) = time.TCUR; 
end

% Potential temperature (Theta)
clim.Theta = clim.T.*(C.Pref./clim.Pres).^(C.Rd/C.Cp);

A.climT = clim.T;
A.climP = clim.P;
A.climC = clim.C;
A.climRH = clim.RH;
A.climWS = clim.WS;
A.climPres = clim.Pres;
A.climsnow = clim.snow;
A.climrain = clim.rain;
A.climTlapse = clim.T_lapse;
A.climPreslapse = clim.Pres_lapse;
A.climPotlapse = clim.Theta_lapse;
A.climq = clim.q; 
A.climVP = clim.VP;

end