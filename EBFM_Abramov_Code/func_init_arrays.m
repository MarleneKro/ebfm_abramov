%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set initial conditions from file or preset values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,clim,insol,OUT] = func_init_arrays(C,grid,io)

A = struct;
clim = struct;
insol = struct;

if (io.readbootfile)
    disp('Initialize from restart file...');
    cd(io.rebootdir);
    load(io.bootfilein);
    A = boot;
    A.snowmass(A.subSOIL(:,1)==1)               = 0.0;
    cd(io.homedir);
else
    disp('Initialize from manually set conditions...');
    A.subSS(1:grid.gpsum,1:grid.nl)             = 0;
    A.snowmass(1:grid.gpsum,1)                  = 0.0;              % snow mass (m w.e.)
    A.mbal_snow(1:grid.gpsum,1)                 = 0.0;              % cumulative snow mass balance (m w.e.)
    A.Tsurf(1:grid.gpsum,1)                     = 273.15;           % surface temperature (K)
    A.subT(1:grid.gpsum,1:grid.nl)              = 273.15;           % vertical temperatures (K)
    A.subW(1:grid.gpsum,1:grid.nl)              = 0.0;              % vertical irreducible water content (kg)
    A.subS(1:grid.gpsum,1:grid.nl)              = 0.0;              % vertical slush water content (kg)
    A.subSOIL(1:grid.gpsum,1:grid.nl)           = 0;                % vertical layer type (soil or ice/snow)
    A.subSOIL(grid.mask_short==2,1:grid.nl)     = 1;
    A.subD(1:grid.gpsum,1:grid.nl)              = C.Dice;           % vertical densities (kg m-3)
    A.subD(A.subSOIL==1)                        = C.soil_density;
    A.subTmean                                  = A.subT;           % vertical annual mean layer temperature (K)
    A.timelastsnow(1:grid.gpsum,1)              = 1.0;              % timestep of last snow fall (days)
    A.ys(1:grid.gpsum,1)                        = 1000;             % annual snow fall (mm w.e.)
    A.subZ(1:grid.gpsum,1:grid.nl)              = grid.max_subZ;
    A.alb_snow(1:grid.gpsum,1)                  = C.alb_fresh;      % snow albedo
    if grid.doubledepth                                             % vertical layer depths (m)
        A.subZ(grid.mask_short==1,1:grid.split(1)-1) = grid.max_subZ;
        if length(grid.split)>1
            for n=2:length(grid.split)
                A.subZ(grid.mask_short==1,grid.split(n-1):grid.split(n)-1) = (2.0^(n-1))*grid.max_subZ;
            end
        end
        A.subZ(grid.mask_short==1,grid.split(end):grid.nl) = (2.0^length(grid.split))*grid.max_subZ;
    end
    
end
A.mbal(1:grid.gpsum,1)                          = 0.0;              % cumulative climatic mass balance (m w.e.)
A.mbal_stake(1:grid.gpsum,1)                    = 0.0;              % cumulative stake mass balance (m w.e.)
A.smb(1:grid.gpsum,1)                           = 0.0;              % climatic mass balance (m w.e.)
A.smb_stake(1:grid.gpsum,1)                     = 0.0;              % stake mass balance (m w.e.)
A.subK(1:grid.gpsum,1:grid.nl)                  = 0.0;              % vertical conductivity (m2 s-1)
A.subCeff(1:grid.gpsum,1:grid.nl)               = 0.0;              % vertical effective heat capacity (J m-3 K)
A.subWvol(1:grid.gpsum,1:grid.nl)               = 0.0;              % vertical volumetric water content (fraction)
A.surfH(1:grid.gpsum,1)                         = 0.0;              % surface height
A.Dfreshsnow(1:grid.gpsum,1)                    = 0.0;              % fresh snow density
A.tstar(1:grid.gpsum,1)                         = 0.0;              % albedo decay timescale

A.runoff_irr_deep_mean(1:grid.gpsum,1)          = 0.0;              % runoff of irreducible water below base of the mdoel

clim.T(1:grid.gpsum,1)                          = 0.0;              % CLIM: temperature
clim.P(1:grid.gpsum,1)                          = 0.0;              % CLIM: precipitation
clim.snow(1:grid.gpsum,1)                       = 0.0;              % CLIM: snow precip
clim.rain(1:grid.gpsum)                         = 0.0;              % CLIM: rain precip
clim.yearsnow(1:grid.gpsum,1:grid.nl)           = 0.0;              % CLIM: annual snow precip
clim.logyearsnow(1:grid.gpsum,1:grid.nl)        = 0.0;              % CLIM: annual snow precip (log)
clim.C(1:grid.gpsum,1)                          = 0.0;              % CLIM: cloud cover
clim.WS(1:grid.gpsum,1)                         = 0.0;              % CLIM: wind speed
clim.RH(1:grid.gpsum,1)                         = 0.0;              % CLIM: relative humidity
clim.q(1:grid.gpsum,1)                          = 0.0;              % CLIM: specific humidity
clim.VP(1:grid.gpsum,1)                         = 0.0;              % CLIM: vapor pressure
clim.Dair(1:grid.gpsum,1)                       = 0.0;              % CLIM: air density
clim.Pres(1:grid.gpsum,1)                       = 0.0;              % CLIM: air pressure

OUT = struct;

end