%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update density, temperature, water content of vertical model after:
%%% - snow fall & riming (+)
%%% - melting & sublimation (-)
%%% - gravitational densification
%%% - heat diffusion
%%% - liquid water percolation, refreezing & storage
%%% - slush runoff and storage
%%% - refreezing stored slush water
%%% - refreezing stored irreducible water
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A] = func_snowmodel(C,A,clim,dt,grid,time,phys)

nl = grid.nl;
A.sumSinit = sum(A.subS,2);
A.subTmean = A.subTmean .* (1-dt./C.yeardays) + A.subT .* dt./C.yeardays;

%% Fresh snow density [Kampenhout et al. 2017]
% Temperature-dependent contribution
A.Dfreshsnow_T(clim.T>C.T0+2) = 50+1.7*17^(3/2);
A.Dfreshsnow_T(clim.T<=C.T0+2 & clim.T>C.T0-15) = 50+1.7.*(clim.T(clim.T<=C.T0+2 & clim.T>C.T0-15)-C.T0+15).^(3/2);
A.Dfreshsnow_T(clim.T<=C.T0-15) = -3.8328*(clim.T(clim.T<=C.T0-15)-C.T0)-0.0333*(clim.T(clim.T<=C.T0-15)-C.T0).^2;

% Wind-dependent contribution
A.Dfreshsnow_W = 266.86.*(0.5.*(1+tanh(clim.WS./5))).^8.8;

A.Dfreshsnow = A.Dfreshsnow_T(:) + A.Dfreshsnow_W(:);

%% Update vertical layers after snow fall & riming
shift_snowfall = clim.snow.*C.Dwater./A.Dfreshsnow;
shift_riming = A.moist_deposition.*C.Dwater./A.Dfreshsnow;
shift_tot(:,1) = shift_snowfall + shift_riming;
A.surfH(:,1) = A.surfH(:,1) + shift_tot(:,1);
runoff_irr_deep = zeros(grid.gpsum,1);

while any(shift_tot>0)
    
    shift = min(shift_tot,grid.max_subZ);
    shift_tot = shift_tot - shift;
        
    subT_old = A.subT;
    subD_old = A.subD;
    subW_old = A.subW;
    subZ_old = A.subZ;
    subSOIL_old = A.subSOIL;
    subSS_old = A.subSS;

    i_noshift = find(subZ_old(:,1)+shift<=grid.max_subZ);
    i_shift = find(subZ_old(:,1)+shift>grid.max_subZ);

    A.subZ(i_noshift,1) = subZ_old(i_noshift,1) + shift(i_noshift,1);
    A.subT(i_noshift,1) = subT_old(i_noshift,1).*subZ_old(i_noshift,1)./A.subZ(i_noshift,1) + ...
        A.Tsurf(i_noshift,1).*shift(i_noshift,1)./A.subZ(i_noshift,1);
    A.subD(i_noshift,1) = subD_old(i_noshift,1).*subZ_old(i_noshift,1)./A.subZ(i_noshift,1) + ...
        A.Dfreshsnow(i_noshift,1).*shift(i_noshift,1)./A.subZ(i_noshift,1);
    A.subW(i_noshift,1) = subW_old(i_noshift,1);

    % shift layers down when first layer exceeds maximum thickness
    A.subZ(i_shift,3:nl) = subZ_old(i_shift,2:nl-1);
    A.subT(i_shift,3:nl) = subT_old(i_shift,2:nl-1);
    A.subD(i_shift,3:nl) = subD_old(i_shift,2:nl-1);
    A.subW(i_shift,3:nl) = subW_old(i_shift,2:nl-1);
    A.subZ(i_shift,2) = grid.max_subZ;
    A.subZ(i_shift,1) = (subZ_old(i_shift,1)+shift(i_shift,1)) - grid.max_subZ;
    A.subT(i_shift,2) = subT_old(i_shift,1).*subZ_old(i_shift,1)./A.subZ(i_shift,2) + ...
        A.Tsurf(i_shift,1).*(A.subZ(i_shift,2)-subZ_old(i_shift,1))./A.subZ(i_shift,2);
    A.subT(i_shift,1) = A.Tsurf(i_shift,1);
    A.subD(i_shift,2) = subD_old(i_shift,1).*subZ_old(i_shift,1)./A.subZ(i_shift,2) + ...
        A.Dfreshsnow(i_shift,1).*(A.subZ(i_shift,2)-subZ_old(i_shift,1))./A.subZ(i_shift,2);
    A.subD(i_shift,1) = A.Dfreshsnow(i_shift,1);
    A.subW(i_shift,2) = subW_old(i_shift,1);
    A.subW(i_shift,1) = 0.0;
    A.subSOIL(i_shift,2:nl) = subSOIL_old(i_shift,1:nl-1);  
    A.subSOIL(i_shift,1) = 0;
    A.subSS(i_shift,2:nl) = subSS_old(i_shift,1:nl-1);  
    A.subSS(i_shift,1) = 1;
    runoff_irr_deep(i_shift,1) = subW_old(i_shift,nl);
end

%% Update vertical layers after melt & sublimation
A.sumWinit = sum(A.subW,2);
mass_removed = (A.melt + A.moist_sublimation)*1d3;
mass_layer = A.subD.*A.subZ;
n = 0;
while any(mass_removed>0)
    n = n+1;
    cond1 = mass_removed(:,1)>mass_layer(:,n);
    cond2 = ~cond1 & mass_removed(:,1)>0;
    mass_removed(cond1,1) = mass_removed(cond1,1) - A.subD(cond1,n).*A.subZ(cond1,n);
    shift_tot(cond1,1) = shift_tot(cond1,1) - A.subZ(cond1,n);
    shift_tot(cond2,1) = shift_tot(cond2,1) - mass_removed(cond2,1)./mass_layer(cond2,n).*A.subZ(cond2,n);
    mass_removed(cond2,1) = 0.0;
end

while any(shift_tot<0)
    shift = max(shift_tot,-A.subZ(:,2));
    shift_tot = shift_tot - shift;
    
    shift(A.subSOIL(:,2)==1,1) = max(shift(A.subSOIL(:,2)==1,1),-A.subZ(A.subSOIL(:,2)==1,1));
    shift(A.subSOIL(:,1)==1,1) = 0;
    
    A.surfH(:,1) = A.surfH(:,1) + shift(:,1); 
    
    subT_old = A.subT;
    subD_old = A.subD;
    subW_old = A.subW;
    subZ_old = A.subZ;
    subSOIL_old = A.subSOIL;
    subSS_old = A.subSS;
    
    i_noshift = find(subZ_old(:,1)+shift>1d-17);
    i_shift = find(subZ_old(:,1)+shift<=1d-17);
    
    A.subZ(i_noshift,1) = subZ_old(i_noshift,1) + shift(i_noshift,1);
    A.subT(i_noshift,1) = subT_old(i_noshift,1);
    A.subD(i_noshift,1) = subD_old(i_noshift,1);
    temp = A.subZ(i_noshift,1)./subZ_old(i_noshift,1);
    A.subW(i_noshift,1) = subW_old(i_noshift,1).*temp;
    
    % shift layers up when first layer thickness is completely removed 
    A.subZ(i_shift,2:nl-1) = subZ_old(i_shift,3:nl);
    A.subT(i_shift,2:nl-1) = subT_old(i_shift,3:nl);
    A.subD(i_shift,2:nl-1) = subD_old(i_shift,3:nl);
    A.subW(i_shift,2:nl-1) = subW_old(i_shift,3:nl);
    A.subZ(i_shift,1) = subZ_old(i_shift,1) + subZ_old(i_shift,2) + shift(i_shift,1);
    A.subT(i_shift,1) = subT_old(i_shift,2);
    A.subD(i_shift,1) = subD_old(i_shift,2);
    temp = A.subZ(i_shift,1)./subZ_old(i_shift,2);
    A.subW(i_shift,1) = subW_old(i_shift,2).*temp;
    A.subSOIL(i_shift,1:nl-1) = subSOIL_old(i_shift,2:nl);
    A.subSS(i_shift,1:nl-1) = subSS_old(i_shift,2:nl);
    for n=1:length(i_shift)
        if grid.mask_short(i_shift(n))==1 && grid.doubledepth==1
            A.subZ(i_shift(n),nl) = 2.0^length(grid.split)*grid.max_subZ;
        else
            A.subZ(i_shift(n),nl) = grid.max_subZ;
        end
        if grid.mask_short(i_shift(n))==1
            A.subSOIL(i_shift(n),nl) = 0;
            A.subD(i_shift(n),nl) = subD_old(i_shift(n),nl);
        else
            A.subSOIL(i_shift(n),nl) = 1;
            A.subD(i_shift(n),nl) = C.soil_density;
        end
    end
    A.subT(i_shift,nl) =  subT_old(i_shift,nl);
    A.subW(i_shift,nl) = 0.0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Snow compaction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subD_old = A.subD;
subZ_old = A.subZ;
mliqmax = zeros(grid.gpsum,grid.nl);


%% Densification of firn [Ligtenberg et al. 2011]
cond = A.subD>=C.Dfirn & A.subSOIL==0;

grav_const = zeros(grid.gpsum,grid.nl);
grav_const(cond) =  0.07 .* max(1.435-0.151*clim.logyearsnow(cond),0.25) .* (A.subD(cond)<550) ...
    + 0.03 .* max(2.366-0.293*clim.logyearsnow(cond),0.25) .* (A.subD(cond)>=550);   

A.subD(cond) = A.subD(cond) + dt./C.yeardays.*grav_const(cond).*clim.yearsnow(cond).*C.g.*(C.Dice-A.subD(cond)).* ...
         exp(-C.Ec./(C.rd.*A.subT(cond)) + C.Eg./(C.rd.*A.subTmean(cond)));
     
A.Dens_firn = zeros(grid.gpsum,grid.nl);     
A.Dens_firn(cond) = dt./C.yeardays.*grav_const(cond).*clim.yearsnow(cond).*C.g.*(C.Dice-A.subD(cond)).* ...
         exp(-C.Ec./(C.rd.*A.subT(cond)) + C.Eg./(C.rd.*A.subTmean(cond)));
     
     
%% Densification of seasonal snow [Kampenhout et al. 2017]
cond = A.subD<C.Dfirn & A.subSOIL==0;

% Destructive metamorphism
CC_tap = 300; %default: 200
CC1 = zeros(grid.gpsum,grid.nl);
CC2 = zeros(grid.gpsum,grid.nl);
CC1(A.subD<CC_tap) = 1;
CC1(A.subD>=CC_tap) = exp(-0.046.*(A.subD(A.subD>=CC_tap)-CC_tap));
CC2(A.subW==0) = 1;
CC2(A.subW~=0) = 2;
CC3 = 2.777d-6;
CC4 = 0.04;

A.subD(cond) = A.subD(cond) + dt*3600*24.*A.subD(cond).*CC3.*CC2(cond).*CC1(cond).*exp(-CC4.*(C.T0-A.subT(cond)));
A.subD(cond) = min(A.subD(cond),C.Dice);

A.Dens_destr_metam = zeros(grid.gpsum,grid.nl);   
A.Dens_destr_metam(cond) = dt*3600*24.*A.subD(cond).*CC3.*CC2(cond).*CC1(cond).*exp(-CC4.*(C.T0-A.subT(cond)));

% Overburden pressure
CC5 = 0.1;
CC6 = 0.023;
CC7 = 4.0.*7.62237d6.*A.subD./358.0.*1./(1+60.*A.subW./(C.Dwater.*A.subZ));
Psload = zeros(grid.gpsum,grid.nl);
Visc = zeros(grid.gpsum,grid.nl);
temp = cumsum(A.subD.*A.subZ,2) - 0.5*A.subD.*A.subZ;
Psload(cond) = temp(cond); 
Visc(cond) = CC7(cond).*exp(CC5.*(C.T0-A.subT(cond))+CC6.*A.subD(cond));

A.subD(cond) = A.subD(cond) + dt*3600*24.*A.subD(cond).*Psload(cond)./Visc(cond);
A.subD(cond) = min(A.subD(cond),C.Dice);

A.Dens_overb_pres = zeros(grid.gpsum,grid.nl); 
A.Dens_overb_pres(cond) = dt*3600*24.*A.subD(cond).*Psload(cond)./Visc(cond);

% Drifting snow
MO = -0.069+0.66.*(1.25-0.0042.*(max(A.subD,50)-50));
SI = -2.868.*exp(-0.085.*repmat(clim.WS,1,grid.nl))+1+MO;
cond_SD = SI>0;
z_i = zeros(grid.gpsum,grid.nl);
z_i(:,2:grid.nl) = cumsum(A.subZ(:,1:grid.nl-1).*(3.25-SI(:,1:grid.nl-1)),2);
gamma_drift = max(0,SI.*exp(-z_i./0.1));
tau = 48*2*3600;
tau_i = tau./gamma_drift;

A.subD(cond & cond_SD) = A.subD(cond & cond_SD) + dt*3600*24.*max(350-A.subD(cond & cond_SD),0)./tau_i(cond & cond_SD);
A.subD(cond & cond_SD) = min(A.subD(cond & cond_SD),C.Dice);

A.Dens_drift = zeros(grid.gpsum,grid.nl); 
A.Dens_drift(cond & cond_SD) = dt*3600*24.*max(350-A.subD(cond & cond_SD),0)./tau_i(cond & cond_SD);

%% Update layer thickness, surface height and stored irreducible water after compaction
cond = A.subD<C.Dice;

A.subZ(cond) = subZ_old(cond).*subD_old(cond)./A.subD(cond);

mliqmax(cond) =  A.subD(cond).*A.subZ(cond).*0.0143.*exp(3.3.*(C.Dice-A.subD(cond))./C.Dice)./(1-0.0143.*exp(3.3.*(C.Dice-A.subD(cond))./C.Dice))...
    .*0.05.*min(C.Dice-A.subD(cond),20);
A.subW = min(mliqmax,A.subW);

shift(:,1) = sum(A.subZ,2)-sum(subZ_old,2);
A.surfH(:,1) = A.surfH(:,1) + shift(:,1);

runoff_irr = A.sumWinit - sum(A.subW,2);

%% Update subsurface temperatures after heat conduction / diffusion
dz1 = ((A.subZ(:,1)+0.5.*A.subZ(:,2)).^2);
dz2 = (0.5.*(A.subZ(:,3:grid.nl)+A.subZ(:,2:grid.nl-1)).^2);

cond = A.subSOIL==0;

kk = zeros(grid.gpsum,grid.nl);
c_eff = zeros(grid.gpsum,grid.nl);
kdTdz = zeros(grid.gpsum,grid.nl);

cond_temp = cond;
if any(any(cond_temp==1)); kk(cond) = 0.138-1.01d-3.*A.subD(cond)+3.233d-6.*A.subD(cond).^2; end
cond_temp = ~cond & A.subT<=C.T0;
if any(any(cond_temp==1)); kk(~cond & A.subT<=C.T0) = C.soil_Kfrozen; end
cond_temp = ~cond & A.subT>C.T0;
if any(any(cond_temp==1)); kk(~cond & A.subT>C.T0) = C.soil_Kthawed; end

cond_temp = cond;
if any(any(cond_temp==1)); c_eff(cond) = A.subD(cond) .* (152.2+7.122.*A.subT(cond)); end
cond_temp = ~cond & A.subT<=C.T0;
if any(any(cond_temp==1))
    c_eff(~cond & A.subT<=C.T0) = C.soil_Cfrozen + ...
        C.Lm.*1d3.*(C.soil_THwmax - C.soil_THwmin).*C.soil_delta./(A.subT(~cond & A.subT<=C.T0)-C.T0-C.soil_delta).^2;
end
cond_temp = ~cond & A.subT>C.T0;
if any(any(cond_temp==1)); c_eff(~cond & A.subT>C.T0) = C.soil_Cthawed; end
    
z_temp = A.subZ(:,2:end);
c_eff_temp = c_eff(:,2:end);
kk_temp = kk(:,2:end);
dt_stab = 0.5*min(c_eff_temp(:)).*min(z_temp(:).^2)./max(kk_temp(:))/3600.0/24.0;

tt = zeros(grid.gpsum,1);
while any(tt<dt)
    subT_old = A.subT;
    
    dt_temp = min(dt_stab,-tt+dt);

    tt = tt + dt_temp;
    
    cond_dt = dt_temp>0;
    
    kdTdz(cond_dt,2) = (kk(cond_dt,1).*A.subZ(cond_dt,1)+0.5.*kk(cond_dt,2).*A.subZ(cond_dt,2)).*(subT_old(cond_dt,2)-A.Tsurf(cond_dt,1)) ./ dz1(cond_dt,1);
    kdTdz(cond_dt,3:nl) = (kk(cond_dt,2:nl-1).*A.subZ(cond_dt,2:nl-1)+kk(cond_dt,3:nl).*A.subZ(cond_dt,3:nl)).*(subT_old(cond_dt,3:nl)-subT_old(cond_dt,2:nl-1)) ./ dz2(cond_dt,:);

    A.subT(cond_dt,2) = subT_old(cond_dt,2) + 24.*3600.*dt_temp(cond_dt,1).*(kdTdz(cond_dt,3)-kdTdz(cond_dt,2))./ ...
                        (c_eff(cond_dt,2) .* ...
                        (0.5*A.subZ(cond_dt,1)+0.5*A.subZ(cond_dt,2)+0.25*A.subZ(cond_dt,3)));
    A.subT(cond_dt,3:nl-1) = subT_old(cond_dt,3:nl-1) + 24.*3600.*repmat(dt_temp(cond_dt,1),1,grid.nl-3).*(kdTdz(cond_dt,4:nl)-kdTdz(cond_dt,3:nl-1))./ ...
                        (c_eff(cond_dt,3:nl-1) .* ...
                        (0.25*A.subZ(cond_dt,2:nl-2)+0.5*A.subZ(cond_dt,3:nl-1)+0.25*A.subZ(cond_dt,4:nl)));
                    
    A.subT(cond_dt,nl) = subT_old(cond_dt,nl) + 24.*3600.*dt_temp(cond_dt,1).*(C.geothermal_flux-kdTdz(cond_dt,nl))./ ...
                        (c_eff(cond_dt,nl) .* ...
                        (0.25*A.subZ(cond_dt,nl-1)+0.75*A.subZ(cond_dt,nl)));
end
A.subT(:,1) = A.Tsurf(:,1) + (A.subT(:,2)-A.Tsurf(:,1))./(A.subZ(:,1)+0.5*A.subZ(:,2)).*0.5.*A.subZ(:,1);
A.subT(A.subT>C.T0 & A.subSOIL==0) = C.T0;
A.subCeff = c_eff;
A.subK = kk;

%% Refreezing of percolating water and irreducible water storage
avail_W =   A.melt*1d3 + ...                                    % melt water
            clim.rain.*1d3 + ...                                % rain fall         
            (A.moist_condensation-A.moist_evaporation).*1d3;    % condensation or evaporation                                   
avail_W = max(avail_W,0);

subW_old = A.subW;
A.cpi = 152.2+7.122.*A.subT;
c1 = A.cpi.*A.subD.*A.subZ.*(C.T0-A.subT)./C.Lm;        % cold content limit on refreezing
c2 = A.subZ.*(1-A.subD./C.Dice).*C.Dice;                % density limit on refreezing
cond1 = c1 >= c2;
cond2 = ~cond1;
cond3 = A.subSOIL==1;
Wlim = zeros(grid.gpsum,nl);                            % maximum potential for refreezing
Wlim(cond1) = c2(cond1);
Wlim(cond2) = c1(cond2);
Wlim(cond3) = 0.0;
Wlim = max(Wlim,0);
mliqmax = zeros(grid.gpsum,grid.nl);                    % maximum irreducible water storage
noice = A.subD<C.Dice-1;
mliqmax(noice) =  A.subD(noice).*A.subZ(noice).*0.0143.*exp(3.3.*(C.Dice-A.subD(noice))./C.Dice)./(1-0.0143.*exp(3.3.*(C.Dice-A.subD(noice))./C.Dice))...
    .*0.05.*min(C.Dice-A.subD(noice),20);
Wirr = mliqmax-subW_old;

RP = zeros(grid.gpsum,grid.nl);                         % amount of refreezing of percolating water
leftW = zeros(grid.gpsum,1);
avail_W_loc = zeros(grid.gpsum,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z0 = C.perc_depth;
percolation = phys.percolation;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zz = cumsum(A.subZ,2)-0.5*A.subZ;

% description of water percolation
carrot = zeros(grid.gpsum,nl);      

% prob. density. distribution funtions
if  percolation == 1
    carrot(:,1) = 1;
elseif   percolation == 2      % normal law
    carrot = 2*exp( - zz.^2/2/(z0/3)^2 )/(z0/3)/sqrt(2*pi);
elseif   percolation == 3      % linear
    carrot = 2*(z0 - zz) / z0.^2;
    carrot = max(carrot,0);    % eliminate negative values
elseif   percolation == 4      % uniform
    [~, ind] = min(abs(zz - z0));
    carrot(1:ind) = 1/z0; clear ind
end

% probability function: take into account thickness of each layer
carrot = carrot.*A.subZ;

% normalize by sum to not loose or gain water
carrot = carrot./repmat(sum(carrot,2),1,grid.nl);

% multiply by available water mass to distribute it along the profile
carrot = carrot.*repmat(avail_W,1,grid.nl);

for n=1:nl
    avail_W_loc = avail_W_loc + carrot(:,n);
    
    % refreeze percolating water and store irreducible water
    cond1 = avail_W_loc>Wlim(:,n);
    
    RP(cond1,n) = Wlim(cond1,n);        % more water than refreezing limit
    leftW(cond1,1) = (avail_W_loc(cond1,1)-Wlim(cond1,n));
    A.subW(cond1,n) = subW_old(cond1,n) + min(leftW(cond1,1),Wirr(cond1,n));
    
    RP(~cond1,n) = avail_W_loc(~cond1); % less water than refreezing limit
    A.subW(~cond1,n) = subW_old(~cond1,n);
    
    avail_W_loc = avail_W_loc - RP(:,n) - (A.subW(:,n)-subW_old(:,n));
    
    % update temperature and density after refreezing
    A.subT(:,n) = A.subT(:,n) + C.Lm.*RP(:,n)./(A.subD(:,n).*A.cpi(:,n).*A.subZ(:,n));
    A.subD(:,n) = A.subD(:,n) + RP(:,n)./A.subZ(:,n);
end
avail_W = avail_W_loc;

%% Runoff of slush water
%  - according to linear-reservoir model
%  - recharged by nonrefrozen melt/rain water
%  - discharge depends on slush water amount and runoff time-scale
slush_old = sum(A.subS,2);
avail_S = 1.0/(1.0+dt/C.Trunoff).*slush_old;
avail_S(avail_S<1d-25) = 0.0;
runoff_slush = slush_old - avail_S;

%% Storage of slush water
avail_W = avail_W + avail_S;
A.subS = zeros(grid.gpsum,grid.nl);
slushspace = zeros(grid.gpsum,grid.nl);

for n=nl:-1:1
    if (any(A.subD(:,n)<C.Dice) && any(avail_W>0))
        
        % compute available space and fill with slush water
        % some space is already occupied by irreducible water
        slushspace(:,n) = max(A.subZ(:,n).*(1-A.subD(:,n)./C.Dice).*C.Dwater - A.subW(:,n),0);
        cond1 = avail_W>slushspace(:,n);
        A.subS(cond1,n) = slushspace(cond1,n);
        A.subS(~cond1,n) = avail_W(~cond1);
        avail_W = avail_W - A.subS(:,n);
    end
end
runoff_surface = avail_W;


%% Refreezing of slush water
A.cpi = 152.2+7.122.*A.subT;                            % specific heat capacity
c1 = A.cpi.*A.subD.*A.subZ.*(C.T0-A.subT)./C.Lm;        % cold content limit on refreezing
c2 = A.subZ.*(1-A.subD./C.Dice).*C.Dice;                % density limit on refreezing
cond1 = c1 >= c2;
cond2 = ~cond1;
cond3 = A.subSOIL==1;
Wlim = zeros(grid.gpsum,nl);                            % maximum potential for refreezing
Wlim(cond1) = c2(cond1);
Wlim(cond2) = c1(cond2);
Wlim(cond3) = 0.0;
RS = zeros(grid.gpsum,grid.nl);                         % amount of slush water refreezing

for n=nl:-1:1
    if (any(A.subS(:,n)>0) && any(A.subT(:,n)<C.T0))
        avail_W = A.subS(:,n);
        Wlim_loc = Wlim(:,n);

        cond = avail_W>Wlim_loc;
        RS(cond,n) = Wlim_loc(cond);
        RS(~cond,n) = avail_W(~cond);

        % update slush water content, temperature and density
        A.subS(:,n) = A.subS(:,n) - RS(:,n);
        A.subT(:,n) = A.subT(:,n) + C.Lm.*RS(:,n)./(A.subD(:,n).*A.cpi(:,n).*A.subZ(:,n));
        A.subD(:,n) = A.subD(:,n) + RS(:,n)./A.subZ(:,n);
    end
end

%% Refreezing of irreducible water
A.cpi = 152.2+7.122.*A.subT;                            % specific heat capacity
c1 = A.cpi.*A.subD.*A.subZ.*(C.T0-A.subT)./C.Lm;        % cold content limit on refreezing
c2 = A.subZ.*(1-A.subD./C.Dice).*C.Dice;                % density limit on refreezing
cond1 = c1 >= c2;
cond2 = ~cond1;
cond3 = A.subSOIL==1;
Wlim = zeros(grid.gpsum,nl);                            % maximum potential for refreezing                        
Wlim(cond1) = c2(cond1);
Wlim(cond2) = c1(cond2);
Wlim(cond3) = 0.0;
RI = zeros(grid.gpsum,grid.nl);                         % amount of irreducible water refreezing

for n=nl:-1:1
    if (any(A.subW(:,n)>0) && any(A.subT(:,n)<C.T0))
        
        % refreeze irreducible water
        avail_W = A.subW(:,n);
        cond1 = avail_W>Wlim(:,n);
        RI(cond1,n) = Wlim(cond1,n);
        RI(~cond1,n) = avail_W(~cond1);
        
        % update irreducible water content, tempeature and density
        A.subW(:,n) = A.subW(:,n) - RI(:,n);
        A.subT(:,n) = A.subT(:,n) + C.Lm.*RI(:,n)./(A.subD(:,n).*A.cpi(:,n).*A.subZ(:,n));
        A.subD(:,n) = A.subD(:,n) + RI(:,n)./A.subZ(:,n);
    end
end

A.refr = 1d-3.*(sum(RP,2) +  sum(RS,2) +  sum(RI,2));
A.refr_P = 1d-3.*sum(RP,2);
A.refr_S = 1d-3.*sum(RS,2);
A.refr_I = 1d-3.*sum(RI,2);
A.slushw = sum(A.subS,2);
A.irrw = sum(A.subW,2);

%% Reset summer surface
if day(datetime(datevec(time.TCUR)))==1 && month(datetime(datevec(time.TCUR)))==9
    A.subSS(:,:) = 0;
end
A.subSS(all(A.subSS==0,2),1) = 1;
A.refr_seasnl_snow = 1d-3.*sum((RP+RS+RI).*(A.subSS==1),2);
A.refr_intacc = 1d-3.*sum((RP+RS+RI).*(A.subSS==0),2);

%% Determine volumetric water content
cond = A.subSOIL==0;
A.subWvol(cond) = A.subW(cond).*1d-3./A.subZ(cond);
A.subWvol(~cond & A.subT<=C.T0) = C.soil_THwmin - (C.soil_THwmax - C.soil_THwmin).*C.soil_delta./(A.subT(~cond & A.subT<=C.T0)-C.T0-C.soil_delta);
A.subWvol(~cond & A.subT>C.T0) = C.soil_THwmax;

%% Merge / split layers to increase layer thickness at depths defined in grid.split
if grid.doubledepth
    for n=1:length(grid.split)
        % Merge layers (accumulation case)
        split = grid.split(n);
        cond = A.subZ(:,split)<=(2.0^(n-1))*grid.max_subZ & grid.mask_short==1;
        subT_old = A.subT;
        subD_old = A.subD;
        subW_old = A.subW;
        subZ_old = A.subZ;
        subS_old = A.subS;
        subSS_old = A.subSS;
        A.subZ(cond,split-1) = subZ_old(cond,split-1) + subZ_old(cond,split);
        A.subW(cond,split-1) = subW_old(cond,split-1) + subW_old(cond,split);
        A.subS(cond,split-1) = subS_old(cond,split-1) + subS_old(cond,split);
        A.subD(cond,split-1) = (subZ_old(cond,split-1).*subD_old(cond,split-1) + subZ_old(cond,split).*subD_old(cond,split))./(subZ_old(cond,split-1) + subZ_old(cond,split));
        A.subT(cond,split-1) = (subZ_old(cond,split-1).*subT_old(cond,split-1) + subZ_old(cond,split).*subT_old(cond,split))./(subZ_old(cond,split-1) + subZ_old(cond,split));
        A.subSS(cond,split-1) = max(subSS_old(cond,split-1),subSS_old(cond,split));
        A.subZ(cond,split:nl-1) = subZ_old(cond,split+1:nl);
        A.subW(cond,split:nl-1) = subW_old(cond,split+1:nl);
        A.subS(cond,split:nl-1) = subS_old(cond,split+1:nl);
        A.subD(cond,split:nl-1) = subD_old(cond,split+1:nl);
        A.subT(cond,split:nl-1) = subT_old(cond,split+1:nl);
        A.subSS(cond,split:nl-1) = subSS_old(cond,split+1:nl);
        A.subZ(cond,nl) = 2.0^length(grid.split)*grid.max_subZ;
        A.subT(cond,nl) = 2.0*subT_old(cond,nl) - subT_old(cond,nl-1);
        A.subD(cond,nl) = subD_old(cond,nl);
        A.subW(cond,nl) = 0.0;
        A.subS(cond,nl) = 0.0;
        A.subSS(cond,nl) = 0;
        
        % Split layers (ablation case)
        cond = A.subZ(:,split-2)>(2.0^(n-1))*grid.max_subZ & grid.mask_short==1;
        subT_old = A.subT;
        subD_old = A.subD;
        subW_old = A.subW;
        subZ_old = A.subZ;
        subS_old = A.subS;
        subSS_old = A.subSS;
        A.subZ(cond,split-2) = 0.5*subZ_old(cond,split-2);
        A.subW(cond,split-2) = 0.5*subW_old(cond,split-2);
        A.subS(cond,split-2) = 0.5*subS_old(cond,split-2);
        A.subSS(cond,split-2) = subSS_old(cond,split-2);
        A.subT(cond,split-2) = subT_old(cond,split-2);
        A.subD(cond,split-2) = subD_old(cond,split-2);
        A.subZ(cond,split-1) = 0.5*subZ_old(cond,split-2);
        A.subW(cond,split-1) = 0.5*subW_old(cond,split-2);
        A.subS(cond,split-1) = 0.5*subS_old(cond,split-2);
        A.subSS(cond,split-1) = subSS_old(cond,split-2);
        A.subT(cond,split-1) = subT_old(cond,split-2);
        A.subD(cond,split-1) = subD_old(cond,split-2);
        A.subZ(cond,split:nl) = subZ_old(cond,split-1:nl-1);
        A.subW(cond,split:nl) = subW_old(cond,split-1:nl-1);
        A.subS(cond,split:nl) = subS_old(cond,split-1:nl-1);
        A.subSS(cond,split:nl) = subSS_old(cond,split-1:nl-1);
        A.subT(cond,split:nl) = subT_old(cond,split-1:nl-1);
        A.subD(cond,split:nl) = subD_old(cond,split-1:nl-1);
        runoff_irr_deep(cond,1) = runoff_irr_deep(cond,1) + subW_old(cond,nl);
        runoff_slush(cond,1) = runoff_slush(cond,1) + subS_old(cond,nl);     
    end
end

%% Runoff
A.runoff_irr_deep_mean = A.runoff_irr_deep_mean .* (1-dt./C.yeardays) + runoff_irr_deep .* dt./C.yeardays;
A.runoff = 1d-3.*(runoff_surface + runoff_slush + runoff_irr + A.runoff_irr_deep_mean);
A.runoff_surf = 1d-3.*runoff_surface;
A.runoff_slush = 1d-3.*runoff_slush;
A.runoff_irr = 1d-3.*runoff_irr;
A.runoff_irr_deep = 1d-3.*A.runoff_irr_deep_mean;

end

