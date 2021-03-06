%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute unattenuated incoming solar radiation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [insol] = func_flux_insol(grid,time,t,insol)

%% Time as a fraction of a year (in radians)
trad = 2*pi*((time.TCUR-datenum(year(datetime(datevec(time.TCUR))),1,1))/365.242);
    
%% Top-of-atmosphere radiation on a surface normal to incident beam
insol.I0 = 1353.0*(1.0+0.034*cos(trad));                                          
    
%% Solar declination
d = 0.322003-22.971*cos(trad)-0.357898*cos(2d0*trad) ...
            - 0.14398*cos(3d0*trad)+3.94638*sin(trad) ...
			+ 0.019334*sin(2d0*trad)+0.05928*sin(3d0*trad);

%% Surface slope
slopebeta = grid.slope_beta;
slopegamma = grid.slope_gamma; 

%% Solar hour angle
B = 360/365.*(time.TCUR-datenum(year(datetime(datevec(time.TCUR))),1,1)-81);
Tcor = 9.87.*sin(2.*B*pi/180)-7.53.*cos(B*pi/180)-1.5*sin(B*pi/180);
Tcor = Tcor/60/24;

h = -180.0+360.0*(hour(datetime(datevec(time.TCUR)))/24+Tcor+minute(datetime(datevec(time.TCUR)))/24/60);

%% Incoming solar radiation on flat surface
insol.TOAflat = insol.I0.*(sin(grid.lat_mask*pi/180.0).*sin(d*pi/180.0) + cos(grid.lat_mask*pi/180.0).*cos(d*pi/180.0).*cos(h*pi/180.0));

%% Incoming solar radiation on inclined surface
insol.TOA = insol.I0.*((sin(grid.lat_mask*pi/180.0).*cos(slopebeta) ...
        - cos(grid.lat_mask*pi/180.0).*sin(slopebeta).*cos(slopegamma)).*sin(d*pi/180.0) ...
        + (cos(grid.lat_mask*pi/180.0).*cos(slopebeta) + sin(grid.lat_mask*pi/180.0).*sin(slopebeta) ...
        .*cos(slopegamma)).*cos(d*pi/180.0).*cos(h*pi/180.0) ...
        + cos(d*pi/180.0).*sin(slopebeta).*sin(slopegamma).*sin(h*pi/180.0));
insol.TOA(insol.TOA<0) = 0.0; 


%% Shading by the surrounding topography
elevationangle_mask = asin(sin(grid.lat_mask.*pi/180.0).*sin(d*pi/180d0)+cos(grid.lat_mask.*pi/180.0) ...
                        .*cos(d*pi/180.0).*cos(h*pi/180.0));
if h<0
    azimuth = acos((cos(h*pi/180.0).*cos(d*pi/180.0).*sin(grid.lat_mask*pi/180.0) ...
                        - sin(d*pi/180.0).*cos(grid.lat_mask*pi/180.0))./cos(elevationangle_mask));
else
    azimuth = -acos((cos(h*pi/180.0).*cos(d*pi/180.0).*sin(grid.lat_mask*pi/180.0) ...
                        - sin(d*pi/180.0).*cos(grid.lat_mask*pi/180.0))./cos(elevationangle_mask));
end
az = real(azimuth);

yl = grid.Ly;
xl = grid.Lx;

ddx_mask(az<=-0.75*pi) = -tan(pi+az(az<=-0.75*pi));
ddx_mask(az<=-0.25*pi & az>-0.75*pi) = -1;
ddx_mask(az<=0.25*pi & az>-0.25*pi) = tan(az(az<=0.25*pi & az>-0.25*pi));
ddx_mask(az<=0.75*pi & az>0.25*pi) = 1;
ddx_mask(az>0.75*pi) = tan(pi-az(az>0.75*pi));

ddy_mask(az<=-0.75*pi) = 1;
ddy_mask(az<=-0.25*pi & az>-0.75*pi) = -tan(0.5*pi+az(az<=-0.25*pi & az>-0.75*pi));
ddy_mask(az<=0.25*pi & az>-0.25*pi) = -1;
ddy_mask(az<=0.75*pi & az>0.25*pi) = -tan(0.5*pi-az(az<=0.75*pi & az>0.25*pi));
ddy_mask(az>0.75*pi) = 1;

azmean = nanmean(az(:));

istart = 1*(azmean<=0) + yl*(azmean>0);
jstart = xl*(azmean<=-0.5*pi || azmean>0.5*pi) + 1*(azmean>-0.5*pi && azmean<=0.5*pi);
iend = yl*(azmean<=0) + 1*(azmean>0);
jend = 1*(azmean<=-0.5*pi || azmean>0.5*pi) + xl*(azmean>-0.5*pi && azmean<=0.5*pi);
ii = 1*(azmean<=0) - 1*(azmean>0);
jj = -1*(azmean<=-0.5*pi || azmean>0.5*pi) + 1*(azmean>-0.5*pi && azmean<=0.5*pi);

%%%%%%%%%%%%%
ddx = zeros(xl,yl);
ddy = zeros(xl,yl);
elevationangle = zeros(xl,yl);
for n=1:length(grid.xind)
    ddx(grid.xind(n),grid.yind(n)) = ddx_mask(n);
    ddy(grid.xind(n),grid.yind(n)) = ddy_mask(n);
    elevationangle(grid.xind(n),grid.yind(n)) = elevationangle_mask(n);
end
%%%%%%%%%%%%%

insol.shade = ones(xl,yl);
for i=istart:ii:iend
    for j=jstart:jj:jend
        if (grid.maskfull(j,i)>0)
            count=1;
            kk=round(i+ddx(j,i)*count);
            ll=round(j+ddy(j,i)*count);
            while ( kk~=0 && kk~=yl+1 && ll~=0 && ll~=xl+1)
                griddist = sqrt((grid.x(ll,kk)-grid.x(j,i))^2 + (grid.y(ll,kk)-grid.y(j,i))^2);
                gridangle = atan((grid.z(ll,kk)-grid.z(j,i))/griddist);
                if (elevationangle(j,i)<=gridangle)
                    insol.shade(j,i)=1;
                    kk=0;
                elseif (insol.shade(ll,kk)==0)
                    insol.shade(j,i)=0;
                    kk=0;
                else
                    insol.shade(j,i)=0;
                    count=count+1;
                    kk=round(i+ddx(j,i)*count);
                    ll=round(j+ddy(j,i)*count);
                end
                if(elevationangle(j,i)<0),insol.shade(j,i)=1; end
            end
        end
    end
end
insol.shade(:,1) = 0;
insol.shade(1,:) = 0;
insol.shade(:,grid.Ly) = 0;
insol.shade(grid.Lx,:) = 0;
insol.shade = reshape(insol.shade,[grid.Lx*grid.Ly,1]);
insol.shade = insol.shade(grid.mask>0);

end
