%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Load DEM and glacier outline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [grid] = func_init_grid(grid,io)

load(io.gridfile);

input = gridinfo;

input.y = input.y;
input.z = input.z;
input.mask = input.mask;
input.dhdt = input.dhdt; 

grid.x = input.x;
grid.y = input.y;
grid.z = input.z;
input.x = input.x;
grid.dhdt = input.dhdt; 

grid.Lx = size(grid.x,1);
grid.Ly = size(grid.x,2);

grid.maskfull = input.mask;

[~,FY] = gradient(grid.y);
if FY(1)<0
    grid.x = flipud(grid.x);
    grid.y = flipud(grid.y);    
    grid.z = flipud(grid.z);
    grid.maskfull = flipud(grid.maskfull);
    grid.dhdt=flipud(grid.dhdt); 
end
[FX,~] = gradient(grid.x);
if FX(1)<0
    grid.x = fliplr(grid.x);
    grid.y = fliplr(grid.y);    
    grid.z = fliplr(grid.z);
    grid.maskfull = fliplr(grid.maskfull);
    grid.dhdt=flipud(grid.dhdt); 
end

grid.mask = grid.maskfull(:);
grid.gpsum = sum(grid.mask>0);
grid.mask_short = grid.mask(grid.mask>0);

[grid.lat,grid.lon] = utm2ll(grid.x,grid.y,grid.utmzone);
grid.lat =  reshape(grid.lat,[numel(grid.x(1,:)),numel(grid.x(:,1))]); 
grid.lon =  reshape(grid.lon,[numel(grid.y(1,:)),numel(grid.y(:,1))]);

grid.x_mask = grid.x(grid.mask>0); 
grid.y_mask = grid.y(grid.mask>0);
grid.z_mask = grid.z(grid.mask>0);
grid.dhdt_mask = grid.dhdt(grid.mask>0);
grid.ind = find(grid.maskfull>0);
[grid.xind, grid.yind] = find(grid.maskfull>0);

%% Grid slope and aspect
spheroid = referenceEllipsoid('wgs84');
[ASPECT, SLOPE, gradN, gradE] = gradientm(grid.lat,grid.lon,grid.z,spheroid);

grid.slope = tan(SLOPE*pi/180);
grid.slope_x = gradE;
grid.slope_y = gradN;
grid.aspect = ASPECT;

grid.slope = grid.slope(grid.mask>0);
grid.slope_x = grid.slope_x(grid.mask>0);
grid.slope_y = grid.slope_y(grid.mask>0);
grid.aspect = grid.aspect(grid.mask>0);
grid.lat_mask = grid.lat(grid.mask>0);
grid.lon_mask = grid.lon(grid.mask>0);

grid.slope_beta = atan(grid.slope);
grid.slope_gamma = atan(-grid.slope_x./grid.slope_y).* (grid.slope_y>=0) ...
                    + (-pi + atan(-grid.slope_x./grid.slope_y)) .* (grid.slope_y<0 & grid.slope_x>0) ...
                    + (pi + atan(-grid.slope_x./grid.slope_y)) .* (grid.slope_y<0 & grid.slope_x<0);
grid.slope_gamma(grid.slope_x==0 & grid.slope_y<0) = pi;
grid.slope_gamma(grid.slope_x==0 & grid.slope_y==0) = 0;
grid.slope_gamma = -grid.slope_gamma; 

end

