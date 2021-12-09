%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Update DEM (MK, 2021.02.03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [grid,A] = func_update_grid(io,grid,time,A)
if io.updateDEM == 1

  grid.z_mask = grid.z_mask  + grid.dhdt_mask*time.dt;
  grid.z = grid.z + grid.dhdt*time.dt;

  
%% Grid slope and aspect

spheroid = referenceEllipsoid('wgs84');
[ASPECT, SLOPE, gradN, gradE] = gradientm(grid.lat,grid.lon,grid.z,spheroid);

grid.slope = tan(SLOPE*pi/180);
grid.slope_x = gradE;
grid.slope_y = gradN;
grid.aspect = ASPECT;

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
end