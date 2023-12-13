function [r,vd] = xyz2rvd(x,y,z,vx,vy,vz)
%
% The function xyz2rvd transforms the carthesian x, y, z
% coordinates to range coordinates
% *** Polar coordinates ***
%========================================================
% v.1.0 - 05.10.2012
%========================================================

r  = sqrt(x.^2 + y.^2 + z.^2); 	% range
vd = -(vx.*x + vy.*y + vz.*z)./r; 	% radial velocity with a minus sign

