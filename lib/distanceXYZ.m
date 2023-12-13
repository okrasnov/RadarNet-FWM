function out = distanceXYZ(XYZ,XYZ0)
%RANGE_XYZ calculate the distance between 3D points XYZ and XYZ0
%Usage:
     %out = distance(XYZ,XYZ0)
%Input:
%   XYZ and XYZ0 - coordinate vectors [x,y,z] of points
%Output:
%   out - distance in units of X,Y,Z
%==================================================================
%  v.1.0 - 27.12.2012, OK@TUD
%==================================================================


out=sqrt((XYZ(1)-XYZ0(1))^2+...
         (XYZ(2)-XYZ0(2))^2+...
         (XYZ(3)-XYZ0(3))^2);

end

