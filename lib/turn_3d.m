function [x,y,z,vx,vy,vz] = turn_3d(r0_vec, rstart_vec, vstart_vec, dphi)
% Constant speed 3D turn
%
% Input 
%
% r0_vec        : Center of the turn, 3 x 1, [m]
% rstart_vec    : Posistion vector a the start of the turn, 3 x 1, [m]
% vstart_vec    : Velocity vector at the start of the turn, 3 x 1, [m/s]
% dphi          : Angle of the turn (positive values),      1 x N, [rad]
%
% Output 
%  x,y,z,vx,vy,vz : Position and velocity coordinates after the turn,
%                   1 x N, [m] and [m/s]
%
%========================================================
% v.1.0 - 11.07.2012
%========================================================

% --- Check input
if any(dphi < 0)
    error('Angle of the turn must be positive values');
end

% --- Get vector from turn origin to the start of the turn
r_vec = rstart_vec - r0_vec; % 3 x 1, [m]

% --- Get radius of the turn
r = norm(r_vec); % scalar, [m]

% --- Get magnitude of the velocity
vstart = norm(vstart_vec); % scalar, [m/s]

% --- Compute the position after the turn
x = r0_vec(1) + r_vec(1)*cos(dphi) + (r/vstart)*vstart_vec(1)*sin(dphi); % [m], 1 x N
y = r0_vec(2) + r_vec(2)*cos(dphi) + (r/vstart)*vstart_vec(2)*sin(dphi); % [m], 1 x N
z = r0_vec(3) + r_vec(3)*cos(dphi) + (r/vstart)*vstart_vec(3)*sin(dphi); % [m], 1 x N

% --- Compute the velocity after the turn
vx = vstart_vec(1)*cos(dphi) - (vstart/r)*r_vec(1)*sin(dphi); % [m/s], 1 x N
vy = vstart_vec(2)*cos(dphi) - (vstart/r)*r_vec(2)*sin(dphi); % [m/s], 1 x N
vz = vstart_vec(3)*cos(dphi) - (vstart/r)*r_vec(3)*sin(dphi); % [m/s], 1 x N
