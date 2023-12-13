function [x,y,z,vx,vy,vz] = targettrajectory_3d(t, tbirth, tdeath, waypoints)
% Generate a target trajectory that is composed of piece-wise CV,
% stop-and-go and turn motions. The movements are restricted to the
% horizontal plane, thus height is fixed.
%
% Input
%  t      : Time validation, [s], 1 x N
%  tbirth : Time at which target appears at the first way point, [s],
%           scalar
%  tdeath : Time at which target disappears, [s], scalar
%
%  
%  The first waypoint must be defined as:
%   waypoints(1).params.v : absolute value of velocity, [m/s]
%   waypoints(1).params.x     : x-coord, scalar, [m]
%   waypoints(1).params.y     : y-coord, scalar, [m]
%   waypoints(1).params.z     : z-coord, scalar, [m]
%
%  The last waypoint must be defined as:
%   waypoints(K).params.x     : x-coord, scalar, [m]
%   waypoints(K).params.y     : y-coord, scalar, [m]
%   waypoints(K).params.z     : z-coord, scalar, [m]
%
%  For the intermediate waypoints two types are allowed. The first type is
%  a turn:
%   waypoints(i).params.a     : absolute value of acceleration, [m/s/s]
%   waypoints(i).params.x     : x-coord, scalar, [m]
%   waypoints(i).params.y     : y-coord, scalar, [m]
%   waypoints(i).params.z     : z-coord, scalar, [m]
% 
%  The second type is a 'stop-and-go'. The target start with
%  decelerating at distance 'dr' from the previous waypoint until the
%  speed is zero. Then it pauses for 'dt' seconds and starts accelerating
%  in the direction of the next waypoint until it reaches speed 'v'.
%  Finally it moves with fixed speed 'v' until it reaches the waypoint. So:
%   waypoints(i).params.dr    : distance from the previous waypoint at
%                               which deceleration starts, [m]
%   waypoints(i).params.decel : absolute value of deceleration, [m/s/s]
%   waypoints(i).params.accel : absolute value of acceleration, [m/s/s]
%   waypoints(i).params.dt    : pause time, [s]
%   waypoints(i).params.v     : absolute value of velocity to which the
%                               target accelerates, [m/s] 
%
%              where i = 2, ..., K-1
%              and K-1 is the number segments
%
% Comments:
%           If tdeath is larger than the time at which the target passes
%           the last waypoint, then the trajectory after the last waypoint
%           is again CV with the speed vector unchanged.
%
%           If t is smaller than tbirth or larger than tdeath, then the
%           position and velocity are NaN.
%
%
% Test script : 
% =============
%
%  tbirth = 5;
%  tdeath = 20;
%  dt = 0.1;
%  t = [0:dt:tdeath];
%  waypoints(1).params.a = 0;
%  waypoints(1).params.x = 3;
%  waypoints(1).params.y = 2;
%  waypoints(1).params.z = 0;
%  waypoints(1).params.v = 1;
%  waypoints(2).params.dr = 1.5;
%  waypoints(2).params.dt = 1;
%  waypoints(2).params.v = 1.2;
%  waypoints(2).params.decel = 2;
%  waypoints(2).params.accel = 0.5;
%  waypoints(3).params.a = 1;
%  waypoints(3).params.x = 0;
%  waypoints(3).params.y = 6;
%  waypoints(3).params.z = 0;
%  waypoints(4).params.a = 10;
%  waypoints(4).params.x = 2;
%  waypoints(4).params.y = 8;
%  waypoints(4).params.z = 0;
%  waypoints(5).params.a = 4;
%  waypoints(5).params.x = -2;
%  waypoints(5).params.y = 9;
%  waypoints(5).params.z = 0;
%  waypoints(6).params.a = 200;
%  waypoints(6).params.x = -2;
%  waypoints(6).params.y = 2;
%  waypoints(6).params.z = 0;
%
%  count = 0;
%  x_wp = [];
%  y_wp = [];
%  z_wp = []
%  for i = 1:length(waypoints)
%      if isfield(waypoints(i).params,'x')
%          count = count + 1;
%          x_wp(count) = waypoints(i).params.x;
%          y_wp(count) = waypoints(i).params.y;
%          z_wp(count) = waypoints(i).params.z;
%      end
%  end
%  figure
%  [x,y,z,vx,vy,vz] = targettrajectory_3d(t, tbirth, tdeath, waypoints);
%  plot3(x,y,z,'.');
%  axis equal;
%  grid;
%  hold on;
%  plot3(x_wp,y_wp,z_wp,'o-');
%  for i = 1:length(x)
%      if ~isnan(x(i))
%           plot3([x(i) x(i)+dt*vx(i)],[y(i) y(i)+dt*vy(i)],[z(i) z(i)+dt*vz(i)],'r');
%      end
%  end
%  xlabel('x'); ylabel('y'); zlabel('z');
%  figure
%  plot(t,sqrt(vx.^2+vy.^2+vz.^2),'.');
%
%========================================================
% v.1.0 - 11.12.2012
%========================================================


% --- Number of time steps
N = length(t); % []

% --- Init the output
x  = nan*ones(size(t)); % [m],   1 x N
y  = nan*ones(size(t)); % [m],   1 x N
z  = nan*ones(size(t)); % [m],   1 x N
vx = nan*ones(size(t)); % [m/s], 1 x N
vy = nan*ones(size(t)); % [m/s], 1 x N
vz = nan*ones(size(t)); % [m/s], 1 x N

% --- Number of segments
Nseg = length(waypoints) - 1; % []

% --- Init
tseg_start = tbirth; % [s]
point1     = waypoints(1).params;
v          = waypoints(1).params.v; % [m/s]

% --- For each segment
for i = 1:Nseg
    
    if (isfield(waypoints(i+1).params,'a')&&(i~=Nseg))
        % If waypoints(i+1) is a horizontal turn and we are not in the
        % last segment
        
        % --- Get a position in the trajectory after the turn.
        %     If waypoint(i+2) is again a horizontal turn then this is just
        %     the position of the waypoint.
        %     If waypoint(i+2) is a 'stop-and-go' then we have to obtain
        %     the first next waypoint that is not a 'stop-and-go'.
        count = i+2;
        while ~isfield(waypoints(count).params,'x')
            count = count + 1;
        end
        point3 = waypoints(count).params;
        
        % --- Get parameters of the CV followed by turn segment
        [v_vec, rend_cv_vec, rend_turn_vec, r0_vec, dt_cv, dt_turn] = get_cv_traj_params(point1, waypoints(i+1).params, point3, v);
        
        % --- Valid time indices for CV segment
        ind = ((t >= tseg_start)&(t <= (tseg_start + dt_cv)));
        
        % --- Compute CV trajectory
        x(ind)  = point1.x + (t(ind)-tseg_start)*v_vec(1); % [m],   1 x Nind
        y(ind)  = point1.y + (t(ind)-tseg_start)*v_vec(2); % [m],   1 x Nind
        z(ind)  = point1.z + (t(ind)-tseg_start)*v_vec(3); % [m],   1 x Nind
        vx(ind) = v_vec(1);                                % [m/s], 1 x Nind
        vy(ind) = v_vec(2);                                % [m/s], 1 x Nind
        vz(ind) = v_vec(3);                                % [m/s], 1 x Nind
        
        % --- Valid time indices for turn segment
        ind = ((t > (tseg_start + dt_cv))&(t < (tseg_start + dt_cv + dt_turn)));
        
        % --- Compute turn trajectory
        dphi                                                = (t(ind) - (tseg_start + dt_cv))*((waypoints(i+1).params.a)/v); % [rad], 1 x Nind
        [x(ind), y(ind), z(ind), vx(ind), vy(ind), vz(ind)] = turn_3d(r0_vec, rend_cv_vec, v_vec, dphi);                     %        1 x Nind
        
        % --- Update start segment time and start point
        tseg_start = tseg_start + dt_cv + dt_turn; % [s]
        point1.x   = rend_turn_vec(1); % [m]
        point1.y   = rend_turn_vec(2); % [m]
        point1.z   = rend_turn_vec(3); % [m]

        
    elseif (isfield(waypoints(i+1).params,'dr')&&(i~=Nseg))
        % If waypoints(i+1) is a 'stop-and-go' and we are not in the
        % last segment
        
        % --- Get a position in the trajectory after the 'stop-and-go'.
        %     If waypoint(i+2) is a horizontal turn then this is just
        %     the position of the waypoint.
        %     If waypoint(i+2) is a 'stop-and-go' then we have to obtain
        %     the first next waypoint that is not a 'stop-and-go'.
        count = i+2;
        while ~isfield(waypoints(count).params,'x')
            count = count + 1;
        end
        point3 = waypoints(count).params;
        
        % --- Get parameters of the 'stop-and-go' which is a chain of:
        %     (1) a piece of CV of length 'dr'
        %     (2) a deceleration part
        %     (3) a pause of 'dt' [s]
        %     (4) an acceleration part
        [vx_cv, vy_cv, vz_cv, xend_cv, yend_cv, zend_cv, x_pause, y_pause, z_pause, xend_accel, yend_accel, zend_accel, dt_cv, dt_decel, dt_accel] = get_stop_and_go_traj_params(point1, waypoints(i+1).params, point3, v);
        
        % --- Valid time indices for CV segment
        ind = ((t >= tseg_start)&(t <= (tseg_start + dt_cv)));
        
        % --- Compute CV trajectory
        x(ind)  = point1.x + (t(ind)-tseg_start)*vx_cv; % [m],   1 x Nind
        y(ind)  = point1.y + (t(ind)-tseg_start)*vy_cv; % [m],   1 x Nind
        z(ind)  = point1.z + (t(ind)-tseg_start)*vz_cv; % [m],   1 x Nind
        vx(ind) = vx_cv;                                % [m/s], 1 x Nind
        vy(ind) = vy_cv;                                % [m/s], 1 x Nind
        vz(ind) = vz_cv;                                % [m/s], 1 x Nind

        % --- Valid time indices for deceleration segment
        ind = ((t > (tseg_start + dt_cv))&(t < (tseg_start + dt_cv + dt_decel)));
        
        % --- Compute deceleration trajectory
        decel = waypoints(i+1).params.decel; % [m/s/s]
        ax = -decel*vx_cv/v; % [m/s/s]
        ay = -decel*vy_cv/v; % [m/s/s]
        az = -decel*vz_cv/v; % [m/s/s]
        t_seg = t(ind) - tseg_start - dt_cv; % [s], 1 x Nind
        vx(ind) = ax*t_seg + vx_cv;                          % [m/s], 1 x Nind
        vy(ind) = ay*t_seg + vy_cv;                          % [m/s], 1 x Nind
        vz(ind) = az*t_seg + vz_cv;                          % [m/s], 1 x Nind
        x(ind)  = 0.5*ax*(t_seg.^2) + vx_cv*t_seg + xend_cv; % [m], 1 x Nind
        y(ind)  = 0.5*ay*(t_seg.^2) + vy_cv*t_seg + yend_cv; % [m], 1 x Nind
        z(ind)  = 0.5*az*(t_seg.^2) + vz_cv*t_seg + zend_cv; % [m], 1 x Nind
        
        % --- Valid time indices for pause segment
        dt = waypoints(i+1).params.dt; % [s]
        ind = ((t >= (tseg_start + dt_cv + dt_decel))&(t <= (tseg_start + dt_cv + dt_decel + dt)));
        
        % --- Compute pause 'trajectory'
        x(ind)  = x_pause; % [m],   1 x Nind
        y(ind)  = y_pause; % [m],   1 x Nind
        z(ind)  = z_pause; % [m],   1 x Nind
        vx(ind) = 0;       % [m/s], 1 x Nind
        vy(ind) = 0;       % [m/s], 1 x Nind 
        vz(ind) = 0;       % [m/s], 1 x Nind
        
        % --- Valid time indices for the acceleration trajectory
        ind = ((t >= (tseg_start + dt_cv + dt_decel + dt))&(t <= (tseg_start + dt_cv + dt_decel + dt + dt_accel)));
        
        % --- Compute acceleration trajectory
        accel = waypoints(i+1).params.accel; % [m/s/s]
        ax = accel*vx_cv/v; % [m/s/s]
        ay = accel*vy_cv/v; % [m/s/s]
        az = accel*vz_cv/v; % [m/s/s]
        t_seg = t(ind) - tseg_start - dt_cv - dt_decel - dt; % [s], 1 x Nind
        vx(ind) = ax*t_seg;                                  % [m/s], 1 x Nind
        vy(ind) = ay*t_seg;                                  % [m/s], 1 x Nind
        vz(ind) = az*t_seg;                                  % [m/s], 1 x Nind
        x(ind)  = 0.5*ax*(t_seg.^2) + x_pause;               % [m], 1 x Nind
        y(ind)  = 0.5*ay*(t_seg.^2) + y_pause;               % [m], 1 x Nind\
        z(ind)  = 0.5*az*(t_seg.^2) + z_pause;               % [m], 1 x Nind

        % --- Update start segment time, start point and speed
        tseg_start = tseg_start + dt_cv + dt_decel + dt + dt_accel; % [s]
        point1.x   = xend_accel;                                    % [m]
        point1.y   = yend_accel;                                    % [m]
        point1.z   = zend_accel;                                    % [m]
        v          = waypoints(i+1).params.v;                       % [m/s]
        
    elseif (i==Nseg) % We are in the last segment
               
        % --- Compute velocity vector
        dx    = waypoints(i+1).params.x - point1.x; % [m]
        dy    = waypoints(i+1).params.y - point1.y; % [m]
        dz    = waypoints(i+1).params.z - point1.z; % [m]
        r_cv  = sqrt(dx^2+dy^2+dz^2); % [m]
        vx_cv = v*dx/r_cv; % [m/s]
        vy_cv = v*dy/r_cv; % [m/s]
        vz_cv = v*dz/r_cv; % [m/s]
        %Tdur  = r_cv/v; % [s]
        
        % Indices of last CV piece
        %ind = ((t >= tseg_start)&(t <= (tseg_start + Tdur))); % 1 x Nind, []
        ind = ((t >= tseg_start)&(t <= tdeath)); % 1 x Nind, []
        
        % --- Compute CV trajectory
        x(ind)  = point1.x + (t(ind)-tseg_start)*vx_cv; % [m], 1 x Nind
        y(ind)  = point1.y + (t(ind)-tseg_start)*vy_cv; % [m], 1 x Nind
        z(ind)  = point1.z + (t(ind)-tseg_start)*vz_cv; % [m], 1 x Nind
        vx(ind) = vx_cv; % [m/s], 1 x Nind
        vy(ind) = vy_cv; % [m/s], 1 x Nind
        vz(ind) = vz_cv; % [m/s], 1 x Nind
    end
    
end

% --- Set target state to NaN for time > tdeath
ind = (t > tdeath);
x(ind)  = nan;
y(ind)  = nan;
z(ind)  = nan;
vx(ind) = nan;
vy(ind) = nan;
vz(ind) = nan;

function [v_vec, rend_cv_vec, rend_turn_vec, r0_vec, dt_cv, dt_turn] = get_cv_traj_params(point1,point2,point3,v)
% Get properties of the CV followed by turn segment
%
% v_vec         : Velocity vector in the cv segment, 3 x 1, [m/s]
% rend_cv_vec   : End point of cv segment, 3 x 1, [m]
% rend_turn_vec : End point of the turn segment, 3 x 1, [m]
% r0_vec        : Origin of the turn in the plane defined by point1,
%                 point2, point3, rend_cv_vec, rend_turn_vec, 3 x 1, [m]
% dt_cv         : Time spent in the cv segment, scalar, [s]
% dt_turn       : Time spent in the turn segment, scalar, [s]
%

% --- Vector d1
d1_vec = [point2.x - point1.x ; point2.y - point1.y ; point2.z - point1.z]; % 3 x 1, [m]

% --- Vector d2
d2_vec = [point3.x - point2.x ; point3.y - point2.y ; point3.z - point2.z]; % 3 x 1, [m]

% --- Distance between point 1 and point 2
d1 = norm(d1_vec); % scalar, [m]

% --- Distance between point 2 and point 3
d2 = norm(d2_vec); % scalar, [m]

% --- Velocity vector
v_vec = v*d1_vec/d1; % [m/s], 3 x 1

% --- Angle (without preserving sign) between d1_vec and d2_vec
cos_theta = dot(d1_vec, d2_vec)/(d1*d2);
theta     = acos(cos_theta); % [rad], [0 and pi]

% --- Radius of turn
radius = (v^2)/(point2.a); % [m]

% --- Length of cv segment
d    = radius*tan(theta/2); % [m]
r_cv = d1 - d; % [m]

% --- Check
if (r_cv <= 0)
    error('The radius of the bend is too large');
end

% --- Time duration of cv segment
dt_cv = r_cv/v; % [s]

% --- Time duration of turn segment
dt_turn = theta*radius/v; % [s]

% --- Compute the end point of the cv trajectory
rend_cv_vec = [point1.x ; point1.y ; point1.z] + dt_cv*v_vec; % [m/s], 3 x 1

% --- Normal vector (in the same direction as rotation axis of the turn)
n_vec = cross(d1_vec, d2_vec);  % [m], 3 x 1
n     = norm(n_vec);            % [m], scalar
n_vec = n_vec/n; % [], 3 x 1

% --- sin_theta
sin_theta = n/(d1*d2); % scalar, []

% --- Unit vector from 'rend_cv_vec' to the origin of the turn
n_v_vec = cross(n_vec, v_vec);    % [m/s], 3 x 1
n_v_vec = n_v_vec/norm(n_v_vec); % [],    3 x 1

% --- Get the origin of the turn in the plane defined by point1, point2 and point3
r0_vec = rend_cv_vec + radius*n_v_vec; % [m], 3 x 1

% --- End point of turn segment
rend_turn_vec = r0_vec - radius*n_v_vec*cos_theta + (radius*v_vec/v)*sin_theta; % [m], 3 x 1

function [vx_cv, vy_cv, vz_cv, xend_cv, yend_cv, zend_cv, x_pause, y_pause, z_pause, xend_accel, yend_accel, zend_accel, dt_cv, dt_decel, dt_accel] = get_stop_and_go_traj_params(point1, point2, point3, v)
% Get properties of the 'stop-and-go' segment
% Which is a chain of:
%     (1) a piece of CV of length 'dr'
%     (2) a deceleration part
%     (3) a pause of 'dt' [s]
%     (4) an acceleration part
%

% --- dx,dy,dz
dx = point3.x-point1.x; % [m]
dy = point3.y-point1.y; % [m]
dz = point3.z-point1.z; % [m]

% --- Distance between point 1 and point 3
r = sqrt(dx^2+dy^2+dz^2); % [m]

% --- Velocity components
vx_cv = v*dx/r; % [m/s]
vy_cv = v*dy/r; % [m/s]
vz_cv = v*dz/r; % [m/s]

% --- Time duration of cv segment
dt_cv = point2.dr/v; % [s]

% --- xend_cv, yend_cv, zend_cv
xend_cv = point1.x + dt_cv*vx_cv; % [m]
yend_cv = point1.y + dt_cv*vy_cv; % [m]
zend_cv = point1.z + dt_cv*vz_cv; % [m]

% --- Time duration of the deceleration
dt_decel = v/(point2.decel); % [s]

% --- Position where the target pauses
ax = -point2.decel*vx_cv/v; % [m/s/s]
ay = -point2.decel*vy_cv/v; % [m/s/s]
az = -point2.decel*vz_cv/v; % [m/s/s]
x_pause  = 0.5*ax*(dt_decel^2) + vx_cv*dt_decel + xend_cv; % [m]
y_pause  = 0.5*ay*(dt_decel^2) + vy_cv*dt_decel + yend_cv; % [m]
z_pause  = 0.5*az*(dt_decel^2) + vz_cv*dt_decel + zend_cv; % [m]

% --- Time duration of the acceleration
dt_accel = (point2.v)/(point2.accel); % [s]

% --- Position where the target comes out of the acceleration
ax = point2.accel*vx_cv/v;                  % [m/s/s]
ay = point2.accel*vy_cv/v;                  % [m/s/s]
az = point2.accel*vz_cv/v;                  % [m/s/s]
xend_accel = 0.5*ax*(dt_accel^2) + x_pause; % [m]
yend_accel = 0.5*ay*(dt_accel^2) + y_pause; % [m]
zend_accel = 0.5*az*(dt_accel^2) + z_pause; % [m]


% --- Check
if (((xend_accel-point1.x)^2 + (yend_accel-point1.y)^2 + (zend_accel-point1.z)^2) >= (r^2))
    error('The stop-and-go piece doesn''t fit');
end

