
%function RAEBELL_Sim

%========================================================
% v.1.0 - 27.12.2012, OK@TUD
%========================================================

%nondirective-radar_model
%
close all
clear all
clc

Rad=pi/180;
Deg=180/pi;

% To simulate or not the range-Doppler planes
Scena.flag_RD_plane=1;
%--- How sensors are working: -----------------------------
% 1 - simultaneously
% !=1 - consequently, 1st, then 2md, 3rd, so on...
% ---------------------------------------------------------------
Scena.Sensors.IsWorkingSimultaneously=1;

Scena.Start_Time=0;    % start simulation time [sec]
Scena.N_bursts = 10;%100;   % total number of simulated bursts per radar

Scena.N_Radars=4;
Scena.N_Times=Scena.N_bursts;
Scena.Time_Interval=1; %sec

%--- Limits for scene visualization -----------------------------
Scena.xlimits=[-8,10]*km;
Scena.ylimits=[-8,10]*km;
Scena.zlimits=[0,0]*km;


%% ---------------------------------------------------------------
% Highway(s) definition
% ---------------------------------------------------------------
Scena.N_Highways=1;

Scena.Highway(1).y_Distance=2*km; %from radar #1
Scena.Highway(1).Targets.Speed=4; %+ to the right, - to the left
Scena.Highway(1).Targets.Speed_Var=0;
Scena.Highway(1).Targets.Interval=1000;
Scena.Highway(1).Targets.RCS=5;

% Scena.Highway(2).y_Distance=2.01*km;
% Scena.Highway(2).Targets.Speed=20; %+ to the right, - to the left
% Scena.Highway(2).Targets.Speed_Var=0;
% Scena.Highway(2).Targets.Interval=100;
%
% Scena.Highway(3).y_Distance=2.03*km;
% Scena.Highway(3).Targets.Speed=-20; %+ to the right, - to the left
% Scena.Highway(3).Targets.Speed_Var=0;
% Scena.Highway(3).Targets.Interval=160;
%
% Scena.Highway(4).y_Distance=2.04*km;
% Scena.Highway(4).Targets.Speed=-10; %+ to the right, - to the left
% Scena.Highway(4).Targets.Speed_Var=0;
% Scena.Highway(4).Targets.Interval=100;



%% =====================================================================
% Setting radar nodes and characteristics
% =====================================================================
%  radar positions
Scena.Radar(1).position=[0,0,0]*km; %x0,y0,z0
Scena.Radar(1).beam.H_beamwidth=120*Rad;
Scena.Radar(1).beam.H_beam_direction=15*Rad;
Scena.Radar(1).Range.Max_Range=8*km;

Scena.Radar(2).position=[-1,5,0]*km; %x0,y0,z0
Scena.Radar(2).beam.H_beamwidth=170*Rad;
Scena.Radar(2).beam.H_beam_direction=160*Rad;
Scena.Radar(2).Range.Max_Range=8*km;

Scena.Radar(3).position=[5,1,0]*km; %x0,y0,z0
Scena.Radar(3).beam.H_beamwidth=170*Rad;
Scena.Radar(3).beam.H_beam_direction=-80*Rad;
Scena.Radar(3).Range.Max_Range=8*km;

Scena.Radar(4).position=[4,4,0]*km; %x0,y0,z0
Scena.Radar(4).beam.H_beamwidth=190*Rad;
Scena.Radar(4).beam.H_beam_direction=-120*Rad;
Scena.Radar(4).Range.Max_Range=8*km;

%  radar characteristics
for i_radar=1:Scena.N_Radars
    Scena.Radar(i_radar).sensorId=i_radar;
    Scena.Radar(i_radar).Carrier=1.3*GHz;         % [Hz]
    Scena.Radar(i_radar).Transmit_Power=1.0;      % [W]
    Scena.Radar(i_radar).Bandwidth=10.0*MHz;      % [Hz]
    Scena.Radar(i_radar).Pulse_Duration=0.436*msec;      % [sec]
    Scena.Radar(i_radar).Sampling_Duration=0.236*msec;   % [sec]
    Scena.Radar(i_radar).Tx_Gain   = 18; %[dB]
    Scena.Radar(i_radar).Rx_Gain   = 18; %[dB]
    Scena.Radar(i_radar).Noise_Fig = 3; %[dB]
    Scena.Radar(i_radar).antenna_coupling = 0.00001; % [times] = -50dB
    
    
    Scena.Radar(i_radar).PRF=2.2*kHz;             % [Hz]
    Scena.Radar(i_radar).Doppler.Bursts=256;      % number of sweeps in burst
    Scena.Radar(i_radar).Interval_Between_Bursts=4; % Time interval between bursts [sec]
    
    Scena.Radar(i_radar).Doppler.Umbiguity=Scena.Radar(i_radar).PRF*speed_of_light/(4*Scena.Radar(i_radar).Carrier);    % [m/sec]
    Scena.Radar(i_radar).Doppler.Resolution=Scena.Radar(i_radar).Doppler.Umbiguity/(Scena.Radar(i_radar).Doppler.Bursts/2); % [m/sec]
    
    Scena.Radar(i_radar).Range.Resolution=speed_of_light/(2*Scena.Radar(i_radar).Bandwidth*Scena.Radar(i_radar).Sampling_Duration/Scena.Radar(i_radar).Pulse_Duration);   % [m]
    Scena.Radar(i_radar).Range.Samples=1024;      % [m]
    Scena.Radar(i_radar).Range.Umbiguity=speed_of_light/(2*Scena.Radar(i_radar).PRF);  % [m]
    Scena.Radar(i_radar).Range.Output_Range_Samples=400;      % [m]
    
    % definition of radar beam polygon
    angles=[Scena.Radar(i_radar).beam.H_beam_direction-(Scena.Radar(i_radar).beam.H_beamwidth/2):0.01:Scena.Radar(i_radar).beam.H_beam_direction+(Scena.Radar(i_radar).beam.H_beamwidth/2)];
    Scena.Radar(i_radar).beam.x_circ=Scena.Radar(i_radar).position(1)+Scena.Radar(i_radar).Range.Max_Range.*sin(angles);
    Scena.Radar(i_radar).beam.y_circ=Scena.Radar(i_radar).position(2)+Scena.Radar(i_radar).Range.Max_Range.*cos(angles);
    Scena.Radar(i_radar).beam.x_circ(end+1)=Scena.Radar(i_radar).position(1);
    Scena.Radar(i_radar).beam.y_circ(end+1)=Scena.Radar(i_radar).position(2);
    
    % definition of radar beam direction line
    Scena.Radar(i_radar).beam.center_line=Coordinates_for_beam_direction(Scena.Radar(i_radar));
    
    % SNR [dB] after signal processing (so in r,vD domain) at 1 m for 1 m^2 rcs at boresight
    tempRadar.Pt = Scena.Radar(i_radar).Transmit_Power;
    tempRadar.frequency = Scena.Radar(i_radar).Carrier;
    tempRadar.B = Scena.Radar(i_radar).Bandwidth*...
        Scena.Radar(i_radar).Sampling_Duration/...
        Scena.Radar(i_radar).Pulse_Duration;
    tempRadar.lin_Gt = 10^(Scena.Radar(i_radar).Tx_Gain/10);
    tempRadar.lin_Gr = 10^(Scena.Radar(i_radar).Rx_Gain/10);
    tempRadar.lin_F  = 10^(Scena.Radar(i_radar).Noise_Fig/10);
    tempRadar.lin_PG = tempRadar.B*...
        Scena.Radar(i_radar).Sampling_Duration*...
        Scena.Radar(i_radar).Doppler.Bursts;
    range=1;
    rcs=1;
    T=300;
    Scena.Radar(i_radar).lin_SNR=radar_SNR(tempRadar,range,rcs,T);
    
    sensor(i_radar).x  = Scena.Radar(i_radar).position(1);                   % Sensor x position [m]
    sensor(i_radar).y  = Scena.Radar(i_radar).position(2);                  	% Sensor y position [m]
    sensor(i_radar).z  = Scena.Radar(i_radar).position(3);                  	% Sensor z position [m]
    sensor(i_radar).tr  = Scena.Radar(i_radar).beam.H_beam_direction;                 	% Training [rad]
    sensor(i_radar).tlt  = 0;                	% Tilt, [rad]
    
    
end


%% ------------------------------------------------------------------------
% Highway targets initialization
% ------------------------------------------------------------------------
disp(['i_radar=',num2str(i_radar)])

if (Scena.N_Highways>0)
    angles_1=Scena.Radar(1).beam.H_beam_direction-(Scena.Radar(1).beam.H_beamwidth/2);
    angles_end=Scena.Radar(1).beam.H_beam_direction+(Scena.Radar(1).beam.H_beamwidth/2);
    for ind=1:Scena.N_Highways
        Scena.Highway(ind).x_min=Scena.Radar(1).position(1)+Scena.Highway(ind).y_Distance*tan(angles_1);
        Scena.Highway(ind).x_max=Scena.Radar(1).position(1)+Scena.Highway(ind).y_Distance*tan(angles_end);
        Scena.Highway(ind).y=Scena.Radar(1).position(2)+Scena.Highway(ind).y_Distance;
        Scena.Highway(ind).z=0;
    end;
    
    for ind=1:Scena.N_Highways
        Road(ind)=0;
        if (Scena.Highway(ind).x_min<Scena.Highway(ind).x_max)
            Road(ind)=Scena.Highway(ind).x_max-Scena.Highway(ind).x_min;
            Scena.Highway(ind).N_Targets=floor(Road(ind)/Scena.Highway(ind).Targets.Interval)+1;
            for i=1:Scena.Highway(ind).N_Targets
                Scena.Highway(ind).Target(i).position=[Scena.Highway(ind).x_min+(i-1)*Scena.Highway(ind).Targets.Interval,Scena.Highway(ind).y,Scena.Highway(ind).z];
                Scena.Highway(ind).Target(i).velocity=Scena.Highway(ind).Targets.Speed;
                if (Scena.Highway(ind).Targets.Speed_Var~=0)
                    Scena.Highway(ind).Target(i).velocity=Scena.Highway(ind).Target(i).velocity+Scena.Highway(ind).Targets.Speed_Var*randn(1,1);
                end
            end;
        end;
    end;
end;

%% ------------------------------------------------------------------------
% waupoint target model
% ------------------------------------------------------------------------
tbirth = Scena.Start_Time;
if (Scena.Sensors.IsWorkingSimultaneously==1)
    tdeath = tbirth +(0.5*Scena.Radar(1).Doppler.Bursts/Scena.Radar(1).PRF+Scena.Radar(i_radar).Interval_Between_Bursts)*Scena.N_bursts;
else
    tdeath = tbirth +(0.5*Scena.Radar(1).Doppler.Bursts/Scena.Radar(1).PRF+Scena.Radar(i_radar).Interval_Between_Bursts)*Scena.N_bursts*Scena.N_Radars;
end
dt = (tdeath-tdeath)/1000;
if (dt<0.1)
    dt=0.1;
end
t = [tbirth:dt:tdeath];

waypoints(1).params.a = 0;
waypoints(1).params.x = 4000;
waypoints(1).params.y = -2000;
waypoints(1).params.z = 500;
waypoints(1).params.v = 70;

waypoints(2).params.a = 10;
waypoints(2).params.x = 1000;
waypoints(2).params.y = 2000;
waypoints(2).params.z = 500;

waypoints(3).params.a = 10;
waypoints(3).params.x = 3000;
waypoints(3).params.y = 4000;
waypoints(3).params.z = 500;

waypoints(4).params.a = 10;
waypoints(4).params.x = 0;
waypoints(4).params.y = 4000;
waypoints(4).params.z = 500;

waypoints(5).params.a = 10;
waypoints(5).params.x = -6000;
waypoints(5).params.y = -2000;
waypoints(5).params.z = 500;

waypoints(6).params.a = 10;
waypoints(6).params.x = 4000;
waypoints(6).params.y = 4000;
waypoints(6).params.z = 500;

waypoints(7).params.a = 10;
waypoints(7).params.x = 0;
waypoints(7).params.y = 8000;
waypoints(7).params.z = 500;

[Scena.ufo.x,Scena.ufo.y,Scena.ufo.z,Scena.ufo.vx,Scena.ufo.vy,Scena.ufo.vz] = targettrajectory_3d(t, tbirth, tdeath, waypoints);
Scena.ufo.t=t;
Scena.ufo.RCS=7; %[m^2]

%% ------------------------------------------------------------------------
% Bursts simulator
% ------------------------------------------------------------------------

Time1=Scena.Start_Time;
burst_counter=0;
for i_Time=1:Scena.N_bursts % loop for total number of simulated bursts per radar
    
    if (Scena.Sensors.IsWorkingSimultaneously==1)
        Time=Time1+0.5*Scena.Radar(i_radar).Doppler.Bursts/Scena.Radar(i_radar).PRF; %Mid burst time [s]
        Time1=Time+Scena.Radar(i_radar).Interval_Between_Bursts; % Time interval between bursts [sec]
    end
    for i_radar=1:Scena.N_Radars % loop for total number of radars
        if (Scena.Sensors.IsWorkingSimultaneously~=1)
            Time=Time1+0.5*Scena.Radar(i_radar).Doppler.Bursts/Scena.Radar(i_radar).PRF; %Mid burst time [s]
            Time1=Time+Scena.Radar(i_radar).Interval_Between_Bursts; % Time interval between bursts [sec]
        end
        burst_counter=burst_counter+1
        
        % Mid burst time [s]
        burstData(burst_counter).t        = Time;
        % Burst identification number, []
        burstData(burst_counter).burst    = burst_counter;
        % Sensor identification number, []
        burstData(burst_counter).sensorId = Scena.Radar(i_radar).sensorId;
        % Range in the center of first range quant, [m]
        burstData(burst_counter).rStart   = Scena.Radar(i_radar).Range.Resolution/2;
        % Number of range quants, []
        burstData(burst_counter).nR       = Scena.Radar(i_radar).Range.Output_Range_Samples;
        % Range bin size, [m]
        burstData(burst_counter).dr       = Scena.Radar(i_radar).Range.Resolution;
        % Doppler speed in center of the first Doppler bin, [m/s]
        burstData(burst_counter).vdStart  = -Scena.Radar(i_radar).Doppler.Umbiguity+Scena.Radar(i_radar).Doppler.Resolution/2;
        % Number of Doppler bins, []
        burstData(burst_counter).nV       = Scena.Radar(i_radar).Doppler.Bursts;
        % Doppler bin size, [m/s]
        burstData(burst_counter).dv       = Scena.Radar(i_radar).Doppler.Resolution;
        % Blind speed, [m/s]
        burstData(burst_counter).vBlind   = 2*Scena.Radar(i_radar).Doppler.Umbiguity;
        
        %======================================
        % TBD
        %======================================
        % -3 dB pulse width [m]
        burstData(burst_counter).pulseWidth = NaN;
        % -3 dB Doppler filter width [m/s]
        burstData(burst_counter).dopplerWidth = NaN;
        % SNR [dB] after signal processing (so in r,vD domain)
        % at 1 m for 1 m^2 rcs at boresight
        % SNR = Pt x Gt x Gr x Lambda^2 x RCS / (4pi)^3 x R^4 x N
        % N = k x T x B x F
        burstData(burst_counter).snrUnitDb = 10*log10(Scena.Radar(i_radar).lin_SNR);
        
        
        
        
        Scena.Radar(i_radar).Time(i_Time).Time=Time;
        
        i_N_targets=0;
        
        if Scena.N_Highways>0
            for ind=1:Scena.N_Highways
                for i=1:Scena.Highway(ind).N_Targets
                    i_N_targets=i_N_targets+1;
                    
                    burstData(burst_counter).target(i_N_targets).m = 1;     	% [0,1] Indicates if the target is present
                    % 1: present, 0: not present, nTargets x 1
                    
                    position(1)=Scena.Highway(ind).Target(i).position(1)+Scena.Highway(ind).Target(i).velocity*Time;
                    position(2)=Scena.Highway(ind).Target(i).position(2);
                    position(3)=Scena.Highway(ind).Target(i).position(3);
                    
                    burstData(burst_counter).target(i_N_targets).x  = position(1);    	% [m],   nTargets x 1
                    burstData(burst_counter).target(i_N_targets).vx = Scena.Highway(ind).Target(i).velocity;   	% [m/s], nTargets x 1
                    burstData(burst_counter).target(i_N_targets).y  = position(2);      	% [m],   nTargets x 1
                    burstData(burst_counter).target(i_N_targets).vy = 0;    	% [m/s], nTargets x 1
                    burstData(burst_counter).target(i_N_targets).z  = position(3);      	% [m],   nTargets x 1
                    burstData(burst_counter).target(i_N_targets).vz = 0;    	% [m/s], nTargets x 1
                    
                    Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Range=distanceXYZ(Scena.Radar(i_radar).position,position);
                    burstData(burst_counter).target(i_N_targets).r = Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Range;     	% [m],   nTargets x 1
                    
                    
                    Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Amplitude=(Scena.Highway(1).y_Distance/Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Range)^4;
                    
                    cosA_RTx=(position(1)-Scena.Radar(i_radar).position(1))/Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Range;
                    cosA_RTy=(position(2)-Scena.Radar(i_radar).position(2))/Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Range;
                    cosA_RTz=(position(3)-Scena.Radar(i_radar).position(3))/Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Range;
                    
                    cosA_TVx=1;
                    cosA_TVy=0;
                    cosA_TVz=0;
                    gamma=cosA_RTx*cosA_TVx+cosA_RTy*cosA_TVy+cosA_RTz*cosA_TVz;
                    %gamma0=gamma*Deg
                    Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Doppler=Scena.Highway(ind).Target(i).velocity*(gamma);
                    
                    burstData(burst_counter).target(i_N_targets).vd = Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Doppler;    	% [m/s], nTargets x 1
                    
                    dx=burstData(burst_counter).target(i_N_targets).x-Scena.Radar(i_radar).position(1);
                    dy=burstData(burst_counter).target(i_N_targets).y-Scena.Radar(i_radar).position(2);
                    direction_angle=atan2(dx,dy);
                    
                    burstData(burst_counter).target(i_N_targets).az = direction_angle-Scena.Radar(i_radar).beam.H_beam_direction;    	% [rad], nTargets x 1
                    
                    if (abs(burstData(burst_counter).target(i_N_targets).az) > (Scena.Radar(i_radar).beam.H_beamwidth/2))
                        burstData(burst_counter).target(i_N_targets).m = 0;
                        Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Visibility=0;
                    else
                        Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Visibility=1;
                    end;
                    RCS=Scena.Highway(ind).Targets.RCS; % here can be added burst-to-burst fluctuations of RCS
                    burstData(burst_counter).target(i_N_targets).snrLin = Scena.Radar(i_radar).lin_SNR*RCS/burstData(burst_counter).target(i_N_targets).r^4;	% [],    nTargets x 1
                    
                end;
            end;
        end; % Scena.N_Highways>0
        
        %% ------------------------------------------------------------------------
        % waupoint target model
        % ------------------------------------------------------------------------
        i_N_targets=i_N_targets+1;
        if (tbirth<Time)&&(tdeath>Time)
            burstData(burst_counter).target(i_N_targets).m = 1;     	% [0,1] Indicates if the target is present
            % 1: present, 0: not present, nTargets x 1
            j_min=find(abs(Scena.ufo.t-Time)==min(abs(Scena.ufo.t-Time)));
            %[j_min,Time]
            
            burstData(burst_counter).target(i_N_targets).x  = Scena.ufo.x(j_min);    	% [m],   nTargets x 1
            burstData(burst_counter).target(i_N_targets).vx = Scena.ufo.vx(j_min);   	% [m/s], nTargets x 1
            burstData(burst_counter).target(i_N_targets).y  = Scena.ufo.y(j_min);      	% [m],   nTargets x 1
            burstData(burst_counter).target(i_N_targets).vy = Scena.ufo.vy(j_min);    	% [m/s], nTargets x 1
            burstData(burst_counter).target(i_N_targets).z  = Scena.ufo.z(j_min);      	% [m],   nTargets x 1
            burstData(burst_counter).target(i_N_targets).vz = Scena.ufo.vz(j_min);    	% [m/s], nTargets x 1
            
            [Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Range,...
                Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Doppler] = xyz2rvd(...
                burstData(burst_counter).target(i_N_targets).x-Scena.Radar(i_radar).position(1),...
                burstData(burst_counter).target(i_N_targets).y-Scena.Radar(i_radar).position(2),...
                burstData(burst_counter).target(i_N_targets).z-Scena.Radar(i_radar).position(3),...
                burstData(burst_counter).target(i_N_targets).vx,...
                burstData(burst_counter).target(i_N_targets).vy,...
                burstData(burst_counter).target(i_N_targets).vz);
            
            burstData(burst_counter).target(i_N_targets).r = Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Range;     	% [m],   nTargets x 1
            burstData(burst_counter).target(i_N_targets).vd = Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Doppler;    	% [m/s], nTargets x 1
            
            Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Amplitude=(Scena.Highway(1).y_Distance/Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Range)^4;
            
            dx=burstData(burst_counter).target(i_N_targets).x-Scena.Radar(i_radar).position(1);
            dy=burstData(burst_counter).target(i_N_targets).y-Scena.Radar(i_radar).position(2);
            direction_angle=atan2(dx,dy);
            burstData(burst_counter).target(i_N_targets).az = direction_angle-Scena.Radar(i_radar).beam.H_beam_direction;    	% [rad], nTargets x 1
            
            if (abs(burstData(burst_counter).target(i_N_targets).az) > (Scena.Radar(i_radar).beam.H_beamwidth/2))
                burstData(burst_counter).target(i_N_targets).m = 0;
                Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Visibility=0;
            else
                Scena.Radar(i_radar).Time(i_Time).Target(i_N_targets).Visibility=1;
            end;
            
            RCS=Scena.ufo.RCS; % here can be added burst-to-burst fluctuations of RCS
            burstData(burst_counter).target(i_N_targets).snrLin = Scena.Radar(i_radar).lin_SNR*RCS/burstData(burst_counter).target(i_N_targets).r^4;	% [],    nTargets x 1
        else
            burstData(burst_counter).target(i_N_targets).m = 0;     	% [0,1] Indicates if the target is present
        end;
        
        Scena.Radar(i_radar).Time(i_Time).N_targets=i_N_targets;
        
        %% ------------------------------------------------------------------------
        % Targets in radar's range-doppler plane/space,
        % taken into accounts the range and doppler ambiguities
        % ------------------------------------------------------------------------
        i_target=0;
        for i_ind=1:Scena.Radar(i_radar).Time(i_Time).N_targets
            i_target=i_target+1;
            Target(i_target).Visibility=Scena.Radar(i_radar).Time(i_Time).Target(i_ind).Visibility;
            Target(i_target).range=mod(Scena.Radar(i_radar).Time(i_Time).Target(i_ind).Range,Scena.Radar(i_radar).Range.Umbiguity);
            Target(i_target).Doppler=mod(Scena.Radar(i_radar).Time(i_Time).Target(i_ind).Doppler+Scena.Radar(i_radar).Doppler.Umbiguity,2*Scena.Radar(i_radar).Doppler.Umbiguity)-Scena.Radar(i_radar).Doppler.Umbiguity;
            
            Target(i_target).amplitude=sqrt(burstData(burst_counter).target(i_N_targets).snrLin);
        end
        N_targets=i_target;
        
        %% ------------------------------------------------------------------------
        % Adding targets with Zero-Doppler
        % ------------------------------------------------------------------------
        %        RCS_zero_doppler=0.1.*log10([1:Scena.Radar(i_radar).Range.Samples])+3.0.*randn(1,Scena.Radar(i_radar).Range.Samples);
        RCS_zero_doppler=0.1.*[1:Scena.Radar(i_radar).Range.Samples].^(4/6).*(1+3.0.*randn(1,Scena.Radar(i_radar).Range.Samples));
        %        RCS_zero_doppler=ones([1,Scena.Radar(i_radar).Range.Samples])*0.1;
        
        
        for index=1:Scena.Radar(i_radar).Range.Samples
            N_targets=N_targets+1;
            Target(N_targets).Visibility=1;
            Target(N_targets).Doppler=0;
            if (i==1)
                % direct coupling target
                Target(N_targets).range=0;
                Target(N_targets).amplitude=sqrt(Scena.Radar(i_radar).lin_SNR*Scena.Radar(i_radar).antenna_coupling);
            else
                Target(N_targets).range=0+Scena.Radar(i_radar).Range.Resolution*index;
                Target(N_targets).amplitude = sqrt(Scena.Radar(i_radar).lin_SNR*RCS_zero_doppler(index)/Target(N_targets).range^4);
            end
        end
        
        %% ------------------------------------------------------------------------
        % Range-Doppler planes
        % ------------------------------------------------------------------------
        if (Scena.flag_RD_plane==1)
            
            N_R_FFT=2*Scena.Radar(i_radar).Range.Samples;
            t_R_FFT=0:(N_R_FFT-1);
            t_R_FFT=(t_R_FFT-length(t_R_FFT)/2);
            w_R=hamming(length(t_R_FFT));
            w_R=w_R(:)';
            N_D_FFT=Scena.Radar(i_radar).Doppler.Bursts;
            t_D_FFT=0:(N_D_FFT-1);
            w_D=hamming(length(t_D_FFT));
            w_D=w_D(:)';
            
            an_R=0;%0.1;
            an_D=0;%0.1;
            
            RD_plane=1+1*randn(Scena.Radar(i_radar).Range.Output_Range_Samples,Scena.Radar(i_radar).Doppler.Bursts);
            
            for i_target=1:N_targets
                if (Target(i_target).Visibility==1)
                    omega0=pi*(Target(i_target).range/(Scena.Radar(i_radar).Range.Resolution*(Scena.Radar(i_radar).Range.Samples)));
                    x_R=cos(omega0.*t_R_FFT);
                    x1_R=x_R.*ones(size(t_R_FFT))+an_R.*randn(size(t_R_FFT));
                    x1_R=Target(i_target).amplitude.*x1_R;
                    s=fft(w_R.*x1_R,length(x_R));
                    s_R=s(1:Scena.Radar(i_radar).Range.Output_Range_Samples)./(Scena.Radar(i_radar).Range.Samples/2);
                    
                    omega1=pi*(Target(i_target).Doppler/Scena.Radar(i_radar).Doppler.Umbiguity);
                    omega1=j*omega1;
                    x_D=exp(omega1.*t_D_FFT);
                    x1_D=x_D.*ones(size(t_D_FFT))+an_D.*randn(size(t_D_FFT));
                    s_D=fft(w_D.*x1_D,length(x_D));
                    s_D=fftshift(s_D);
                    s_D=s_D./N_D_FFT;
                    
                    RD_plane=RD_plane+s_R'*s_D;
                end;
            end
            
            % Filling Output RD plane variable
            Scena.Radar(i_radar).Time(i_Time).RD_plane=RD_plane;
            % nR(id) x nV(id), [lin power], for an fftshifted
            % Doppler grid, normalized with mean noise power
            % For synthetic data only:
            burstData(burst_counter).vid = RD_plane;
        else
            burstData(burst_counter).vid = [];
        end;
    end
end;%i_radar

RadarNet_Data_FileName='RadarNet-4-radars.mat';
save(RadarNet_Data_FileName,'burstData','sensor')

%%
visualize_model(Scena);


