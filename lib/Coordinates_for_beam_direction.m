function out=Coordinates_for_beam_direction(Radar)

%==================================================================
%  v.1.0 - 10.12.2012, OK@TUD
%==================================================================

x0=Radar.position(1);
y0=Radar.position(2);


x1=x0+1.0*Radar.Range.Max_Range*sin(Radar.beam.H_beam_direction);
y1=y0+1.0*Radar.Range.Max_Range*cos(Radar.beam.H_beam_direction);


out=[x0,y0;x1,y1];

% ========== old version without structures ===============
%function out=Coordinates_for_beam_direction(position,line_direction,xlimits,ylimits)

% if (abs(line_direction)<pi/2)
%     tx=(xlimits(2)-x0)/cos(pi/2-line_direction);
% else
%     tx=(xlimits(1)-x0)/cos(pi/2-line_direction);
% end
% 
% if (line_direction>0)
%     ty=(ylimits(2)-y0)/cos(line_direction);
% else
%     ty=(ylimits(1)-y0)/cos(line_direction);
% end
% 
% if (abs(tx)>abs(ty))
% %if (tx>ty)
%     t=abs(ty);
% else
%     t=abs(tx);
% end

% x1=x0+cos(pi/2-line_direction)*t;
% y1=y0+cos(line_direction)*t;
