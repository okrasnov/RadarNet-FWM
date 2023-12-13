function out = radar_SNR(Radar,range,RCS,Temperature)
%RADAR_SNR Calculate signal-to-noise ration
%     SNR = Pt x Gt x Gr x Lambda^2 x RCS / (4pi)^3 x R^4 x N
%       N = k x T x B x F
%Usage:
%   out = radar_SNR(Radar,range,RCS,Temperature)
%Input:
% Radar - structure, which describes the radar characteristics:
%      Radar.Pt - transmit power, [W]
%      Radar.frequency - currier frequency, [Hz]
%      Radar.B - bandwidth, [Hz}
%      Radar.lin_Gt - transmit antenna gain, [times]
%      Radar.lin_Gr - receive antenna gain, [times]
%      Radar.lin_F - noise figure of the radar receiver. [times]
%      Radar.lin_PG - radar's processing gain, [times]
% range - distance to target, [m]
% RCS - target's radar cross-section, [m]
% Temperature - ambient/system temperature, [K]
%Output:
%  out - radar's SNR in linear scale
%========================================================
% v.1.0 - 27.12.2012, OK@TUD
%========================================================


   N = Temperature.*Radar.B.*Radar.lin_F.*Boltzmann_const();
   10*log10(N)
   out=(4*pi)^3;
   out=1/out;
   out = out*Radar.Pt*Radar.lin_Gt*Radar.lin_Gr*wavelength(Radar.frequency)^2;
   out = out*RCS;
   out = out/(range)^4;
   out = out/N;
   out = out*Radar.lin_PG;
   
end

