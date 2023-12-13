function out = wavelength(frequency)
%WAVELENGTH Convert the frequency of EM wave into Wavelength
%Usage:
% out = wavelength(frequency)
%Input:
%  frequency, Hz
%Output:
% out - wavelength, [m]
%========================================================
% v.1.0 - 15.07.2012, OK@TUD
%========================================================

  out = speed_of_light/frequency;

end

