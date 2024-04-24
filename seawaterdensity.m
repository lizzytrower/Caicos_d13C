function [rho_sw] = seawaterdensity(T,S)
%SEAWATERDENSITY as a function of temperature and salinity
%   Function from F. J. Millero, and A. Poisson, International 
% one-atmosphere equation of state of seawater. Deep-Sea Research, 28A (6),
% 625 – 629, 1981.
%   Inputs:
%   T : temperature in degrees Celsius
%   S : practical salinity in g/kg
%   Valid for T = -2 to 40 degrees C and S = 0 to 42 g/kg

A = 0.824496 - 4.0899*10^-3*T + 7.6438*10^-5*T.^2 - 8.2467*10^-7*T.^2 ...
    + 5.3875*10^-9*T.^4;
B = -5.72466*10^-3 + 1.0227*10^-4*T - 1.6546*10^-6*T.^2;
C = 4.8314*10^-4;
rho_w = 999.842594 + 6.793952*10^-2*T - 9.09529*10^-3*T.^2 ...
    + 1.001685*10^-4*T.^3 - 1.120083*10^-6*T.^4 + 6.536336*10^-9*T.^5;
rho_sw = rho_w + A*S + B*S.^1.5 + C*S.^2;

end