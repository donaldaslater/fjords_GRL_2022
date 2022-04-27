%**********************************************************************
% Equation of State of Sea Water At High Pressure
%**********************************************************************
function [density_seawater]=rho(temperature,salinity,depth)

% temperature in degrees celsius, salinity in ppt, depth in metres
t= temperature;
S= salinity;
rho_0=1027;
g=9.81;
P= rho_0*g*abs(depth)*1*10^(-5); % get pressure from depth and convert to bars

%======================================================================
% Calculating Secant Bulk Modulus
%======================================================================

kw= 19652.21+ 148.4206*t- 2.327105*t.^2+ 1.360477e-2*(t.^3)-5.155288e-5*(t.^4);
Aw= 3.239908+ 1.43713e-3*t+ 1.16092e-4*t.^2- 5.77905e-7*t.^3;
Bw= 8.50935e-5- 6.12293e-6*t + 5.2787e-8*(t.^2);
k0= kw + (54.6746- 0.603459*t+ 1.09987e-2*(t.^2)- 6.1670e-5*(t.^3)).*S +(7.944e-2 + 1.6483e-2*t- 5.3009e-4*(t.^2)).*(S.^1.5);
A= Aw+ (2.2838e-3- 1.0981e-5*t- 1.6078e-6*(t.^2)).*S+ 1.91075e-4*(S.^1.5);
B= Bw+ (-9.9348e-7+ 2.0816e-8*t+ 9.1697e-10*t.^2).*S;
bulk_modulus= k0+ A*P+ B*P.^2;

%======================================================================
% One atmoSphere International Equation of State [1980]
%======================================================================

A= 8.24493e-1- 4.0899e-3*t+ 7.6438e-5*t.^2- 8.2467e-7*t.^3+5.3875e-9*t.^4;
B= -5.72466e-3 + 1.0227e-4*t- 1.6546e-6*t.^2;
C= 4.8314e-4;
rho_w = 999.842594 + 6.793952e-2*t - 9.095290e-3*t.^2 + 1.001685e-4*t.^3 - 1.120083e-6*t.^4 + 6.536336e-9*t.^5;
rho_zero= rho_w+ A.*S + B.*(S.^1.5)+ C.*(S.^2);

%======================================================================
% The High Pressure International Equation of State of
% Seawater,1980
%======================================================================
density_seawater = rho_zero./(1-(P./bulk_modulus));