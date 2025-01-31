clear all
close all
clc

rho = 988;
rhodot= 4.86;
A = 36;
Adot= 0.59;
a = 36.6;
adot = -0.263;
theta= 40;
phi= 35;
H = 0;
mu= 398600;
[r,v]=rv_from_observe(rho,rhodot,A,Adot,a,adot,theta,phi,H);
r_mag=norm(r);
v_mag=norm(v);
% OEs=OEfromRVmu(r,v,mu);

semi_major = 1/(2/r_mag - v_mag^2/mu);  %SEMI-MAJOR AXIS
hbar = cross(r, v);   %ANGULAR MOMENTTUM VECTOR
h = norm(hbar);
p = h^2/mu; %PARAMETER
e = sqrt(1 - p/semi_major);  %ECCENTRICITY

%orbit is hyberbolic with a=7.3E4 km