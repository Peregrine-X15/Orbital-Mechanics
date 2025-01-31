clear all
close all
clc

J0=2460643.000;
T0=(J0-2451545)/36525;

theta_Go=100.4606184+36000.77004*T0+0.000387933*T0^2-2.583*(10^-8)*T0^3; %degrees
%double check HST to UT conversion.
% UT=10;
UT=22;
theta_G=theta_Go+360.98564724*(UT/24);

theta_Go=wrapTo360(theta_Go);
theta_G=wrapTo360(theta_G);

theta=theta_G-156.25;
% theta=deg2rad(theta);
phi=20.71;
% phi=deg2rad(phi);
Re=6378.14; %km
f=1/298.25;
R_phi=Re/(sqrt(1-(2*f-f^2)*(sind(phi))^2));
H=3055;
Rc=R_phi+H;
Rs=(1-f)^2*R_phi+H;

R_Pos=[Rc*cosd(phi)*cosd(theta);Rc*cosd(phi)*sind(theta);Rs*sind(phi)];
% R_Pos=[Rc*cos(phi)*cos(theta);Rc*cos(phi)*sin(theta);Rs*sin(phi)];
%double check hard coded numbers.
