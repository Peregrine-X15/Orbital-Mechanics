clear all
close all
clc

r_Earth_o=[1.525E11 9.242E10 2.2236E7];
r1=r_Earth_o;
v_Earth_o=[-1.9163E4 2.3196E4 -4.2957E-1];
r_Earth_o_mag=norm(r_Earth_o);
v_Earth_o_mag=norm(v_Earth_o);

r_Mars_f=[-2.2242E11 -9.7241E10 3.441E9];
r2=r_Mars_f;
v_Mars_f=[1.0608E4 -2.0205E4 -6.8343E2];
r_Mars_f_mag=norm(r_Mars_f);
v_Mars_f_mag=norm(v_Mars_f);

delta_theta=acos(dot(r_Earth_o,r_Mars_f)/(r_Mars_f_mag*r_Earth_o_mag));
delta_theta=2*pi-delta_theta;

M_Earth=5.97219*10^24;
M_Mars=6.39E23;
G_univ = 6.6743*10^-11;
r_EM=norm(r_Earth_o-r_Mars_f);

M_sun=1.989E30;
mu = G_univ*(M_sun);

%orbital elements from given
epsilon_Earth = (v_Earth_o_mag^2)/2 - mu/r_Earth_o_mag; % Specific orbital energy
a_Earth = -mu/(2*epsilon_Earth); % Semi-major axis
h_Earth = norm(cross(r_Earth_o, v_Earth_o)); % Specific angular momentum
e_Earth = sqrt(1 + (2*epsilon_Earth*h_Earth^2)/mu^2); % Eccentricity

epsilon_Mars = (v_Mars_f_mag^2)/2 - mu/r_Mars_f_mag; % Specific orbital energy
a_Mars = -mu/(2*epsilon_Mars); % Semi-major axis
h_Mars = norm(cross(r_Mars_f, v_Mars_f)); % Specific angular momentum
e_Mars = sqrt(1 + (2*epsilon_Mars*h_Mars^2)/mu^2); % Eccentricity


chord=norm(r2-r1);
s=(norm(r1)+norm(r2)+chord)/2;
% c=2*s-r1-r2;


one_div_p=(chord*sin(delta_theta))/((1-cos(delta_theta))*r_Earth_o_mag*r_Mars_f_mag*sin(delta_theta));
p=1/one_div_p;

e = (r_Mars_f_mag - r_Earth_o_mag) / chord;

 a=0.5*(r_Earth_o_mag+r_Mars_f_mag+chord)/2;
 h=sqrt(mu*p);

 
true_anomaly = linspace(0, 2*pi, 500);

% Parametric equation for the transfer orbit
r_MinOrb = p ./ (1 + e * cos(true_anomaly));


x_minOrb = r_MinOrb .* cos(true_anomaly);
y_minOrb = r_MinOrb .* sin(true_anomaly);
z_minOrb=zeros(size(x_minOrb));
% Plot the orbits and points
%earth and mars
% earth_orbit_x = a_Earth * cos(true_anomaly);
% earth_orbit_y = a_Earth * sin(true_anomaly);
% earth_orbit_z = zeros(size(earth_orbit_x)); % 2D plot for Earth

earth_orbit_x=r_Earth_o_mag*cos(true_anomaly);
earth_orbit_y=r_Earth_o_mag*sin(true_anomaly);
earth_orbit_z=zeros(size(earth_orbit_x)); 
% Plot Mars orbit
% mars_orbit_x = a_Mars * cos(true_anomaly) ;
% mars_orbit_y = a_Mars * sin(true_anomaly);
% mars_orbit_z = zeros(size(mars_orbit_x)); % 2D plot for Mars
mars_orbit_x =r_Mars_f_mag * cos(true_anomaly) ;
mars_orbit_y = r_Mars_f_mag * sin(true_anomaly);
mars_orbit_z =zeros(size(mars_orbit_y));
figure;
hold on;


plot3(r1(1), r1(2),r1(3), 'bo', 'MarkerSize', 10, 'DisplayName', 'Earth Position');
plot3(mars_orbit_x,mars_orbit_y,mars_orbit_z,'r-','LineWidth',2,'DisplayName','Mars Orbit')
plot3(earth_orbit_x,earth_orbit_y,earth_orbit_z,'b-','LineWidth',2,'DisplayName','Earth Orbit')
plot3(r2(1), r2(2),r2(3), 'ro', 'MarkerSize', 10, 'DisplayName', 'Mars Position');
plot3(x_minOrb, y_minOrb,z_minOrb, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Minimum Energy Orbit');

% Plot settings
xlabel('X Position (km)');
ylabel('Y Position (km)');
title('Minimum Energy Transfer Orbit from Earth to Mars');
legend('show');
grid on;
axis equal;
hold off;