clear all
close all
clc

% Constants
mu = 3.986004e14; % Gravitational parameter [m^3/s^2]
Re = 6378.14e3;   % Earth's equatorial radius [m]
f = 1 / 298.25;   % Earth's flattening factor
H = 3055;         % Observer altitude above sea level [m]
phi = deg2rad(20.71); % Observer latitude [radians]
lambda_G = deg2rad(156.25); % Observer longitude [radians]

% Julian Dates and Times (HST to UTC conversion)
% HST is UTC-10. Convert local solar noon (Thanksgiving Day) to UTC.
J0 = [2453458.0, 2453458.00071187, 2453458.00142371];
T0 = (J0 - 2451545) / 36525; % Julian centuries since J2000

UT = [12, 12 + 1/60 + 1.505/3600, 12 + 2/60 + 3.009/3600] - 10; % Convert HST to UTC
theta_G0 = 100.4606184 + 36000.77004*T0 + 0.000387933*T0.^2 - 2.583e-8*T0.^3; % Sidereal time at J0
theta_G = theta_G0 + 360.98564724 * (UT / 24); % Greenwich Sidereal Time [degrees]
theta_G = mod(theta_G, 360); % Wrap to [0, 360] degrees

% Local Sidereal Time at observer's longitude
LST = deg2rad(theta_G - 156.25); % Local Sidereal Time [radians]

% Observer's position (geodetic model)
R_phi = Re / sqrt(1 - (2*f - f^2) * sin(phi)^2); % Radius of curvature
Rc = R_phi + H; % Observer's radius from Earth's center [m]
Rs = (1 - f)^2 * R_phi + H; % Polar radius
R_Pos_obs = [
    Rc * cos(phi) * cos(LST); % X
    Rc * cos(phi) * sin(LST); % Y
    Rs * sin(phi) * ones(size(LST)); % Z
];

% RSO's orbital elements and anomalies
a = 9200e3; % Semi-major axis [m]
e = 0.12;   % Eccentricity
inc = deg2rad(27); % Inclination [radians]
Omega = deg2rad(180); % RA of ascending node [radians]
omega = deg2rad(42); % Argument of periapsis [radians]
f1 = deg2rad(19);    % True anomaly at first observation [radians]

% Compute RSO's position in orbital and inertial frames at first observation
r1_mag = a * (1 - e^2) / (1 + e * cos(f1)); % Magnitude of position vector
r1_O = [r1_mag * cos(f1); r1_mag * sin(f1); 0]; % Position in perifocal frame
ROI = FRE(3, omega) * FRE(1, inc) * FRE(3, Omega); % Rotation matrix (orbital -> inertial)
r1_I = ROI * r1_O; % Position in inertial frame

% Anomaly propagation using Kepler's Equation
E1 = 2 * atan(tan(f1 / 2) * sqrt((1 - e) / (1 + e))); % Eccentric anomaly
M1 = E1 - e * sin(E1); % Mean anomaly
n = sqrt(mu / a^3);    % Mean motion
M2 = M1 + n * 61.505;  % Mean anomaly at t2
M3 = M1 + n * 123.009; % Mean anomaly at t3

% Solve Kepler's equation for E2, E3
[E2, ~, ~] = solvekepler(M2, e);
[E3, ~, ~] = solvekepler(M3, e);
f2 = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E2 / 2));
f3 = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E3 / 2));

% Compute positions at second and third observations
r2_mag = a * (1 - e^2) / (1 + e * cos(f2));
r3_mag = a * (1 - e^2) / (1 + e * cos(f3));
r2_O = [r2_mag * cos(f2); r2_mag * sin(f2); 0];
r3_O = [r3_mag * cos(f3); r3_mag * sin(f3); 0];
r2_I = ROI * r2_O;
r3_I = ROI * r3_O;

% Range vectors (rho) in inertial frame
rho1 = r1_I - R_Pos_obs(:, 1);
rho2 = r2_I - R_Pos_obs(:, 2);
rho3 = r3_I - R_Pos_obs(:, 3);

% Convert rho to topocentric equatorial frame
R = [
    -sin(lambda_G), cos(lambda_G), 0;
    -sin(phi)*cos(lambda_G), -sin(phi)*sin(lambda_G), cos(phi);
    cos(phi)*cos(lambda_G), cos(phi)*sin(lambda_G), sin(phi)
];
rho1_topo = R * rho1;
rho2_topo = R * rho2;
rho3_topo = R * rho3;

% Unit vectors for rho
rho1_unit = rho1_topo / norm(rho1_topo);
rho2_unit = rho2_topo / norm(rho2_topo);
rho3_unit = rho3_topo / norm(rho3_topo);

% Azimuth and elevation
az1 = atan2(rho1_unit(2), rho1_unit(1));
az2 = atan2(rho2_unit(2), rho2_unit(1));
az3 = atan2(rho3_unit(2), rho3_unit(1));

el1 = asin(rho1_unit(3));
el2 = asin(rho2_unit(3));
el3 = asin(rho3_unit(3));

% Debugging outputs
fprintf('Azimuths (degrees): %.2f, %.2f, %.2f\n', rad2deg(az1), rad2deg(az2), rad2deg(az3));
fprintf('Elevations (degrees): %.2f, %.2f, %.2f\n', rad2deg(el1), rad2deg(el2), rad2deg(el3));

% Gauss method for orbit determination
[r_IOD, v_IOD, ~, ~] = gauss(rho1_unit, rho2_unit, rho3_unit, ...
    R_Pos_obs(:, 1), R_Pos_obs(:, 2), R_Pos_obs(:, 3), UT(1), UT(2), UT(3));

% Orbital elements from R and V
TB_orbit = OEfromRV_GPT(r_IOD, v_IOD, mu);
% fprintf('Semi-major axis: %.2f km\n', TB_orbit(1) / 1e3);
% fprintf('Eccentricity: %.5f\n', TB_orbit(2));
% fprintf('Inclination: %.2f deg\n', rad2deg(TB_orbit(3)));
% fprintf('RA of Asc. Node: %.2f deg\n', rad2deg(TB_orbit(4)));
% fprintf('Argument of Periapsis: %.2f deg\n', rad2deg(TB_orbit(5)));
% fprintf('True Anomaly: %.2f deg\n', rad2deg(TB_orbit(6)));