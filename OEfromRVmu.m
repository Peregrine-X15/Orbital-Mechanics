function TBorbit = OEfromRVmu(rbar, vbar, mu)

%this function calculates the classical two-body orbital parameters 
%using R and V vectors at one instant as input
erad = 6378.14e3;
dtr = pi/180;

rx = rbar(1); ry = rbar(2); rz = rbar(3);
vx = vbar(1); vy = vbar(2); vz = vbar(3);
R = norm(rbar);
V = norm(vbar);
a = 1/(2/R - V^2/mu);  %SEMI-MAJOR AXIS
hbar = cross(rbar, vbar);   %ANGULAR MOMENTTUM VECTOR
h = norm(hbar);
p = h^2/mu; %PARAMETER
e = sqrt(1 - p/a);  %ECCENTRICITY
E = -mu/2/a;
sig = dot(rbar, vbar)/sqrt(mu); %QUANTITY 'SIGMA'
cbar = cross(vbar, hbar) - mu*rbar/R;
ihhat = hbar/h;
iehat = cbar/(mu*e);
imhat = cross(ihhat, iehat);
i = acos(ihhat(3)); %INCLINATION
OMEGA = atan2(ihhat(1), -ihhat(2));
if OMEGA < 0
    OMEGA = OMEGA + 2*pi;
end
omega = atan2(iehat(3), imhat(3));   %ARGUMENT OF PERIAPSIS
if omega < 0
    omega = omega + 2*pi;
end
E2 = atan2((sig/sqrt(a)), (1 - R/a));  %ECCENTRIC ANOMALY
if E2 < 0
    E2 = E2 + 2*pi;
end
%Alternate way to find eccentric anomaly
fpa = asin(dot(rbar, vbar)/(R*V));
sfpa = sin(fpa);
rdot = V*sfpa;
sinf = rdot*p/(h*e);
cosf = -(1 - p/R)/e;
f = atan2(sinf, cosf);
if f < 0
    f = f + 2*pi;
end
fac = sqrt((1-e)/(1+e));
E = 2*atan(fac*tan(f/2));
if E < 0
    E = E + 2*pi;
end
%M = E - e*sin(E);
M = E2 - e*sin(E2);   %MEAN ANOMALY

%%%return
TBorbit = [a; e; i; OMEGA; omega; f];