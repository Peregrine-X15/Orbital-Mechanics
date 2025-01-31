clear all
close all
clc

%Given
r_Earth_o=[1.1525E11 9.242E10 2.2236E7];
v_Earth_o=[-1.9163E4 2.3196E4 -4.2957E-1];
r_Earth_o_mag=norm(r_Earth_o);
v_Earth_o_mag=norm(v_Earth_o);



r_Mars_f=[-2.2242E11 -9.7241E10 3.4410E9];
v_Mars_f=[1.0608E4 -2.0205E4 -6.8343E2];
r_Mars_f_mag=norm(r_Mars_f);
v_Mars_f_mag=norm(v_Mars_f);


% delta_theta=acos(dot(r_Earth_o,r_Mars_f)/(r_Mars_f_mag*r_Earth_o_mag));
% delta_theta=2*pi-delta_theta;

% sin(delta_theta)
% cross(r_Earth_o,r_Mars_f)
M_Earth=5.97219*10^24;
M_Mars=6.39E23;
M_sun=1.989E30;
G_univ = 6.6743*10^-11;
r_EM=norm(r_Earth_o-r_Mars_f);

dot_product = dot(r_Earth_o, r_Mars_f);
% delta_theta = acos(dot_product / (r_Mars_f_mag * r_Earth_o_mag));

% cdth = dot(br1, br2)/(r1*r2);
% r1cr2 = cross(br1, br2);
% delta_theta = acos(cdth);
% if sign(r1cr2(end)) == -1
%     delta_theta = 2*pi - delta_theta;
% end
r1 = norm(r_Earth_o);
r2 = norm(r_Mars_f);
cdth = dot(r_Earth_o, r_Mars_f)/(r1*r2);
r1cr2 = cross(r_Earth_o, r_Mars_f);
delta_theta = acos(cdth);
if sign(r1cr2(end)) == -1
    delta_theta = 2*pi - delta_theta;
end


mu = G_univ*(M_sun);
K1=mu*(1-cos(delta_theta));
K2=sqrt(mu*r_Earth_o_mag*r_Mars_f_mag/K1)*sin(delta_theta);
% K2 = sqrt(mus*r1*r2/K1)*sin(th);
t=273;
t_secs=t*24*3600; %in seconds

%orbital elements from given
epsilon_Earth = (v_Earth_o_mag^2)/2 - mu/r_Earth_o_mag; % Specific orbital energy
a_Earth = -mu/(2*epsilon_Earth); % Semi-major axis
h_Earth = norm(cross(r_Earth_o, v_Earth_o)); % Specific angular momentum
e_Earth = sqrt(1 + (2*epsilon_Earth*h_Earth^2)/mu^2); % Eccentricity

epsilon_Mars = (v_Mars_f_mag^2)/2 - mu/r_Mars_f_mag; % Specific orbital energy
a_Mars = -mu/(2*epsilon_Mars); % Semi-major axis
h_Mars = norm(cross(r_Mars_f, v_Mars_f)); % Specific angular momentum
e_Mars = sqrt(1 + (2*epsilon_Mars*h_Mars^2)/mu^2); % Eccentricity

sigma_0=dot(r_Earth_o, v_Earth_o)/sqrt(mu);
tol=1E-9;
tlen=500; %number of iterations for z
% del_t=linspace(0,t_secs,tlen);
b=sqrt(mu)*t_secs;
tol_cal=1;
r1=r_Earth_o;
r2=r_Mars_f;
% z=zeros(1,tlen);
z=delta_theta;
z_curr(1)=z(1);
for i = 2:tlen+1
    [C,S]= Stumpf_fns_2_0(z_curr(i-1),tol);
    L=(1-z_curr(i-1).*S)/sqrt(C);
    Y=(norm(r1)+norm(r2))-K2*L;
    H=(Y/C)^(3/2)*S+K2*sqrt(Y);
    % z_curr = z(i-1);  % Use previous z value
% tol_cal=1;

% while tol_cal>tol
   % [C,S]= Stumpf_fns_2_0(z_curr,tol);
  

    %after this is for the lambert's
    % L=(1-z_curr*S)/sqrt(C);
    % Y=(norm(r1)+norm(r2))-K2*L;
    %     if Y < 0
    %         warning('Y is negative. Check input parameters.');
    %         break;  % Avoids complex values if Y is negative
    %     end
        
    % H=(Y/C)^(3/2)*S+K2*sqrt(Y);

    % kk=(C-3/2*S/C);
    % rr=(3*S/C*sqrt(Y)+K2*sqrt(C/Y));
  
    if abs(z_curr)<eps
        C=1/2; S=1/6;
        kk=(C-3/2*S/C);
        rr=(3*S/C*sqrt(Y)+K2*sqrt(C/Y));
        Y=r_Mars_f_mag+r_Earth_o_mag-K2/sqrt(C);

        dH_dz= sqrt(2)*Y^1.5/40 + K2*(sqrt(Y) + K2*sqrt(1/(2*Y)))/8;
    else
        kk=(C-3/2*S/C);
        rr=(3*S/C*sqrt(Y)+K2*sqrt(C/Y));
        dH_dz=(Y/C)^(3/2)*(1/(2*z_curr(i-1))*kk+3/4*S^2/C)+K2/8*rr;
    end
    % dH_dz=(Y/C)^(3/2)*(1/(2*z_curr)*kk+3/4*S^2/C)+K2/8*rr;

    delta_z(i)=(b-H)/dH_dz;
    z_curr(i)=z_curr(i-1)+delta_z(i);
    % tol_cal=abs(delta_z);
    
% end     
if abs(delta_z(i))<tol
    % z_curr=z(i);
    z_good=z_curr(i);
    % z_good=z(i);
end

end
%filtering the z
test=222;
% z=z(find(z<=tol));
% for j=1:length(z)
%     if abs(delta_z)>tol
%         z(j)=NaN;
%     else
%         z(j)=z(j);
%     end
% end

[C,S]= Stumpf_fns_2_0(z_good,tol);
% Semi-major axis for transfer orbit
% a_transfer = Y / (z_good .* C); % Calculate semi-major axis for transfer orbit
r1=r_Earth_o_mag;
r2=r_Mars_f_mag;
yz = r1 + r2 + K2*(z_good*S - 1)/sqrt(C);
chi = sqrt(yz/C);
alpha = z_good/chi^2;
a = 1/alpha;
U1 = chi*(1 - z_good*S);
U2 = chi^2*C;
U3 = chi^3*S;
h = sqrt(K1*r1*r2/yz);

smu = sqrt(mu);
FU = 1 - U2/r1;
GU = t_secs - U3/smu;
FdU = -smu*U1/(r1*r2);
GdU = 1 - U2/r2;

%F-G from dth SAVE FOR LATER!
% Tth = 1 - cos(delta_theta);
% sth = sin(delta_theta);
% Ft = 1 - mus*r2*Tth/h^2;
% Gt = r1*r2*sth/h;
% Fdt = mus*Tth*(mus*Tth/h^2 - 1/r1 - 1/r2)/(h*sth);
% Gdt = 1 - mus*r1*Tth/h^2;

bv1 = (r_Mars_f - FU*r_Earth_o)/GU;
bv2 = (GdU*r_Mars_f - r_Earth_o)/GU;

bh = cross(r_Mars_f, bv2);
uh = bh/h;
p = h^2/mu;
e = sqrt(1 - p/a);
i = acos(uh(3));
Om = atan2(uh(1), -uh(2));
bce = cross(bv2, bh) - mu*r_Mars_f/r2;
uc = bce/norm(bce);
uy = cross(uh, uc);
om = atan2(uc(3), uy(3));

%transfer orbit info
torbit1 = OEfromRVmu(r_Earth_o, bv1, mu);
torbit2 = OEfromRVmu(r_Mars_f, bv2, mu);
%true anomalies
f1 = torbit1(end);
f2 = torbit2(end);
if f1 > f2
    f1 = f1 - 2*pi;
end
Eorbit=OEfromRVmu(r_Earth_o, v_Earth_o, mu);
Morbit=OEfromRVmu(r_Mars_f, v_Mars_f, mu);
%for plotting where AU is put 1 to keep in metric
figure(1)
plot3(0,0,0, 'yo', 'Markersize', 18, 'Markerfacecolor', 'y');
hold on
grid on
xlabel('x-inertial (m)');
ylabel('y-inertial (m)');
zlabel('z-inertial (m)');
set(gca, 'fontweight', 'bold', 'fontsize', 10)
draworbit(Eorbit(1:5), 0, 2*pi, 1, 2, 'b', 'b', 2, 1, 1, 0);
draworbit(Morbit(1:5), 0, 2*pi, 1, 2, 'r', 'r', 2, 1, 1, 0);

orbit_te = [a, e, i, Om, om];
draworbit(orbit_te, 0, 2*pi, 1, 2, 'c', 'g', 2, 1, 1, 0);
draworbit(orbit_te, f1, f2, 1, 2, '#7E2F8E', 'g', 2, 1,1, 0);
axis([-3 3 -3 3 -3 3]*1e11)
legend('Sun', 'Earth Orbit', 'Mars Orbit', 'Transfer Ellipse', 'Transfer Arc');
% axis equal
% figure(2)