% clear; close all; clc

hr = 3600;
a=12000e3;
e=0.25;
Omega=deg2rad(260);
omega=deg2rad(120);
in = deg2rad(30);
f_o = deg2rad(90);
M_Earth=5.97219*10^24;
%m_sat=1000;
G_univ = 6.6743*10^-11;
mu = G_univ*M_Earth;
alpha = 1/a;
n=sqrt(mu/a^3); %rads/sec
param=a*(1-e^2);
period=2*pi/n; %secs

E_o = 2*atan(tan(f_o/2)/sqrt((1+e)/(1-e)));
Edot_o = n/(1 - e*cos(E_o));
del_t_total = period;
% tlen = 100;
tlen=500;
del_t=linspace(0,period,tlen);
r_o_Mag=param/(1+e*cos(f_o));
r_o = [r_o_Mag*cos(f_o); r_o_Mag*sin(f_o); 0];
v_o = [-a*sin(E_o); a*sqrt(1-e^2)*cos(E_o); 0];
v_o = v_o*Edot_o;
v_o_Mag = sqrt(v_o(1)^2+v_o(2)^2);

ROI = [cos(omega)*cos(Omega)-sin(omega)*cos(in)*sin(Omega) cos(omega)*sin(Omega)+sin(omega)*cos(in)*cos(Omega) sin(omega)*sin(in);
    -(sin(omega)*cos(Omega)+cos(omega)*cos(in)*sin(Omega)) -sin(omega)*sin(Omega)+cos(omega)*cos(in)*cos(Omega) cos(omega)*sin(in);
    sin(in)*sin(Omega) -sin(in)*cos(Omega) cos(in)];
r_I = ROI'*r_o;
v_I = ROI'*v_o;

tol = 1e-9;
CHI = zeros(1, tlen);
F = ones(tlen, 1); G = zeros(tlen, 1);
r_t_save=zeros(1,tlen);
% brI = zeros(3, tlen); brI(:,1) = r_I;
% bvI = zeros(3, tlen); bvI(:,1) = v_I;
for i = 2:tlen
    delt_curr = del_t(i);
    CHI_o = sqrt(mu)/r_o_Mag*delt_curr;
    sigma_o = dot(r_o, v_o)/sqrt(mu);
    tol_cal = 1;
    CHI(i)=CHI_o;
    chi_curr = CHI(i);
    ctr = 1;
    while tol_cal > tol        
        [U_o ,U_1, U_2, U_3] = UnivFns_iii(alpha, chi_curr, tol);
        del_chi = (sqrt(mu)*delt_curr - (r_o_Mag*U_1 + sigma_o*U_2 + U_3))/(r_o_Mag*U_o + sigma_o*U_1 + U_2);
        chi_curr = chi_curr + del_chi;
        tol_cal =  abs(del_chi);
        ctr = ctr + 1;
    end
    CHI(i) = chi_curr;
    %%%COMPUTE Universal funcitons at converged chi
    %%%EVAL LAGRANGE COEFFICIENTS: F, G, FD, GD
    F(i) = U_o+(alpha-1/r_o_Mag)*U_2;
    G(i) = r_o_Mag/sqrt(mu)*U_1+sigma_o/sqrt(mu)*U_2;
    %%%EVALUATE position and velocity vectors at each t (each chi)
    r_t = F(i).*r_I + G(i).*v_I;
    r_t_mag = norm(r_t);
    F_D = -sqrt(mu)/(r_t_mag*r_o_Mag)*U_1;
    G_D = 1-1/r_t_mag*U_2;
    v_t = F_D.*r_I + G_D.*v_I;
    %%%PLOT orbit
    plot3(r_t(1), r_t(2), r_t(3), 'ko', 'Markersize', 4, 'linewidth', 2);
    xlabel('x (inertial) in m');
    ylabel('y (inertial) in m');
    zlabel('z (inertial) in m');
    % plot3(r_t(1), r_t(2), r_t(3), 'k-', 'linewidth', 2);
    grid on
    hold on;
end
hold off
for p=2:tlen
 r_t_save([1 2 3],p)=F(p).*r_I+G(p).*v_I;
end
%zeroith and first time step? for b

