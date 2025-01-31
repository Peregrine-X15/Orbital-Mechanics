clear all
close all
clc

load brvCloud24; 

hr=3600;
a=12000e3;
M_Earth=5.97219*10^24;
G_univ = 6.6743*10^-11;
mu = G_univ*M_Earth;
n=sqrt(mu/a^3); %rads/sec
period=2*pi/n; %secs
alpha=1/a;

r_o=brvC(:,[1 2 3]);
v_o=brvC(:,[4 5 6]);
r_I=r_o;
v_I=v_o;

r_o_Mag=sqrt(r_o(:,1).^2+r_o(:,2).^2+r_o(:,3).^2);
% del_t_total = period*(1/3600); %in hours
% tlen = 500;
% del_t=[0 0.1 0.25 0.5 1 2 5 10 50 100]*3600;
del_t=[0 0.1 0.25 0.5 1 2 5 10]*3600;
tol=1e-6;
% CHI = zeros(1, tlen);
for i = 1:length(del_t)
    delt_curr = del_t(i);
    for j=2:length(r_o_Mag)
    CHI_o = sqrt(mu)/r_o_Mag(j,:)*delt_curr;
    sigma_o = dot(r_o(j,:), v_o(j,:))/sqrt(mu);
    tol_cal = 1;
    CHI(i)=CHI_o;
    chi_curr = CHI(i);
    % ctr = 1;
    while tol_cal > tol        
        [U_o ,U_1, U_2, U_3] = UnivFns_iii(alpha, chi_curr, tol);
        del_chi = (sqrt(mu)*delt_curr - (r_o_Mag(j)*U_1 + sigma_o*U_2 + U_3))/(r_o_Mag(j)*U_o + sigma_o*U_1 + U_2);
        chi_curr = chi_curr + del_chi;
        tol_cal =  abs(del_chi);
        % ctr = ctr + 1;
    end
    CHI(i) = chi_curr;
    %%%COMPUTE Universal funcitons at converged chi
    %%%EVAL LAGRANGE COEFFICIENTS: F, G, FD, GD
    F(i) = U_o+(alpha-1/r_o_Mag(j))*U_2;
    G(i) = r_o_Mag(j)/sqrt(mu)*U_1+sigma_o/sqrt(mu)*U_2;
    %%%EVALUATE position and velocity vectors at each t (each chi)
    r_t_d([1 2 3],j,i) = F(i).*r_I(j,:) + G(i).*v_I(j,:);
    r_t_mag = norm(r_t_d(i));
    F_D = -sqrt(mu)/(r_t_mag*r_o_Mag(j))*U_1;
    G_D = 1-1/r_t_mag*U_2;
    v_t = F_D.*r_I(j,:) + G_D.*v_I(j,:);
    %%%PLOT orbit
    % plot3(r_t(1), r_t(2), r_t(3), 'ko', 'Markersize', 4, 'linewidth', 2);
    % plot3(r_t(1), r_t(2), r_t(3), 'k-', 'linewidth', 2);
    grid on
    hold on;
    end
end
CHI_find_iii;
figure (1)
plot3(r_t_d(1,:,1),r_t_d(2,:,1) ,r_t_d(3,:,1),'ro','MarkerSize',4,'MarkerFaceColor','r')
hold on
plot3(r_t_save(1,1),r_t_save(2,1),r_t_save(3,1),'ko', 'Markersize', 15, 'MarkerFaceColor', 'b')
title('debris at 0 hr')
xlabel('x (inertial) in m');
ylabel('y (inertial) in m');
zlabel('z (inertial) in m');
hold off
grid on

figure(2)
plot3(r_t_d(1,:,2),r_t_d(2,:,2) ,r_t_d(3,:,2),'ro','MarkerSize',4,'MarkerFaceColor','r')
hold on
plot3(r_t_save(1,14),r_t_save(2,14),r_t_save(3,14),'ko', 'Markersize', 15, 'MarkerFaceColor', 'b')
title('debris at 0.1 hr')
xlabel('x (inertial) in m');
ylabel('y (inertial) in m');
zlabel('z (inertial) in m');
hold off
grid on

figure(3)
plot3(r_t_d(1,:,3),r_t_d(2,:,3) ,r_t_d(3,:,3),'ro','MarkerSize',4,'MarkerFaceColor','r')
hold on
plot3(r_t_save(1,35),r_t_save(2,35),r_t_save(3,35),'ko', 'Markersize', 15,'MarkerFaceColor', 'b')
title('debris at 0.25 hr')
xlabel('x (inertial) in m');
ylabel('y (inertial) in m');
zlabel('z (inertial) in m');
hold off
grid on

figure(4)
plot3(r_t_d(1,:,4),r_t_d(2,:,4) ,r_t_d(3,:,4),'ro','MarkerSize',4,'MarkerFaceColor','r')
hold on
plot3(r_t_save(1,69),r_t_save(2,69),r_t_save(3,69),'ko', 'Markersize', 15, 'MarkerFaceColor', 'b')
title('debris at 0.5 hr')
xlabel('x (inertial) in m');
ylabel('y (inertial) in m');
zlabel('z (inertial) in m');
hold off
grid on

figure(5)
plot3(r_t_d(1,:,5),r_t_d(2,:,5) ,r_t_d(3,:,5),'ro','MarkerSize',4,'MarkerFaceColor','r')
hold on
plot3(r_t_save(1,137),r_t_save(2,137),r_t_save(3,137),'ko', 'Markersize', 15,'MarkerFaceColor', 'b')
title('debris at 1 hr')
xlabel('x (inertial) in m');
ylabel('y (inertial) in m');
zlabel('z (inertial) in m');
hold off
grid on

figure(6)
plot3(r_t_d(1,:,6),r_t_d(2,:,6) ,r_t_d(3,:,6),'ro','MarkerSize',4,'MarkerFaceColor','r')
hold on
plot3(r_t_save(1,276),r_t_save(2,276),r_t_save(3,276),'ko', 'Markersize', 15,'MarkerFaceColor', 'b')
title('debris at 2 hr')
xlabel('x (inertial) in m');
ylabel('y (inertial) in m');
zlabel('z (inertial) in m');
hold off
grid on

figure(7)
plot3(r_t_d(1,:,7),r_t_d(2,:,7) ,r_t_d(3,:,7),'ro','MarkerSize',4,'MarkerFaceColor','r')
hold on
plot3(r_t_save(1,187),r_t_save(2,187),r_t_save(3,187),'ko', 'Markersize', 15, 'MarkerFaceColor', 'b')
title('debris at 5 hr')
xlabel('x (inertial) in m');
ylabel('y (inertial) in m');
zlabel('z (inertial) in m');
hold off
grid on

figure(8)
plot3(r_t_d(1,:,8),r_t_d(2,:,8) ,r_t_d(3,:,8),'ro','MarkerSize',4,'MarkerFaceColor','r')
hold on
plot3(r_t_save(1,375),r_t_save(2,375),r_t_save(3,375),'ko', 'Markersize',15,'MarkerFaceColor', 'b')
title('debris at 10 hr')
xlabel('x (inertial) in m');
ylabel('y (inertial) in m');
zlabel('z (inertial) in m');
hold off
grid on