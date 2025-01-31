clear all
close all
clc

mu = 3.986004e14; % Gravitational parameter [m^3/s^2]
dataTable=load("odata.mat");

brv=dataTable.brv;
br=brv(:,[1 2 3]); brx=br(:,1); bry=br(:,2); brz=br(:,3);
bv=brv(:,[4 5 6]); bvx=bv(:,1);bvy=bv(:,2); bvz=bv(:,3);
time=dataTable.tv;
time_s=time(1:125);
% time_s=time(1:80);
%shortend for samples
brx_S=brx(1:125); bry_S=bry(1:125); brz_S=brz(1:125);
bvx_S=bvx(1:125); bvy_S=bvy(1:125);bvz_S=bvz(1:125);
% brx_S=brx(1:80); bry_S=bry(1:80); brz_S=brz(1:80);
% bvx_S=bvx(1:80); bvy_S=bvy(1:80);bvz_S=bvz(1:80);

figure(1)
plot3(bvx_S,bry_S,brz_S)
legend('Position Vector')
xlabel('Position (m)')
ylabel('Position (m)')
zlabel('Position (m)')
grid on
hold off
testvar=1;
figure(2)
plot(time_s,bvx_S)
hold on
grid on
plot(time_s,bvy_S)
plot(time_s,bvz_S)
legend('v_{x}','v_{y}','v_{z}')
xlabel('time(minutes)')
ylabel('velocity (km/hr)')
%to perfrom the LSQs, there must be 3 that are done, 1 for each component.
%since the t_0 is 0, just use t for delta_t. use eqns for epsilon and
%lambda. 



yx=brx_S(1:80);
yy=bry_S(1:80);
yz=brz_S(1:80);
H_x=[ones(length(time_s(1:80)),1), time_s(1:80)',time_s(1:80)'.^2,time_s(1:80)'.^3];
H_y=H_x;
H_z=H_x;

% Solve the least squares for each component
coeffs_x = (H_x' * H_x) \ (H_x' * yx);
coeffs_y = (H_y' * H_y) \ (H_y' * yy);
coeffs_z = (H_z' * H_z) \ (H_z' * yz);
% testvar=1;
% time_norm = (time(1:250) - mean(time(1:250))) / std(time(1:250));

time_norm = (time_s - mean(time_s)) / std(time_s);

% y_tild_x = H_x * coeffs_x;
% y_tild_y = H_y * coeffs_y;
% y_tild_z = H_z * coeffs_z;
H_x=[ones(length(time_s),1), time_s',time_s'.^2,time_s'.^3];
H_y=H_x;
H_z=H_x;
time_norm = (time_s - mean(time_s)) / std(time_s);
y_tild_x = H_x * coeffs_x;
y_tild_y = H_y * coeffs_y;
y_tild_z = H_z * coeffs_z;
%real coeffs stuff
r_o=sqrt(brx_S(1)^2+bry_S(1)^2+brz_S(1)^2); c_o=r_o;
v_o=sqrt(bvx_S(1)^2+bvy_S(1)^2+bvz_S(1)^2); c_1=v_o;
eps_o = mu / r_o^3;
lambda_o = dot([brx(1) bry(1) brz(1)], [bvx(1) bvy(1) bvz(1)]) / r_o^2;
c_2=-0.5*(eps_o*r_o); c_3=1/6*(3*eps_o*lambda_o*r_o-eps_o*v_o);

coeff_Table = table(coeffs_x, coeffs_y, coeffs_z, [c_o; c_1; c_2; c_3], ...
    'VariableNames', {'Coeff_X', 'Coeff_Y', 'Coeff_Z', 'Expected'});
% figure(2)
% plot(time_norm,y_tild_x)
% hold on
% plot(time_norm,y_tild_y)
% plot(time_norm,y_tild_z)
% plot(time(1:125),y_real)
figure(3);
subplot(3, 1, 1);
plot(time_norm, brx_S, 'b', 'LineWidth', 1.5); hold on;
plot(time_norm, y_tild_x, 'r--', 'LineWidth', 1.5);
title('X Component'); legend('True Data', 'Fitted Data'); xlabel('normalized time');ylabel('position (m)');grid on;

subplot(3, 1, 2);
plot(time_norm, bry_S, 'b', 'LineWidth', 1.5); hold on;
plot(time_norm, y_tild_y, 'r--', 'LineWidth', 1.5);
title('Y Component'); legend('True Data', 'Fitted Data'); xlabel('normalized time');ylabel('position (m)');grid on;

subplot(3, 1, 3);
plot(time_norm, brz_S, 'b', 'LineWidth', 1.5); hold on;
plot(time_norm, y_tild_z, 'r--', 'LineWidth', 1.5);
title('Z Component'); legend('True Data', 'Fitted Data');xlabel('normalized time');ylabel('position (m)');grid on;