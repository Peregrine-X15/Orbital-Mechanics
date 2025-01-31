clear all
close all
clc

dataTable_1=readtable("cdata.xlsx");
dataTable_2=load("cfuture.mat");
Time=table2array(dataTable_1(:,1));
Concentration=table2array(dataTable_1(:,2));
% Time_f=table2array(dataTable_1(:,4)); Time_f=Time_f(1:61,1);
% Concentration_f=table2array(dataTable_1(:,5)); Concentration_f=Concentration_f(1:61,1);
Time_f=dataTable_2.t_future;
Concentration_f=dataTable_2.conc_future;

y=Concentration;
H=[Time, ones(size(Time))];

H_t=H'*H;
H_ty=H'*y;
coeffs=H_t\H_ty; %gives m and b (y=mx+b) 
H_inv=inv(H_t);
coeffs=H_inv*H_ty;
m=coeffs(1);
b=coeffs(2);
LSQ_lin=m*Time+b;
LSQ_lin68=m*Time_f+b; %predicting 6-8 months
err_G1b=Concentration-LSQ_lin;
err_G1_d=Concentration_f-LSQ_lin68;
%Training the model
% sumer=sum(err_G1_d)+sum(err_G1b);
% sizer=length(err_G1b)+length(err_G1_d);
sumer=sum(err_G1b);
sizer=length(err_G1b);
err_avg=sumer/sizer;
figure(1)
plot(Time,Concentration,'k')
hold on
grid on
plot(Time,LSQ_lin,'g')
plot(Time,err_G1b,'r')
legend('Original','LSQ','error')
xlabel('Time (months)')
ylabel('Concentration')
hold off

figure(2)
plot(Time,Concentration,'k')
xlabel('Time (months)')
ylabel('Concentration')
hold on
grid on
plot(Time,LSQ_lin,'g')
plot(Time_f,LSQ_lin68,'g')
plot(Time,err_G1b,'r')
plot(Time_f,err_G1_d,'r')
plot(Time_f,Concentration_f,'k')
legend('Original','LSQ Linear fit','','error')
hold off

figure(3)
plot(Time,Concentration,'k')
hold on
grid on
xlabel('Time (months)')
ylabel('Concentration')
plot(Time_f,Concentration_f,'k')
plot(Time,LSQ_lin,'b')
plot(Time_f,LSQ_lin68,'b')
plot(Time_f,LSQ_lin68+err_avg,'g')
plot(Time,LSQ_lin+err_avg,'g')
legend('Original','','LSQ baseline','','LSQ with correction','')

