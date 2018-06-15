clear all;
close all;
clc;

lambda=0.6:0.001:1; %micron

%% Sellmeier's equation of BK7

B1=1.03961212;
C1=6.00069867*10^(-3);
B2=0.231792344;
C2=2.00179144*10^(-2);
B3=1.01046945;
C3=1.03560653*10^2;

n=(1+(B1*lambda.^2)./((lambda.^2)-C1)+(B2*lambda.^2)./((lambda.^2)-C2)+(B3*lambda.^2)./((lambda.^2)-C3)).^0.5;



%% Lensmaker's equation

R1=50*1000;  %mm
R2=-50*1000;   %mm
d=2*1000;      %mm


f=((n-1).*(1/R1-1/R2+(n-1)./n/R1/R2*d)).^-1;

f_shift=f-f(1);

plot(f_shift,lambda);