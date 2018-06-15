clear all;
close all;
clc;


%% Data pre-loading for Reference selection


cd('D:\110915\2Hz - Green 2\');     %4 (2Hz)
n_all=importdata('n_final_1.txt');
k_all=importdata('k_final_1.txt');
wavelength=importdata('wavelength_micron.txt');

n=mean(n_all,2);

k=mean(k_all,2);

thickness=2.68E-6;

frequency=3E8./(wavelength*1E-6);


n1=1;
n2(1:length(frequency),1)=1.5;        %n2 is Bk7 here

C1 = 1.03961212; 
C2 = 0.00600069867; 
C3 = 0.231792344; 
C4 = 0.0200179144; 
C5 = 1.01046945; 
C6 = 103.560653;

c=3E8;

t1_check=2*n1./(n+n1+i*k);
t1_1_check=2*(n+i*k)./(n+n1+i*k);           %_1 means reverse direction
t2_check=2*(n+i*k)./(n+n2+i*k);
t2_1_check=2*n2./(n+n2+i*k);
d=exp(i*2*pi.*frequency/c.*(-1*n+i*k).*thickness);   %ª`·N! -1*n!

T=abs(t1_check.*t1_1_check.*t2_check.*t2_1_check.*d.^2);

plot(wavelength(500:end),T(500:end));

dlmwrite('T.txt',T,'delimiter','\t','newline','pc');

