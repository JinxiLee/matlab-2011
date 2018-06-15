clear all

%% Data Loading

cd('D:\');
%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
data_T_R=importdata('111129_T mirror.txt');
data_T_S=importdata('111129_T sam.txt');
data_R_R=importdata('111129_R glass.txt');
data_R_S=importdata('111129_R sam.txt');

Wavelength=data_T_R(:,1);

R=(data_R_S(:,2)./(data_R_R(:,2)));

T=data_T_S(:,2)./data_T_R(:,2);

plot(Wavelength,T);

%R: 527
%T: 532