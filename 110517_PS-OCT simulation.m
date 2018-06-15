clear all;
close all;
clc;

c=3E8;      %m/sec


%% spectrum

center_wavelength=0.76;     %micron
bandwidth=0.18;              %micron


center_frequency=c/center_wavelength*1E6;   %Hz
bandwidth_f=c*bandwidth/(center_wavelength^2)*1E6;    %Hz

frequency=0:(center_wavelength+2*bandwidth_f)/4999:(center_frequency+2*bandwidth_f);

spectrum_f=gaussmf(frequency,[bandwidth_f/2/(2*log(2))^0.5 center_frequency]);

plot(frequency, spectrum_f);


%% sample

axial_position0=100;        %initial position of following array, micron

axial_position=0:500:       %micron
axial_profile(1:length(axial_position))=0;
axial_profile(axial_position==1)=1;
axial_profile(axial_position==10)=1;
axial_profile(axial_position==100)=1;


%
