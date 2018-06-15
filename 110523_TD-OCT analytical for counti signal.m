clear all;
close all;
clc;

c=3E8;      %m/sec


%% spectrum

center_wavelength=0.56;     %micron
bandwidth=0.098;              %micron         !!! PROBLEM!
%bandwidth=0.13;


center_frequency=c/center_wavelength*1E6;   %Hz
f_small=c/(center_wavelength+(0.5*bandwidth))*1E6;
f_large=c/(center_wavelength-(0.5*bandwidth))*1E6;
bandwidth_f=f_large-f_small;    %Hz

frequency=(center_frequency-2*bandwidth_f):(4*bandwidth_f)/19999:(center_frequency+2*bandwidth_f);

%spectrum_f=gaussmf(frequency,[bandwidth_f/(2*(2*log(2))^0.5) center_frequency]);

spectrum_f=gaussmf(frequency,[c*1E6*bandwidth/2/((2*log(2))^0.5)/center_wavelength^2 center_frequency]);       %for intensity!! ��field�n�}�ڸ�

wavefactor=2*pi*frequency/c;    %1/m, a array

scanning_position=140:0.1:160;    %micron

Analytical(1:length(scanning_position))=0;

for j=1:length(scanning_position)
    Analytical(j)=sum((spectrum_f.^0.5).*2./wavefactor.*sin(wavefactor.*(155-145)/2*(1E-6)).*cos(wavefactor.*((155+145)/2-scanning_position(j))*(1E-6)));
end
    
plot(scanning_position,Analytical);