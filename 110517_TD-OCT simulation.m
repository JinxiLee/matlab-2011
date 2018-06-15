clear all;
close all;
clc;

c=3E8;      %m/sec


%% spectrum

center_wavelength=0.76;     %micron
bandwidth=0.18;              %micron


center_frequency=c/center_wavelength*1E6;   %Hz
bandwidth_f=c*bandwidth/(center_wavelength^2)*1E6;    %Hz

frequency=(center_frequency-2*bandwidth_f):(4*bandwidth_f)/9999:(center_frequency+2*bandwidth_f);

spectrum_f=gaussmf(frequency,[bandwidth_f/2/(2*log(2))^0.5 center_frequency]);

plot(frequency, spectrum_f);


%% sample

axial_position=10000:0.01:10010;       %micron
axial_profile(1:length(axial_position))=0;
%axial_profile(axial_position==1)=1;
%axial_profile(axial_position==10)=1;
axial_profile(axial_position==250)=1;


%% reference (scanning)

scanning_position=10000:0.01:10010;    %micron


%% interference term (等一下.... 不同波長間該不會是out-of-phase吧......................)

wavefactor=2*pi*frequency/c;    %1/m, a array
Esz(1:length(axial_position))=0;

%for j=1:length(axial_position)
    %Esz(j)=sum(axial_profile(j).*spectrum_f.*exp(i*2.*wavefactor*(axial_position(j))/1E6));
%end

    Esz=sum(spectrum_f.*exp(i*2*wavefactor*(10010)/1E6));

Erz(1:length(scanning_position))=0;

for j=1:length(scanning_position)
    Erz(j)=sum(spectrum_f.*exp(i*2*wavefactor*(scanning_position(j))/1E6));
end

inter(1:length(scanning_position))=0;

for j=1:length(scanning_position)
inter(j)=sum((Esz+Erz(j))*conj(Esz+Erz(j)));
end

plot(scanning_position,abs(inter));