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

axial_position=100000:0.01:100010;       %micron
axial_profile(1:length(axial_position))=0;
%axial_profile(axial_position==1)=1;
%axial_profile(axial_position==10)=1;
axial_profile(axial_position==100005)=1;


%% reference (scanning)

scanning_position=100000:0.01:100010;    %micron


%% interference term (等一下.... 不同波長間該不會是out-of-phase吧......................:但把initial position設高一點應該也就沒問題吧)

wavefactor=2*pi*frequency/c;    %1/m, a array
Esk(1:length(wavefactor))=0;

for j=1:length(wavefactor)
    %Esk(j)=sum(axial_profile.*spectrum_f(j).*exp(i*2.*wavefactor(j).*(axial_position)/1E6));
     Esk(j)=spectrum_f(j).*exp(i*2.*wavefactor(j).*(1000)/1E6);
end

 %   Esz=sum(spectrum_f.*exp(i*2*wavefactor*(10005)/1E6));

Erz(1:length(scanning_position))=0;

inter(1:length(scanning_position))=0;

for j=1:length(scanning_position)
    Erk=sum(spectrum_f.*exp(i*2*wavefactor*(scanning_position(j))/1E6));    %array
    
inter(j)=sum((Erk+Esk).*conj(Erk+Esk));
        
end


plot(scanning_position,abs(inter));