clear all;
close all;
clc;

c=3E8;      %m/sec


%% spectrum

center_wavelength=0.56;     %micron
%bandwidth=0.098;              %micron         !!! PROBLEM!
bandwidth=0.13;


center_frequency=c/center_wavelength*1E6;   %Hz
f_small=c/(center_wavelength+(0.5*bandwidth))*1E6;
f_large=c/(center_wavelength-(0.5*bandwidth))*1E6;
bandwidth_f=f_large-f_small;    %Hz

frequency=(center_frequency-2*bandwidth_f):(4*bandwidth_f)/9999:(center_frequency+2*bandwidth_f);

%spectrum_f=gaussmf(frequency,[bandwidth_f/(2*(2*log(2))^0.5) center_frequency]);

spectrum_f=gaussmf(frequency,[c*1E6*bandwidth/2/((2*log(2))^0.5)/center_wavelength^2 center_frequency]);

plot(frequency, spectrum_f);


%% sample

axial_position=140:0.1:160;       %micron
axial_profile(1:length(axial_position))=0;

    
    %axial_profile=0.1*gaussmf(axial_position,[10/2/(2*log(2))^0.5 150]);
    
    axial_profile(1:length(axial_position))=0.1;
    axial_profile(axial_position==148)=1;
    axial_profile(axial_position==152)=1;
    %axial_profile(axial_position==151)=0.1;
axial_profile(axial_position>152)=0;
axial_profile(axial_position<148)=0;
%axial_profile(axial_position==1)=1;
%axial_profile(axial_position==150)=0.1;


%% reference (scanning)

scanning_position=140:0.01:160;    %micron


%% interference term 

wavefactor=2*pi*frequency/c;    %1/m, a array
Eszk(1:length(axial_position),1:length(wavefactor))=0;

for j=1:length(wavefactor)
    %Esk(j)=sum(axial_profile.*spectrum_f(j).*exp(i*2.*wavefactor(j).*(axial_position)/1E6));
     Eszk(:,j)=axial_profile.*spectrum_f(j).*exp(i*2.*wavefactor(j).*(axial_position)/1E6);
end
     Esk=sum(Eszk,1);

 %   Esz=sum(spectrum_f.*exp(i*2*wavefactor*(10005)/1E6));

Erzk(1:length(scanning_position),1:length(wavefactor))=0;

inter(1:length(scanning_position))=0;

for j=1:length(wavefactor)
    Erzk(:,j)=spectrum_f(j).*exp(i*2*wavefactor(j)*(scanning_position)/1E6);    %array
    
        
end

interz(1:length(scanning_position))=0;


for j=1:length(scanning_position)
    interz(j)=sum((Erzk(j,:)+Esk).*conj(Erzk(j,:)+Esk));
end
    
interz_DC_subtracted=(interz-interz(1))';
interz_hilbert=abs(hilbert(interz_DC_subtracted));
scanning_position=scanning_position';

%plot(scanning_position,log(interz_hilbert));
plot(scanning_position,interz);