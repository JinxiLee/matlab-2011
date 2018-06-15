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

frequency=(center_frequency-2*bandwidth_f):(4*bandwidth_f)/9999:(center_frequency+2*bandwidth_f);

%spectrum_f=gaussmf(frequency,[bandwidth_f/(2*(2*log(2))^0.5) center_frequency]);

spectrum_f=gaussmf(frequency,[c*1E6*bandwidth/2/((2*log(2))^0.5)/center_wavelength^2 center_frequency]);       %for intensity!! 對field要開根號

plot(frequency, spectrum_f);


%% sample

axial_position=100:1:200;       %micron
axial_profile(1:length(axial_position))=sin(0.000);

    
    %axial_profile=1*gaussmf(axial_position,[1/2/(2*log(2))^0.5 150]);
    %axial_profile=(axial_position-145)*1;
    %axial_profile(axial_position==150-center_wavelength/4)=1;  
    %axial_profile(axial_position==150+center_wavelength/4)=1;
for j=(1):(length(axial_position))
    axial_profile(length(axial_position)-j+1)=sin((j)*0.0002); 
end
    %axial_profile(axial_position>155)=0;
    %axial_profile(1:length(axial_position))=0.1;
    %axial_profile(axial_position>150.9)=0;
    %axial_profile(axial_position<149.1)=0;
    %axial_profile(axial_position==150)=1;
    %axial_profile(axial_position<145)=0;
    %axial_profile(axial_position==149)=1;
    %axial_profile(axial_position==150)=1;
    %axial_profile(axial_position==151)=1;
    %axial_profile(axial_position==151)=0.1;
%axial_profile(axial_position==1)=1;
%axial_profile(axial_position==150)=0.1;


%% reference (scanning)

scanning_position=140:0.01:160;    %micron


%% interference term 

wavefactor=2*pi*frequency/c;    %1/m, a array
Eszk(1:length(axial_position),1:length(wavefactor))=0;

for j=1:length(wavefactor)
    %Esk(j)=sum(axial_profile.*spectrum_f(j).*exp(i*2.*wavefactor(j).*(axial_position)/1E6));
     Eszk(:,j)=axial_profile.*(spectrum_f(j)^0.5).*exp(i*2.*wavefactor(j).*(axial_position)/1E6);
end
     Esk=sum(Eszk,1);

 %   Esz=sum(spectrum_f.*exp(i*2*wavefactor*(10005)/1E6));

Erzk(1:length(scanning_position),1:length(wavefactor))=0;

inter(1:length(scanning_position))=0;

for j=1:length(wavefactor)
    Erzk(:,j)=(spectrum_f(j)^0.5).*exp(i*2*wavefactor(j)*(scanning_position)/1E6);    %array    
end
       Erzk_DC=(spectrum_f.^0.5);

interz(1:length(scanning_position))=0;

for j=1:length(scanning_position)
    interz(j)=sum((Erzk(j,:)+Esk).*conj(Erzk(j,:)+Esk));
end
    interz_DC=sum((Erzk_DC+Esk).*conj(Erzk_DC+Esk));


%% Hilbert transform

interz_DC_subtracted=(interz-interz_DC)';
interz_hilberted=abs(hilbert(interz_DC_subtracted));


%% Analytic

Analytical(1:length(scanning_position))=0;

for j=1:length(scanning_position)
    Analytical(j)=sum((spectrum_f.^0.5).*2./wavefactor.*sin(wavefactor.*(155-145)/2*(1E-6)).*cos(wavefactor.*((155+145)/2-scanning_position(j))*(1E-6)));
end
    
interz=interz';
interz_hilbert=abs(hilbert(interz_DC_subtracted));
scanning_position=scanning_position';

%plot(scanning_position,log(interz_hilbert));
plot(scanning_position,interz_DC_subtracted,scanning_position,interz_hilberted);