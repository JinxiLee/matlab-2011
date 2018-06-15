clear all;
close all;
clc;

total_OPD=114;        %micron
axial_resolution=1.5;   %micron

assumed_r=0.077;

ratiop=[0:0.1:1];   
ratios=[0:0.1:1];


%ratiop=[1.2:0.01:1.3];   
%ratios=[1.7:0.01:1.8];
 
QQ=0.9:0.01:1.1;

start_position=61;
end_position=67;


Is_o=importdata('D:\110808\110808_2Hz R-LCD');

PSF=importdata('D:\110808\110808_2Hz mirror A-scan');
%Is_o=importdata('D:\110808\110808_2Hz mirror A-scan');
Is_o=Is_o(:,10);
PSF=PSF(:,10);
%Is_o=importdata('D:\110426_Original PSOCT data\PALC_DC_high
%Z_1V_1_s.txt');

position=[0:total_OPD/(length(Is_o)-1):total_OPD]';  

%plot(position,Is_o);


%Is_background=[Is(1):((Is(length(Is))-Is(1))/(length(Is)-1)):Is(end)]';

%Is=Is-Is_background;

%% Background substration
%Starting_pixel=1000;


[maxvalue maxindex]=max(Is_o,[],1);
[maxvalue_PSF maxindex_PSF]=max(PSF,[],1);
needshift=round(length(Is_o)/2)-maxindex;
needshift_PSF=round(length(PSF)/2)-maxindex_PSF;
Is_o=circshift(Is_o,needshift);
PSF=circshift(PSF,needshift_PSF);

Starting_pixel=4700;             %4700; RLCD            %4000; %BLCD
Ending_pixel=5600;               %5600;               %5100;

mean_in_ROI=mean(Is_o(Starting_pixel:Ending_pixel));
mean_in_ROI_PSF=mean(PSF(Starting_pixel:Ending_pixel));
%Starting_pixel=2000; 
%Ending_pixel=3000;

Is_o=Is_o-mean_in_ROI;
PSF=PSF-mean_in_ROI_PSF;
Is_o(1:Starting_pixel,:)=0;

Is_o(Ending_pixel:end,:)=0;

PSF(1:Starting_pixel,:)=0;

PSF(Ending_pixel:end,:)=0;

Is=[zeros(45000,1); Is_o; zeros(45000,1)];

PSF=[zeros(45000,1); PSF; zeros(45000,1)];

[maxvalue maxindex]=max(Is,[],1);
[maxvalue_PSF maxindex_PSF]=max(PSF,[],1);
needshift=-maxindex;
needshift_PSF=-maxindex_PSF;
Is=circshift(Is,needshift);
PSF=circshift(PSF,needshift_PSF);
%% Zero padding and filtering

Starting_pixel_f=4000;

Ending_pixel_f=12000;


Is_o_ftemp=fft(Is,[],1);

Is_fPSF=fft(PSF,[],1);

Is_o_ftemp(1:Starting_pixel_f,:)=0;

Is_o_ftemp((length(Is_o_ftemp)-Starting_pixel_f+1):length(Is_o_ftemp),:)=0;

Is_o_ftemp(Ending_pixel_f:(length(Is_o_ftemp)-Ending_pixel_f),:)=0;


Is_fPSF(1:Starting_pixel_f,:)=0;

Is_fPSF((length(Is_o_ftemp)-Starting_pixel_f+1):length(Is_o_ftemp),:)=0;

Is_fPSF(Ending_pixel_f:(length(Is_o_ftemp)-Ending_pixel_f),:)=0;


Is_o_ftemp_amp=abs(Is_o_ftemp);
Is_o_ftemp_phase=angle(Is_o_ftemp);


Is_fPSF_amp=abs(Is_fPSF);
Is_fPSF_phase=angle(Is_fPSF);

c=3E8;

dx=(total_OPD*2*(1E-6))/(length(Is_o)-1);    %m

dt=2*dx/c;    %this is for real OPD

f_total=1/dt;

frequency=1:f_total/(length(Is)-1):f_total;

wavelength=c./frequency'*1E6;


%X=Is_o_ftemp_amp*0.1;       
X=(Is_o_ftemp_amp./Is_fPSF_amp)/min(Is_o_ftemp_amp(7000:9000)./Is_fPSF_amp(7000:9000))*assumed_r;       %min!!!
Y=tan(angle(Is_o_ftemp));



%% n and k calculation

k_ar=[0:0.001:1];
for j=7000:9000
    n_ar=((Y(j)*(1.5^2)+1.5*2*k_ar-Y(j)*(k_ar.^2))./Y(j)).^0.5;                %the phase relation for 1 and 1.5 is also different
    F=((((n_ar-1.5).^2)+k_ar.^2)./(((n_ar+1.5).^2)+k_ar.^2)).^0.5;      %1.5 for glass
    %F(n_ar<1.5)=-999;
    [minvalue minindex]=min(abs(F-X(j)));
    n(j)=n_ar(minindex);
    k(j)=k_ar(minindex);
end

n=n';
k=k';

plot(wavelength(7000:9000),n(7000:9000),wavelength(7000:9000),k(7000:9000),wavelength(7000:9000),Is_o_ftemp_amp(7000:9000),wavelength(7000:9000),Is_fPSF_amp(7000:9000));

plot(wavelength(7000:9000),n(7000:9000),wavelength(7000:9000),k(7000:9000));

dlmwrite('D:\n.txt',n(7000:9000),'delimiter','\t','newline','pc');

dlmwrite('D:\k.txt',k(7000:9000),'delimiter','\t','newline','pc');

%plot(wavelength(7000:9000),Is_o_ftemp_amp(7000:9000),wavelength(7000:9000),Is_fPSF_amp(7000:9000));

%for j=1:length(Is_o_ftemp_phase)
%    if j<length(Is_o_ftemp_phase)
%    if Is_o_ftemp_phase(j+1) > Is_o_ftemp_phase(j) + 0.5*2*pi
%        Is_o_ftemp_phase(j+1:end)=Is_o_ftemp_phase(j+1:end)-2*pi;
%    elseif Is_o_ftemp_phase(j+1) < Is_o_ftemp_phase(j) - 0.5*2*pi       
%        Is_o_ftemp_phase(j+1:end)=Is_o_ftemp_phase(j+1:end)+2*pi;
%    end
%    end
%end



plot(wavelength(500:end),Is_o_ftemp_phase(500:end));


dlmwrite('D:\wavelength.txt',wavelength(5000:12000),'delimiter','\t','newline','pc');
dlmwrite('D:\Is_o_ftemp_phase.txt',Is_o_ftemp_phase(5000:12000),'delimiter','\t','newline','pc');

%dlmwrite('D:\Is_o_ftemp_amp.txt',Is_o_ftemp_amp,'delimiter','\t','newline','pc');

%dlmwrite('D:\Is_o_ftemp_phase.txt',Is_o_ftemp_phase,'delimiter','\t','newline','pc');

%% lambda genaration


%plot(wavelength(500:end),abs(Is_o_ftemp_amp(500:end,:)));

%dlmwrite('D:\spectrumPower_new.txt',spectrumPower_new,'delimiter','\t','newline','pc');
%plot(position,Is_o_new,position,Is_o_new_envelope);

%dlmwrite('D:\wavelength.txt',wavelength,'delimiter','\t','newline','pc');

