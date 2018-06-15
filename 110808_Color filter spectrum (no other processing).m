clear all;
close all;
clc;

total_OPD=115.8;        %micron
axial_resolution=1.5;   %micron


ratiop=[0:0.1:1];   
ratios=[0:0.1:1];


%ratiop=[1.2:0.01:1.3];   
%ratios=[1.7:0.01:1.8];
 
QQ=0.9:0.01:1.1;

start_position=61;
end_position=67;


Is_o=importdata('D:\110808\110808_2Hz R-LCD');
%Is_o=importdata('D:\110808\110808_2Hz mirror A-scan');
Is_o=Is_o(:,10);
%Is_o=importdata('D:\110426_Original PSOCT data\PALC_DC_high
%Z_1V_1_s.txt');

position=[0:total_OPD/(length(Is_o)-1):total_OPD]';  

%plot(position,Is_o);


%Is_background=[Is(1):((Is(length(Is))-Is(1))/(length(Is)-1)):Is(end)]';

%Is=Is-Is_background;

%% Background substration
%Starting_pixel=1000;

Starting_pixel=2300; %GLCD
Ending_pixel=4000;

mean_in_ROI=mean(Is_o(Starting_pixel:Ending_pixel));
%Starting_pixel=2000; 
%Ending_pixel=3000;

Is_o=Is_o-mean_in_ROI;

Is_o(1:Starting_pixel,:)=0;

Is_o(Ending_pixel:end,:)=0;


Is=[zeros(45000,1); Is_o; zeros(45000,1)];


[maxvalue maxindex]=max(Is,[],1);
needshift=-maxindex;

Is=circshift(Is,needshift);


%% Zero padding and filtering

Starting_pixel_f=4000;

Ending_pixel_f=12000;


Is_o_ftemp=fft(Is,[],1);

Is_o_ftemp(1:Starting_pixel_f,:)=0;

Is_o_ftemp((length(Is_o_ftemp)-Starting_pixel_f+1):length(Is_o_ftemp),:)=0;

Is_o_ftemp(Ending_pixel_f:(length(Is_o_ftemp)-Ending_pixel_f),:)=0;

Is_o_ftemp_amp=abs(Is_o_ftemp);
Is_o_ftemp_phase=angle(Is_o_ftemp);


Y=tan(angle(Is_o_ftemp));


%for j=1:length(Is_o_ftemp_phase)
 %   if j<length(Is_o_ftemp_phase)
  %  if Is_o_ftemp_phase(j+1) > Is_o_ftemp_phase(j) + 0.5*2*pi
   %     Is_o_ftemp_phase(j+1:end)=Is_o_ftemp_phase(j+1:end)-2*pi;
    %elseif Is_o_ftemp_phase(j+1) < Is_o_ftemp_phase(j) - 0.5*2*pi       
    %    Is_o_ftemp_phase(j+1:end)=Is_o_ftemp_phase(j+1:end)+2*pi;
    %end
    %end
%end







%dlmwrite('D:\Is_o_ftemp_amp.txt',Is_o_ftemp_amp,'delimiter','\t','newline','pc');

%dlmwrite('D:\Is_o_ftemp_phase.txt',Is_o_ftemp_phase,'delimiter','\t','newline','pc');

%% lambda genaration

c=3E8;

dx=(total_OPD*2*(1E-6))/(length(Is_o)-1);    %m

dt=2*dx/c;    %this is for real OPD

f_total=1/dt;

frequency=1:f_total/(length(Is)-1):f_total;

wavelength=c./frequency'*1E6;

plot(wavelength(500:end),abs(Is_o_ftemp_amp(500:end,:)));

dlmwrite('D:\spectrumPower_new.txt',spectrumPower_ew,'delimiter','\t','newline','pc');
%plot(position,Is_o_new,position,Is_o_new_envelope);

dlmwrite('D:\wavelength.txt',wavelength,'delimiter','\t','newline','pc');

