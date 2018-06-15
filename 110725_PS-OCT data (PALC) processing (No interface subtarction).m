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


Is_o=importdata('D:\PALC\average100\AC1V30Hz_20Hz_116um_100_1\Cut Data\s_Lee.txt');
%Is_o=importdata('D:\110426_Original PSOCT data\PALC_DC_high Z_1V_1_s.txt');

Is_o(1:500)=0;
Is_o(4500:end)=0;

Is=sum(Is_o,2)/size(Is_o,2);

%Is_background=[Is(1):((Is(length(Is))-Is(1))/(length(Is)-1)):Is(end)]';

%Is=Is-Is_background;

Ip_o=importdata('D:\PALC\average100\AC1V30Hz_20Hz_116um_100_1\Cut Data\p_Lee.txt');
%Ip_o=importdata('D:\110426_Original PSOCT data\PALC_DC_high Z_1V_1_p.txt');


Ip_o(1:500)=0;
Ip_o(4500:end)=0;

Ip=sum(Ip_o,2)/size(Ip_o,2);


position=[0:total_OPD/(length(Ip)-1):total_OPD]';  
%% background substraction
Starting_pixel=2500;

Ending_pixel=3300;

saparation_index=2900;

Is_background_slope=(Is(Ending_pixel)-Is(Starting_pixel))/(Ending_pixel-Starting_pixel-1);

Is_background=[Is(Starting_pixel)-Is_background_slope*Starting_pixel:Is_background_slope:Is(Ending_pixel)+Is_background_slope*(length(Is)-Ending_pixel)]';

%Is=Is-Is_background;

%Is(1:Starting_pixel)=0;

%Is(Ending_pixel:end)=0;

Ip_background_slope=(Ip(Ending_pixel)-Ip(Starting_pixel))/(Ending_pixel-Starting_pixel-1);

Ip_background=[Ip(Starting_pixel)-Ip_background_slope*Starting_pixel:Ip_background_slope:Ip(Ending_pixel)+Ip_background_slope*(length(Ip)-Ending_pixel)]';

%Ip=Ip-Ip_background;

%Ip(1:Starting_pixel)=0;

%Ip(Ending_pixel:end)=0;

%Is(1:Starting_pixel)=0;

%Is(Ending_pixel:end)=0;


%Ip(1:Starting_pixel)=0;

%Ip(Ending_pixel:end)=0;

%% to find first peak (for Is and Ip)

%space: position (micron)
%power: Is_carrier

[Is_peakvalue Is_peakindex]=max(Is(1:saparation_index));
Is_peakposition=position(Is_peakindex);
FWHM_Ispeak_right=position(find(Is>0.5*Is_peakvalue, 1, 'last'));   %只是大概抓, 因為這階段還沒取到envelope
FWHM_Ispeak_left=position(find(Is>0.5*Is_peakvalue, 1, 'first'));
FWHM_Ispeak=FWHM_Ispeak_right-FWHM_Ispeak_left;

position_Is_new=position-Is_peakposition;      %for interpolation of known PSF substraction


[Ip_peakvalue Ip_peakindex]=max(real(Ip(1:Is_peakindex+100)));


Ip=circshift(Ip,Is_peakindex-Ip_peakindex);    %s-pol shift by 200 pixel to minimize total phase retardation 275 for mid min 12 for minmin, 5 for poeriod min

%% first peak generation

%% to use good ratio

%% filtering and manual hilbert

Is_f_o=(fft(Is));
Is_f=Is_f_o;

Ip_f_o=(fft(Ip));
Ip_f=Ip_f_o;

start_index_of_spectrum=200;

end_index_of_spectrum=550;

Is_f(1:start_index_of_spectrum)=0;

Ip_f(1:start_index_of_spectrum)=0;


Is_f(round(length(Is_f)/2):end)=0;

Is_f(end_index_of_spectrum:round(length(Is_f)/2))=0;


Ip_f(round(length(Ip_f)/2):end)=0;

Ip_f(end_index_of_spectrum:round(length(Ip_f)/2))=0;


Is_hilberted=abs(ifft(Is_f));

Is_hilberted_real=real(ifft(Is_f));

Is_carrier=ifft(Is_f);

angs=angle(ifft(Is_f));

Ip_hilberted=abs(ifft(Ip_f));

Ip_hilberted_real=real(ifft(Ip_f));

Ip_carrier=ifft(Ip_f);

angp=angle(ifft(Ip_f));

position_new=[0:total_OPD/(length(Is_hilberted)-1):total_OPD]';  

%% phase (and to filter

phase=atan(Ip_hilberted./Is_hilberted)*180/pi;

Is_3=Is*1E3;

Ip_3=Ip*1E3;

%% Phase handling

for j=1:length(angs)
    if j<length(angs)
    if angs(j+1) > angs(j) + 0.5*2*pi
        angs(j+1:end)=angs(j+1:end)-2*pi;
    elseif angs(j+1) < angs(j) - 0.5*2*pi       
        angs(j+1:end)=angs(j+1:end)+2*pi;
    end
    
    if angp(j+1) > angp(j) + 0.5*2*pi
        angp(j+1:end)=angp(j+1:end)-2*pi;
    elseif angp(j+1) < angp(j) - 0.5*2*pi       
        angp(j+1:end)=angp(j+1:end)+2*pi;
    end
    end
end

ang_diff=angs-angp;
%ang_diff(length(ang_diff):length(ang_diff))=0;
ang_diff(length(ang_diff)-30:length(ang_diff))=0;

%plot(position_new,angs,position_new,angp,position_new,ang_diff,position_new,Ip_new*10000000); figure(gcf)

plot(position_new,Is_hilberted,position_new,Ip_hilberted,position_new,phase);

plot(position_new,phase);

dlmwrite('D:\position.txt',position,'delimiter','\t','newline','pc');
dlmwrite('D:\position_new.txt',position_new,'delimiter','\t','newline','pc');
dlmwrite('D:\Is.txt',Is,'delimiter','\t','newline','pc');
dlmwrite('D:\Ip.txt',Ip,'delimiter','\t','newline','pc');
dlmwrite('D:\Is_hilberted.txt',Is_hilberted,'delimiter','\t','newline','pc');
dlmwrite('D:\Ip_hilberted.txt',Ip_hilberted,'delimiter','\t','newline','pc');
dlmwrite('D:\Is_hilberted_real.txt',Is_hilberted_real,'delimiter','\t','newline','pc');
dlmwrite('D:\Ip_hilberted_real.txt',Ip_hilberted_real,'delimiter','\t','newline','pc');
%dlmwrite('D:\Angs.txt',angs,'delimiter','\t','newline','pc');
%dlmwrite('D:\Angp.txt',angp,'delimiter','\t','newline','pc');
%dlmwrite('D:\Ang_diff.txt',ang_diff,'delimiter','\t','newline','pc');

%plot(position,real(Is_carrier),position,real(Is_carrier_sub),position,CS_real); figure(gcf)


%plot(position,Is_carrier_sub,position,Ip_carrier_sub); figure(gcf)
