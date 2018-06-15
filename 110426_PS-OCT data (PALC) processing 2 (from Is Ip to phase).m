clear all;
close all;
clc;

total_OPD=115.8;        %micron
axial_resolution=1.5;   %micron


Is_o=importdata('D:\110426_Original PSOCT data\PALC_DC_high Z_1V_1_s.txt');

Is=sum(Is_o,2)/size(Is_o,2);

%Is_background=[Is(1):((Is(length(Is))-Is(1))/(length(Is)-1)):Is(end)]';

%Is=Is-Is_background;

Ip_o=importdata('D:\110426_Original PSOCT data\PALC_DC_high Z_1V_1_p.txt');

Ip=sum(Ip_o,2)/size(Ip_o,2);

%% background substraction
Starting_pixel=2220;

Ending_pixel=2981;

Is_background_slope=(Is(Ending_pixel)-Is(Starting_pixel))/(Ending_pixel-Starting_pixel-1);

Is_background=[Is(Starting_pixel)-Is_background_slope*Starting_pixel:Is_background_slope:Is(Ending_pixel)+Is_background_slope*(length(Is)-Ending_pixel)]';

Is=Is-Is_background;

Ip_background_slope=(Ip(Ending_pixel)-Ip(Starting_pixel))/(Ending_pixel-Starting_pixel-1);

Ip_background=[Ip(Starting_pixel)-Ip_background_slope*Starting_pixel:Ip_background_slope:Ip(Ending_pixel)+Ip_background_slope*(length(Ip)-Ending_pixel)]';

Ip=Ip-Ip_background;

Is(1:Starting_pixel)=0;

Is(Ending_pixel:end)=0;


Ip(1:Starting_pixel)=0;

Ip(Ending_pixel:end)=0;

%% to normalize the position with interpolation?

plot(1:length(Ip),Ip,1:length(Is),Is);
position=[0:total_OPD/(length(Ip)-1):total_OPD]';    %position in OPD

%Is_second_peak_position=2805;
%Ip_second_peak_position=2798;                                       %to interpolate Ip to the same position scale as Is

%ratio=(Ip_second_peak_position-Starting_pixel)/(Is_second_peak_position-Starting_pixel);

%Ip_interpoled=interp1(position/ratio,Ip,position);

%plot(position,Is,position/ratio,Ip);

%plot(1:length(Is),Is,1:length(Ip),Ip);

%% to filter spectrum in FD! (and spatially saparate two peaks)

saparation_index=2670;

Is1=Is;

Is2=Is;


Ip1=Ip;

Ip2=Ip;


Is1(saparation_index:end)=0;

Ip1(saparation_index:end)=0;

Is2(1:saparation_index)=0;

Ip2(1:saparation_index)=0;


Is1_f=(fft(Is1));


Ip1_f=(fft(Ip1));


Is2_f=(fft(Is2));


Ip2_f=(fft(Ip2));

%plot(1:length(Is),log(abs(Is1_f)),1:length(Is),log(abs(Is2_f)),1:length(Is),log(abs(Ip1_f)),1:length(Is),log(abs(Ip2_f)));
%% filtering and manual hilbert

Is_f=(fft(Is));

Ip_f=(fft(Ip));

start_index_of_spectrum=200;

end_index_of_spectrum=550;

Is_f(1:start_index_of_spectrum)=0;

Ip_f(1:start_index_of_spectrum)=0;

Is_f(end_index_of_spectrum:end)=0;
Ip_f(end_index_of_spectrum:end)=0;

%Is_f(end_index_of_spectrum:length(Is_f)-end_index_of_spectrum)=0;

%Ip_f(end_index_of_spectrum:length(Ip_f)-end_index_of_spectrum)=0;



%Is_f(length(Is_f)-start_index_of_spectrum:length(Is_f))=0;

%Ip_f(length(Is_f)-start_index_of_spectrum:length(Ip_f))=0;

Is=abs(ifft(Is_f));

Ip=abs(ifft(Ip_f));

Is_f_3=real(Is_f)*1E3;

Ip_f_3=real(Ip_f)*1E3;
%% hilbert transform

%Is_h=abs(hilbert(Is));

%Ip_h=abs(hilbert(Ip));




%plot(1:length(Is),Is_f,1:length(Ip),Ip_f);

%% phase (and to filter

phase=atan(Ip./Is)*180/pi;

%plot(1:length(Is),Is,1:length(Ip),Ip);
%Ip_background=[Ip(1):(Ip(length(Is))-Ip(1))/(length(Ip)-1):Ip(end)]';

%Ip=Ip-Ip_background;

%FFTIs=fft(Is);

%FFTIp=fft(Ip);

%FFTIs(1:240)=0;

%FFTIs(4800-240:4800)=0;

%FFTIp(1:240)=0;

%FFTIp(4800-240:4800)=0;

%Is_filtered=abs(ifft(FFTIs));

%Ip_filtered=abs(ifft(FFTIp));

clear Is_o Ip_o;

%phase=(Ip-Ip(1))./(Is-Is(1));


pixel_size=total_OPD/length(Ip);   %micron

spatial_PSF=gaussmf(position,[axial_resolution 0]);

spectrum=fft(spatial_PSF);

%plot(1:length(Is),phase);

%plot(position,Is,position,Ip);

plot(position,phase);

Is_3=Is*1E3;

Ip_3=Ip*1E3;

[first_max_s first_max_index_s]=max(Is(1:saparation_index));
first_max_position_s=position(first_max_index_s);

[first_max_p first_max_index_p]=max(Ip(1:saparation_index));
first_max_position_p=position(first_max_index_p);

%plot(position,Is_filtered,position,Ip_filtered);

%plot(position,log(abs(FFTIs)),position,log(abs(FFTIp)));