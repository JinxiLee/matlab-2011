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




%% Pauli spin matrices

sig1=[0 1;1 0];
sig2=[0 -i;i 0];
sig3=[1 0;0 -1];


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

position=[0:total_OPD/(length(Ip)-1):total_OPD]';    %position in OPD

%Is_second_peak_position=2805;
%Ip_second_peak_position=2798;                                       %to interpolate Ip to the same position scale as Is

%ratio=(Ip_second_peak_position-Starting_pixel)/(Is_second_peak_position-Starting_pixel);

%Ip_interpoled=interp1(position/ratio,Ip,position);

%plot(position,Is,position/ratio,Ip);

%plot(1:length(Is),Is,1:length(Ip),Ip);

%% stoke
s0(1:length(position))=0;
s1(1:length(position))=0;
s2(1:length(position))=0;
s3(1:length(position))=0;

for j=1:length(position)
    A=[]
    s0=
