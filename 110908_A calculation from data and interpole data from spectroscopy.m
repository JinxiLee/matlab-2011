clear all;
close all;
clc;

%%
%%      The First part: A
%%


total_OPD=114;

%% Data Loading

Signal=importdata('D:\110808\110808_2Hz B-LCD');
Signal=Signal(:,10);
position=[0:total_OPD/(length(Signal)-1):total_OPD]';  

%% Pre-shifting

[maxvalue maxindex]=max(Signal,[],1);
needshift=round(length(Signal)/2)-maxindex;

Signal=circshift(Signal,needshift);

%% DC substracting (to avoid gap bwtween set zero and raw-data DC)

Starting_pixel=4800;
Ending_pixel=5600;

mean_in_ROI=mean(Signal(Starting_pixel:Ending_pixel));

Signal=Signal-mean_in_ROI;

Signal(1:Starting_pixel,:)=0;
Signal(Ending_pixel:end,:)=0;


%% Add points

%Signal=[zeros(45000,1); Signal; zeros(45000,1)];

%Reference=[zeros(45000,1); Reference; zeros(45000,1)];

%% Shifting again

[maxvalue maxindex]=max(Signal,[],1);
needshift=round(length(Signal)/2)-maxindex;

Signal=circshift(Signal,needshift);


%% Separate the lower interface

% Upper ROI

Starting_pixel_Upper=4900;
Ending_pixel_Upper=5100;

Signal_Upper(1:length(Signal),1)=0;
Signal_Upper(Starting_pixel_Upper:Ending_pixel_Upper,1)=Signal(Starting_pixel_Upper:Ending_pixel_Upper);

% Lower ROI

Starting_pixel_Lower=5100;
Ending_pixel_Lower=5500;

Signal_Lower(1:length(Signal),1)=0;
Signal_Lower(Starting_pixel_Lower:Ending_pixel_Lower,1)=Signal(Starting_pixel_Lower:Ending_pixel_Lower);


%% Re-shift to center (not really mater)

%[maxvalue maxindex]=max(Signal_Upper);
%needshift=round(length(Signal_Upper)/2)-maxindex;

%Signal_Upper=circshift(Signal_Upper,needshift);

%[maxvalue maxindex]=max(Signal_Lower);
%needshift=round(length(Signal_Lower)/2)-maxindex;

%Signal_Lower=circshift(Signal_Lower,needshift);


%% To frequency domain

Starting_pixel_f=400;
Ending_pixel_f=1200;

Signal_Upper_f=fft(Signal_Upper,[],1);
Signal_Lower_f=fft(Signal_Lower,[],1);

Signal_Upper_f(1:Starting_pixel_f,:)=0;
Signal_Upper_f((length(Signal_Upper_f)-Starting_pixel_f+1):length(Signal_Upper_f),:)=0;
Signal_Upper_f(Ending_pixel_f:(length(Signal_Upper_f)-Ending_pixel_f),:)=0;


Signal_Lower_f(1:Starting_pixel_f,:)=0;
Signal_Lower_f((length(Signal_Lower_f)-Starting_pixel_f+1):length(Signal_Lower_f),:)=0;
Signal_Lower_f(Ending_pixel_f:(length(Signal_Lower_f)-Ending_pixel_f),:)=0;


Signal_Upper_f_amp=abs(Signal_Upper_f);
%Signal_f_phase=angle(Signal_f);

Signal_Lower_f_amp=abs(Signal_Lower_f);
%Reference_f_phase=angle(Reference_f);

A=Signal_Lower_f./Signal_Upper_f;

A_abs=abs(A);

A_phase=angle(A);


c=3E8;

dx=(total_OPD*2*(1E-6))/(length(Signal)-1);
dt=2*dx/c;
f_total=1/dt;
frequency=1:f_total/(length(Signal)-1):f_total;
frequency=frequency';
wavelength=c./frequency'*1E6;

plot(wavelength(500:end),A_abs(500:end));

%%
%%      The Second part: T
%%

%% Data Loading


%Spextroscopy_Data=importdata('D:\110621 Color filter_RGB\blue2.jws.txt');



%% To find the envelope of signal in time domain


