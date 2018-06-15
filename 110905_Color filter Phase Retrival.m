clear all;
close all;
clc;

total_OPD=114;

%% Data Loading

Signal=importdata('D:\110808\110808_2Hz B-LCD');
Reference=importdata('D:\110808\110808_2Hz mirror A-scan');
Signal=Signal(:,10);
Reference=Reference(:,10);
position=[0:total_OPD/(length(Signal)-1):total_OPD]';  

%% Pre-shifting

[maxvalue maxindex]=max(Signal,[],1);
[maxvalue_Reference maxindex_Reference]=max(Reference,[],1);

needshift=round(length(Signal)/2)-maxindex;
needshift_Reference=round(length(Reference)/2)-maxindex_Reference;

Signal=circshift(Signal,needshift);
Reference=circshift(Reference,needshift_Reference);

%% DC substracting (to avoid gap bwtween set zero and raw-data DC)

Starting_pixel=4700;
Ending_pixel=5600;

mean_in_ROI=mean(Signal(Starting_pixel:Ending_pixel));
mean_in_ROI_Reference=mean(Reference(Starting_pixel:Ending_pixel));

Signal=Signal-mean_in_ROI;
Reference=Reference-mean_in_ROI_Reference;

Signal(1:Starting_pixel,:)=0;
Signal(Ending_pixel:end,:)=0;

Reference(1:Starting_pixel,:)=0;
Reference(Ending_pixel:end,:)=0;

%% Add points

%Signal=[zeros(45000,1); Signal; zeros(45000,1)];

%Reference=[zeros(45000,1); Reference; zeros(45000,1)];

%% Shifting again for safe

%[maxvalue maxindex]=max(Signal,[],1);
%[maxvalue_Reference maxindex_Reference]=max(Reference,[],1);

%needshift=-maxindex;
%needshift_Reference=-maxindex_Reference;

%Signal=circshift(Signal,needshift);
%Reference=circshift(Reference,needshift_Reference);

%% To frequency domain

Starting_pixel_f=400;
Ending_pixel_f=1200;

Signal_f=fft(Signal,[],1);
Reference_f=fft(Reference,[],1);

Signal_f(1:Starting_pixel_f,:)=0;
Signal_f((length(Signal_f)-Starting_pixel_f+1):length(Signal_f),:)=0;
Signal_f(Ending_pixel_f:(length(Signal_f)-Ending_pixel_f),:)=0;


Reference_f(1:Starting_pixel_f,:)=0;
Reference_f((length(Signal_f)-Starting_pixel_f+1):length(Signal_f),:)=0;
Reference_f(Ending_pixel_f:(length(Signal_f)-Ending_pixel_f),:)=0;


Signal_f_amp=abs(Signal_f);
%Signal_f_phase=angle(Signal_f);

Reference_f_amp=abs(Reference_f);
%Reference_f_phase=angle(Reference_f);

X_exp=Signal_f_amp./Reference_f_amp;

c=3E8;

dx=(total_OPD*2*(1E-6))/(length(Signal)-1);
dt=2*dx/c;
f_total=1/dt;
frequency=1:f_total/(length(Signal)-1):f_total;
wavelength=c./frequency'*1E6;

%% To find the envelope of signal in time domain

Signal_f_for_Hilbert=Signal_f;
Signal_f_for_Hilbert(Ending_pixel_f:end)=0;
Signal_Hilbert=abs(ifft(Signal_f_for_Hilbert));             %% This is target, C1

%% Thory
n(1:1000,1:1000)=0;        %df is the same as frequency, but only take a portion of it (1000*df)
k(1:1000,1:1000)=0;


R_upper=( ( ( ((n-1.46).^2)  +  k.^2       )  / ( ((n+1.46).^2)  +  k.^2       )  ).^0.5 ) * exp(atan(2*k*1.46/((n.^2)+(k.^2)-1.46^2)));

R_lower=( ( ( ((n-1).^2)  +  k.^2       )  / ( ((n+1).^2)  +  k.^2       )  ).^0.5 ) * exp(atan(-2*k*1/((n.^2)+(k.^2)-1^2)));

%dx=(total_OPD*2*(1E-6))/(length(Signal)-1);

%% Initial Guess

X=X_exp;
Y=0;

%% Target

E_t_exp=Signal_Hilbert;

%% Iteration starts

iteration=0;    % count of full iteration (to C1 and to C1)

while (iteration < 10)
% to C1

R_upper=( ( ( ((n-1.46).^2)  +  k.^2       )  / ( ((n+1.46).^2)  +  k.^2       )  ).^0.5 ) * exp(atan(2*k*1.46/((n.^2)+(k.^2)-1.46^2)));
R_lower=( ( ( ((n-1).^2)  +  k.^2       )  / ( ((n+1).^2)  +  k.^2       )  ).^0.5 ) * exp(atan(-2*k*1/((n.^2)+(k.^2)-1^2)));
E_upper=Reference_f.*R_upper;          % upper
E_lower=E_upper.*exp(-2*(n+i*k)*distance_variable).*R_lower;
E_t=ifft(E_upper+E_lower);
E_t_phase=angle(E_t);
E_t=E_t_exp.*exp(E_t_phase);
E=fft(E_t);

        
        
    
        

    
    
    
end