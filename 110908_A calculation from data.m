clear all;
close all;
clc;

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

[maxvalue maxindex]=max(Signal_Upper);
needshift=round(length(Signal_Upper)/2)-maxindex;

Signal_Upper=circshift(Signal_Upper,needshift);

[maxvalue maxindex]=max(Signal_Lower);
needshift=round(length(Signal_Lower)/2)-maxindex;

Signal_Lower=circshift(Signal_Lower,needshift);


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

c=3E8;

dx=(total_OPD*2*(1E-6))/(length(Signal)-1);
dt=2*dx/c;
f_total=1/dt;
frequency=1:f_total/(length(Signal)-1):f_total;
frequency=frequency';
wavelength=c./frequency'*1E6;

plot(wavelength(500:end),A_abs(500:end));

%% To find the envelope of signal in time domain




%% Thory
n(1:length(frequency),1)=1.75;        %df is the same as frequency, but only take a portion of it (1000*df)
k(1:length(frequency),1)=0.001;


%R_upper=( ( ( ((n-1.46).^2)  +  k.^2       )  / ( ((n+1.46).^2)  +  k.^2       )  ).^0.5 ) * exp(atan(2*k*1.46/((n.^2)+(k.^2)-1.46^2)));

%R_lower=( ( ( ((n-1).^2)  +  k.^2       )  / ( ((n+1).^2)  +  k.^2       )  ).^0.5 ) * exp(atan(-2*k*1/((n.^2)+(k.^2)-1^2)));

%dx=(total_OPD*2*(1E-6))/(length(Signal)-1);

%% Initial Guess

X=X_exp;
Y=0;

%% Target

E_t_exp=Signal_Hilbert/max(Signal_Hilbert);


%% Variable

n_temp=1.6:0.001:1.8;
k_temp=0:0.001:0.2;
ww=(ones(size(k_temp,1),size(k_temp,2)))';
n_temp=ww*n_temp;
ww=(ones(size(n_temp,1),1))';
k_temp=k_temp'*ww;

%% Theoritical template

distance_variable=0.0000038;   % fixed now
R_upper=( ( ( ((n_temp-1.46).^2)  +  k_temp.^2       )  ./ ( ((n_temp+1.46).^2)  +  k_temp.^2       )  ).^0.5 ) .* exp(atan(2*k_temp*1.46./((n_temp.^2)+(k_temp.^2)-1.46^2)));
R_lower=( ( ( ((n_temp-1).^2)  +  k_temp.^2       )  ./ ( ((n_temp+1).^2)  +  k_temp.^2       )  ).^0.5 ) .* exp(atan(-2*k_temp*1./((n_temp.^2)+(k_temp.^2)-1^2)));

%% Iteration starts

iteration=0;    % count of full iteration (to C1 and to C1)

while (iteration < 10)
    % to C1
    E=ifft(E_t_exp*max(abs(fft(Reference_f_for_Hilbert_amp.* ( ( ( ((n-1.46).^2)  +  k.^2       )  ./ ( ((n+1.46).^2)  +  n.^2       )  ).^0.5 ) .* exp(atan(2*k*1.46./((n.^2)+(k.^2)-1.46^2)))+Reference_f_for_Hilbert_amp.*exp(2*i*(n+i*k)*2*pi.*frequency/c*distance_variable).* ( ( ( ((n-1).^2)  +  k.^2       )  ./ ( ((k+1).^2)  +  k.^2       )  ).^0.5 ) .* exp(atan(-2*k*1./((n.^2)+(k.^2)-1^2)))))).*exp(i*angle(fft(Reference_f_for_Hilbert_amp.* ( ( ( ((n-1.46).^2)  +  k.^2       )  ./ ( ((n+1.46).^2)  +  n.^2       )  ).^0.5 ) .* exp(atan(2*k*1.46./((n.^2)+(k.^2)-1.46^2)))+Reference_f_for_Hilbert_amp.*exp(2*i*(n+i*k)*2*pi.*frequency/c*distance_variable).* ( ( ( ((n-1).^2)  +  k.^2       )  ./ ( ((k+1).^2)  +  k.^2       )  ).^0.5 ) .* exp(atan(-2*k*1./((n.^2)+(k.^2)-1^2)))))));                   % auto normalize the new E to the original one; means C1 cannot change its absolute amplitude
    % to C2
    for p=1:801            %the specified frequency range
        E_upper=Reference_f_amp(Starting_pixel_f+p-1)*R_upper;
        E_lower=Reference_f_amp(Starting_pixel_f+p-1)*exp(2*2*pi*frequency(Starting_pixel_f+p-1)/c.*(n_temp+i*k_temp).*distance_variable).*R_lower;
        E_temp=abs(E_upper+E_lower-E(p));                   % to find the smallest E_temp
        
        
        
        
        
    
        

    
    
    
end