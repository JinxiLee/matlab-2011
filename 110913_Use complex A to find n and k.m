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

A_exp=Signal_Lower_f./Signal_Upper_f;

A_exp_abs=abs(A_exp);

A_exp_phase=angle(A_exp);


A_exp(isnan(A_exp))=0;

A_exp_abs(isnan(A_exp_abs))=0;
 
A_exp_phase(isnan(A_exp_phase))=0;

c=3E8;

dx=(total_OPD*2*(1E-6))/(length(Signal)-1);
dt=2*dx/c;
f_total=1/dt;
frequency=1:f_total/(length(Signal)-1):f_total;
frequency=frequency';
wavelength=c./frequency*1E6;

%plot(wavelength(500:end),A_exp_abs(500:end));




dlmwrite('D:\frequency.txt',frequency,'delimiter','\t','newline','pc');
dlmwrite('D:\Signal_Upper_f_amp.txt',Signal_Upper_f_amp,'delimiter','\t','newline','pc');
dlmwrite('D:\Signal_Lower_f_amp.txt',Signal_Lower_f_amp,'delimiter','\t','newline','pc');


%%
%%      Theoritical Calculation and fitting
%%

thickness_temp=0.0000035:0.0000005:0.0000055;
n(1:length(frequency),1:length(thickness_temp))=1.75;        %df is the same as frequency, but only take a portion of it (1000*df)
k(1:length(frequency),1:length(thickness_temp))=0;



%% Variables

n_o=1.50:0.002:2;
k_o=0:0.0002:0.05;
k_o=k_o';

n_empty=n_o;
n_empty(:)=1;
k_empty=k_o;
k_empty(:)=1;
n_temp=k_empty*n_o;
k_temp=k_o*n_empty;

n1=1;
n2=1.5;
%n=1.75;
%k=0;


%% Defination

t1=2*n1./(n_temp+n1+i*k_temp);
t1_1=2*(n_temp+i*k_temp)./(n_temp+n1+i*k_temp);           %_1 means reverse direction
t2=2*(n_temp+i*k_temp)./(n_temp+n2+i*k_temp);
t2_1=2*n2./(n_temp+n2+i*k_temp);
r1=(n1-(n_temp+i*k_temp))./(n_temp+n1+i*k_temp);
r1_1=((n_temp+i*k_temp)-n1)./(n_temp+n1+i*k_temp);
r2=((n_temp+i*k_temp)-n2)./(n_temp+n2+i*k_temp);
r2_1=(n2-(n_temp+i*k_temp))./(n_temp+n2+i*k_temp);

    % to set in loop

%% to generate a false exp value

n_false=1.55:(1.75-1.55)/(length(frequency)-1):1.75;
n_false=n_false';
k_false=0.03:(0.02-0.03)/(length(frequency)-1):0.02;
k_false=k_false';
thickness_false=0.0000041;

r1_false=(n1-(n_false+i*k_false))./(n_false+n1+i*k_false);
r2_false=((n_false+i*k_false)-n2)./(n_false+n2+i*k_false);

d_false=exp(i*2*pi./wavelength.*(n_false+i*k_false).*thickness_false); 

A_false=(d_false.^2).*r2_false./r1_false;
A_false_abs=abs(A_false);
A_false_phase=angle(A_false);

A_exp_abs=A_false_abs;
A_exp_phase=A_false_phase;

%% Iteration starts


%E_check=fft(Reference_f_for_Hilbert_amp.* ( ( ( ((n-1.46).^2)  +  k.^2
%)  ./ ( ((n+1.46).^2)  +  n.^2       )  ).^0.5 ) .* exp(i*atan(2*k*1.46./((n.^2)+(k.^2)-1.46^2)))+Reference_f_for_Hilbert_amp.*exp(2*i*(n+i*k)*2*pi.*frequency/c*distance_variable).* ( ( ( ((n-1).^2)  +  k.^2       )  ./ ( ((k+1).^2)  +  k.^2       )  ).^0.5 ) .* exp(i*atan(-2*k*1./((n.^2)+(k.^2)-1^2))));                   % auto normalize the new E to the original one; means C1 cannot change its absolute amplitude
    % to C1
value_temp=0;
for q=1:length(thickness_temp)
    value_total=0;
    for p=1:801            %the specified frequency range  
        d=exp(i*2*pi./wavelength(Starting_pixel_f+p-1).*(n_temp+i*k_temp).*thickness_temp(q)); 
        A_temp=(d.^2).*r2./r1;
        A_temp_abs=abs(A_temp);
        A_temp_phase=angle(A_temp);
        DD=((A_temp_abs-A_exp_abs(Starting_pixel_f+p-1)).^2)+((A_temp_phase-A_exp_phase(Starting_pixel_f+p-1)).^2);
        [value index_1]=min(DD);
        [value index_2]=min(value);
        n(Starting_pixel_f+p-1,q)=n_temp(index_1(index_2), index_2);
        k(Starting_pixel_f+p-1,q)=k_temp(index_1(index_2), index_2);
        value_total=value_total+abs(value)^2;
    end
    if value_temp == 0
        value_temp=value_total;
        thickness_finalindex=q;
    elseif value_total < value_temp
        thickness_finalindex=q;                % the same value for all wavelength!!
        value_temp=value_total;
    end
end

n_final=n(:,q);
k_final=k(:,q);
thickness_final=thickness_temp(thickness_finalindex);                % the same value for all wavelength!!


%% check


t1_check=2*n1./(n_final+n1+i*k_final);
t1_1_check=2*(n_final+i*k_final)./(n_final+n1+i*k_final);           %_1 means reverse direction
t2_check=2*(n_final+i*k_final)./(n_final+n2+i*k_final);
t2_1_check=2*n2./(n_final+n2+i*k_final);
r1_check=(n1-(n_final+i*k_final))./(n_final+n1+i*k_final);
r1_1_check=((n_final+i*k_final)-n1)./(n_final+n1+i*k_final);
r2_check=((n_final+i*k_final)-n2)./(n_final+n2+i*k_final);
r2_1_check=(n2-(n_final+i*k_final))./(n_final+n2+i*k_final);

d_check=exp(i*2*pi./wavelength.*(n_final+i*k_final).*thickness_temp(q)); 


A_check=(d_check.^2).*r2_check./r1_check;

A_check_abs=abs(A_check);

A_check_phase=angle(A_check);

plot(frequency,Signal_Upper_f_amp,frequency,Signal_Lower_f_amp);
plot(frequency,n_false,frequency,n_final,frequency,k_false,frequency,k_final);
