clear all;
close all;
clc;

range_specified=1;
thickness_temp=2.68E-6;
thickness_temp=2.68E-6;


n_should=1.6;

for jj=1:10

%% 在n上似乎有ambiguity (~0.1)


%%
%%      The First part: A
%%


total_OPD=43;       %for 0.2 Hz                          %2Hz = 153.4;

%% Data Loading


cd('D:\111003\Green 0.2Hz 1V 100000points\');     %4 (2Hz)
Signal=importdata('1');

Signal=Signal(:,jj);
Signal(50000:end)=Signal(50000);
position=[0:2*total_OPD/(length(Signal)-1):2*total_OPD]';  

%% Pre-shifting (以後shift還是shift到0點附近, 雖然這樣定範圍比較麻煩, 但phase比較單純, 不會跟想像中異號)
%% (結果改了之後還是一樣 不知道為什麼) (可能跟MatLabs的fft定義有關, 接下來直接把本程式的d中的phase加個負號)

[maxvalue maxindex]=max(Signal,[],1);
needshift=-maxindex;

Signal=circshift(Signal,needshift);

%% DC filtering (to avoid gap bwtween set zero and raw-data DC)


pixel_1=4650;               % the end of upper (the larger peak)
pixel_2=11000;               % the end of lower (the larger peak)
pixel_3=95000;               % the start of upper (the larger peak)


Starting_pixel_f=100;
Ending_pixel_f=600;


Starting_pixel_f_considered=260;
Ending_pixel_f_considered=360;


Signal_temp_f=fft(Signal);

Signal_temp_f(1:Starting_pixel_f)=0;
Signal_temp_f((length(Signal_temp_f)-Starting_pixel_f+1):length(Signal_temp_f),:)=0;

Signal=real(ifft(Signal_temp_f));

Signal(pixel_2:pixel_3,:)=0;

%% Separate the lower interface

% Upper ROI

Signal_Upper=Signal;

Signal_Upper(pixel_1:pixel_3,1)=0;

Upper_max_index=1;

Upper_max_value=Signal_Upper(Upper_max_index);

% Lower ROI

Signal_Lower=Signal;

Signal_Lower(1:pixel_1,1)=0;
Signal_Lower(pixel_2:end)=0;

[Lower_max_value Lower_max_index]=max(Signal_Lower);

%% calculate the OPD of 2 interfaces (as additional constraint)

OPD=abs(Upper_max_index-Lower_max_index)*total_OPD/(length(Signal)-1)*1E-6;
center_wavelength=0.56; %micron

c=3E8;

center_frequency=c/(center_wavelength*1E-6);


%% To frequency domain

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
Signal_Upper_f_phase=angle(Signal_Upper_f);


Signal_Lower_f_amp=abs(Signal_Lower_f);
%Reference_f_phase=angle(Reference_f);
Signal_Lower_f_phase=angle(Signal_Lower_f);


A_exp=Signal_Lower_f./Signal_Upper_f;

A_exp_abs=abs(A_exp);

A_exp_phase=angle(A_exp);


A_exp(isnan(A_exp))=0;

A_exp_abs(isnan(A_exp_abs))=0;
 
A_exp_phase(isnan(A_exp_phase))=0;


Signal=circshift(Signal,round(length(Signal)/2));

Signal_Upper=circshift(Signal_Upper,round(length(Signal_Upper)/2));

Signal_Lower=circshift(Signal_Lower,round(length(Signal_Lower)/2));


%dlmwrite(sprintf('Signal_%i.txt',jj),Signal,'delimiter','\t','newline','pc');
%dlmwrite(sprintf('Signal_Upper_%i.txt',jj),Signal_Upper,'delimiter','\t','newline','pc');
%dlmwrite(sprintf('Signal_Lower_%i.txt',jj),Signal_Lower,'delimiter','\t','newline','pc');

%dlmwrite(sprintf('A_exp_abs_%i.txt',jj),A_exp_abs,'delimiter','\t','newline','pc');
%dlmwrite(sprintf('A_exp_phase_%i.txt',jj),A_exp_phase,'delimiter','\t','newline','pc');

c=3E8;

dx=(total_OPD*2*(1E-6))/(length(Signal)-1);
dt=2*dx/c;
f_total=1/dt;
frequency=1:f_total/(length(Signal)-1):f_total;
frequency=frequency';
wavelength=c./frequency;

wavelength_micron=wavelength*1E6;

%plot(wavelength(500:end),A_exp_abs(500:end));

center_frequency_index=find(frequency-center_frequency>0,1,'first');

%dlmwrite('D:\frequency.txt',frequency,'delimiter','\t','newline','pc');
%dlmwrite('D:\Signal_Upper_f_amp.txt',Signal_Upper_f_amp,'delimiter','\t','newline','pc');
%dlmwrite('D:\Signal_Lower_f_amp.txt',Signal_Lower_f_amp,'delimiter','\t','newline','pc');


%%
%%      Theoritical Calculation and fitting
%%




%% Variables

n_o=(n_should-0.2):0.001:(n_should+0.2);
k_o=0:0.001:0.1;
k_o=k_o';

n_empty=n_o;
n_empty(:)=1;
k_empty=k_o;
k_empty(:)=1;
n_temp=k_empty*n_o;
k_temp=k_o*n_empty;

n1=1;
%n2=1.46;        %n2 is Bk7 here

%n2 - 1 = C1λ2/(λ2-C2) + C3λ2/(λ2-C4) + C5λ2/(λ2-C6)


C1 = 1.03961212; 
C2 = 0.00600069867; 
C3 = 0.231792344; 
C4 = 0.0200179144; 
C5 = 1.01046945; 
C6 = 103.560653;


wavelength_index_400nm=find(wavelength_micron<0.4,1,'first');

wavelength_index_700nm=find(wavelength_micron<0.7,1,'first');

n_bk7=(C1*(wavelength_micron.^2)./((wavelength_micron.^2)-C2)+C3*(wavelength_micron.^2)./((wavelength_micron.^2)-C4)+C5*(wavelength_micron.^2)./((wavelength_micron.^2)-C6)+1).^0.5;



%dlmwrite('n_bk7.txt',n_bk7,'delimiter','\t','newline','pc');
%dlmwrite(sprintf('wavelength_micron_%i.txt',jj),wavelength_micron,'delimiter','\t','newline','pc');
%dlmwrite(sprintf('frequency_%i.txt',jj),frequency,'delimiter','\t','newline','pc');

%plot(wavelength_micron(wavelength_index_1000nm:wavelength_index_400nm),n_bk7(wavelength_index_1000nm:wavelength_index_400nm));

n2=n_bk7;
%n=1.75;
%k=0;
%% Loading Spectroscopy data



%% Defination

%t1=2*n1./(n_temp+n1+i*k_temp);
%t1_1=2*(n_temp+i*k_temp)./(n_temp+n1+i*k_temp);           %_1 means reverse direction
%t2=2*(n_temp+i*k_temp)./(n_temp+n2+i*k_temp);
%t2_1=2*n2./(n_temp+n2+i*k_temp);
r1=(n1-(n_temp+i*k_temp))./(n_temp+n1+i*k_temp);
%r1_1=((n_temp+i*k_temp)-n1)./(n_temp+n1+i*k_temp);
%r2=((n_temp+i*k_temp)-n2)./(n_temp+n2+i*k_temp);
%r2_1=(n2-(n_temp+i*k_temp))./(n_temp+n2+i*k_temp);


%thickness_temp=2E-6:0.1E-6:3E-6;             %0.0000020:0.0000001:0.000003;

n_should_index=find(n_o-n_should>0,1,'first');


n(1:length(frequency),1:length(thickness_temp))=1.75;        %df is the same as frequency, but only take a portion of it (1000*df)
k(1:length(frequency),1:length(thickness_temp))=0;

    % to set in loop

%% to generate a false exp value

%n_false=1.5:(2-1.5)/(length(frequency)-1):2;
%n_false=n_false';
%k_false=0.2:(0.00-0.2)/(length(frequency)-1):0.00;
%k_false=k_false';
%thickness_false=0.0000001;

%r1_false=(n1-(n_false+i*k_false))./(n_false+n1+i*k_false);
%r2_false=((n_false+i*k_false)-n2)./(n_false+n2+i*k_false);

%d_false=exp(i*2*pi.*frequency/c.*(-1*n_false+i*k_false).*thickness_false);  %注意! -1*n!

%A_false=(d_false.^2).*r2_false./r1_false;
%A_false_abs=abs(A_false);
%A_false_phase=angle(A_false);

%A_exp_abs=A_false_abs;
%A_exp_phase=A_false_phase;

%% Iteration starts
%E_check=fft(Reference_f_for_Hilbert_amp.* ( ( ( ((n-1.46).^2)  +  k.^2
%)  ./ ( ((n+1.46).^2)  +  n.^2       )  ).^0.5 ) .* exp(i*atan(2*k*1.46./((n.^2)+(k.^2)-1.46^2)))+Reference_f_for_Hilbert_amp.*exp(2*i*(n+i*k)*2*pi.*frequency/c*distance_variable).* ( ( ( ((n-1).^2)  +  k.^2       )  ./ ( ((k+1).^2)  +  k.^2       )  ).^0.5 ) .* exp(i*atan(-2*k*1./((n.^2)+(k.^2)-1^2))));                   % auto normalize the new E to the original one; means C1 cannot change its absolute amplitude
    % to C1


%plot(wavelength(500:end),Signal_Upper_f_amp(500:end),wavelength(500:end),Signal_Lower_f_amp(500:end));
%plot(frequency,Signal_Upper_f_amp,frequency,Signal_Lower_f_amp);
%plot(frequency,n_false,frequency,n_final,frequency,k_false,frequency,k_final);


%plot(frequency,A_check_abs,frequency,A_check_phase,frequency,A_false_abs,frequency,A_false_phase);

%plot(frequency,A_check_abs,frequency,A_check_phase,frequency,A_exp_abs,frequency,A_exp_phase);
%plot(wavelength(wavelength_index_700nm:wavelength_index_400nm),Signal_Upper_f_amp(wavelength_index_700nm:wavelength_index_400nm),wavelength(wavelength_index_700nm:wavelength_index_400nm),Signal_Lower_f_amp(wavelength_index_700nm:wavelength_index_400nm));

%plot(wavelength(wavelength_index_700nm:wavelength_index_400nm),n_final(wavelength_index_700nm:wavelength_index_400nm),wavelength(wavelength_index_700nm:wavelength_index_400nm),k_final(wavelength_index_700nm:wavelength_index_400nm));


dlmwrite(sprintf('Signal_Upper_f_amp_for_wavelength_%i.txt',jj),Signal_Upper_f_amp.*(frequency.^2)/c,'delimiter','\t','newline','pc');
%dlmwrite(sprintf('Signal_Upper_f_phase_%i.txt',jj),Signal_Upper_f_phase,'delimiter','\t','newline','pc');

dlmwrite(sprintf('Signal_Lower_f_amp_for_wavelength_%i.txt',jj),Signal_Lower_f_amp.*(frequency.^2)/c,'delimiter','\t','newline','pc');
%dlmwrite(sprintf('Signal_Lower_f_phase_%i.txt',jj),Signal_Lower_f_phase,'delimiter','\t','newline','pc');

end

dlmwrite('wavelength_micron.txt',wavelength_micron,'delimiter','\t','newline','pc');

plot(position,Signal,position,Signal_Upper,position,Signal_Lower);
plot(wavelength_micron(wavelength_micron<1),n_final(wavelength_micron<1),wavelength_micron(wavelength_micron<1),k_final(wavelength_micron<1));

plot(wavelength_micron(wavelength_micron<1),Signal_Upper_f_amp(wavelength_micron<1),wavelength_micron(wavelength_micron<1),Signal_Lower_f_amp(wavelength_micron<1));