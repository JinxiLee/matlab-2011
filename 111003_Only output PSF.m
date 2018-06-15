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


cd('D:\111003\green 0.2Hz 1V 100000points\');     %4 (2Hz)
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


dlmwrite(sprintf('Signal_%i.txt',jj),Signal,'delimiter','\t','newline','pc');

dlmwrite(sprintf('Signal_Upper_%i.txt',jj),Signal_Upper,'delimiter','\t','newline','pc');
%dlmwrite(sprintf('Signal_Upper_f_phase_%i.txt',jj),Signal_Upper_f_phase,'delimiter','\t','newline','pc');

dlmwrite(sprintf('Signal_Lower_%i.txt',jj),Signal_Lower,'delimiter','\t','newline','pc');
%dlmwrite(sprintf('Signal_Lower_f_phase_%i.txt',jj),Signal_Lower_f_phase,'delimiter','\t','newline','pc');

end

dlmwrite('position.txt',position,'delimiter','\t','newline','pc');

plot(position,Signal,position,Signal_Upper,position,Signal_Lower);
plot(wavelength_micron(wavelength_micron<1),n_final(wavelength_micron<1),wavelength_micron(wavelength_micron<1),k_final(wavelength_micron<1));

plot(wavelength_micron(wavelength_micron<1),Signal_Upper_f_amp(wavelength_micron<1),wavelength_micron(wavelength_micron<1),Signal_Lower_f_amp(wavelength_micron<1));