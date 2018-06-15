clear all

%% Setting


thickness_temp=2.68E-6;
n_should=1.5;

Starting_pixel_f_considered=250;
Ending_pixel_f_considered=350;

range_specified=1;

Max_Wavelength=800;             %nm
Min_Wavelength=300;             %nm
N_f=8192;
N_t=8192*16;
ROI_ratio=1/16;                  %only consider the first ROI_ratio data in TD

DC_cutoff=1000;                 %in TD

%% Data Loading

cd('D:\');
%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Data=importdata('111127_TEST.txt');

Wavelength=Data(:,1);           %nm
Spectrum=Data(:,2);

C=3E8;

Frequency=C./(Wavelength*1E-9);
Spectrum_Frequency=(Spectrum.*((Wavelength*1E-9).^2)/C)/max(Spectrum.*((Wavelength*1E-9).^2)/C);

Max_Frequency=C/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=C/(Max_Wavelength*1E-9);             %Hz

Frequency_New=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency_New=Frequency_New';
Spectrum_New=interp1(Frequency,Spectrum_Frequency,Frequency_New);

Spectrum_New(isnan(Spectrum_New))=0;
Spectrum_New(Frequency_New<Min_Frequency)=0;
%plot(Frequency_New,Spectrum_New);

%% To time domain

Spectrum_New((N_f+1):N_t)=0;

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=C*Time;
Position_micron=Position*1E6;

Signal=fft(Spectrum_New);
Spectrum_New=Spectrum_New(1:N_f);
Signal=Signal(1:round(length(Signal)*ROI_ratio));
Signal(1:DC_cutoff)=0;
Signal_Carrier=real(Signal);
Signal_Envelope=abs(Signal);

%% ROI
Time=Time(1:round(length(Time)*ROI_ratio));
Position=Position(1:round(length(Position)*ROI_ratio));
Position_micron=Position_micron(1:round(length(Position_micron)*ROI_ratio));

%% Shift

[maxvalue maxindex]=max(Signal_Envelope,[],1);
needshift=round(length(Signal_Envelope)/2)-maxindex;

Signal=circshift(Signal,needshift);
Signal_Carrier=circshift(Signal_Carrier,needshift);
Signal_Envelope=circshift(Signal_Envelope,needshift);

plot(Position_micron,Signal_Envelope);

dlmwrite('Position_micron.txt',Position_micron,'delimiter','\t','newline','pc');

dlmwrite('Signal_Carrier.txt',Signal_Carrier,'delimiter','\t','newline','pc');

dlmwrite('Signal_Envelope.txt',Signal_Envelope,'delimiter','\t','newline','pc');