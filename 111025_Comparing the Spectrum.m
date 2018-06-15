clear all

%% Setting

Load_Previous_Result=1;

n_should=1.6;
Wavelength_Center=530;

Number_of_loop=1;

Wavelength_Considered_Min=500;          %nm
Wavelength_Considered_Max=600;


Max_Wavelength=800;             %nm
Min_Wavelength=300;             %nm
N_f=8192;
N_t=8192*16;
ROI_ratio=1/16;                  %only consider the first ROI_ratio data in TD
%% Data Loading

cd('D:\111019\');
%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Data=importdata('111010_Green (2500micros 325mA 100ave)r1.txt');
Data_Reference=importdata('111010_Green (2500micros 325mA 100ave)rref.txt');

Wavelength=Data(:,1);           %nm
Spectrum_Old=Data(:,2);
Spectrum_Reference_Old=Data_Reference(:,2);

C=3E8;

Frequency_Old=C./(Wavelength*1E-9);
Spectrum_Frequency=(Spectrum_Old.*((Wavelength*1E-9).^2)/C);        %/max(Spectrum_Old.*((Wavelength*1E-9).^2)/C);  N more normalization
Spectrum_Reference_Frequency=(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);

Max_Frequency=C/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=C/(Max_Wavelength*1E-9);             %Hz

Frequency_Center=C/(Wavelength_Center*1E-9);

Frequency_Considered_Min=C/(Wavelength_Considered_Max*1E-9);             %Hz
Frequency_Considered_Max=C/(Wavelength_Considered_Min*1E-9);             %Hz

Frequency=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency=Frequency';


Spectrum=interp1(Frequency_Old,Spectrum_Frequency,Frequency);
Spectrum_Reference=interp1(Frequency_Old,Spectrum_Reference_Frequency,Frequency);

Spectrum(isnan(Spectrum))=0;
Spectrum(Frequency<Min_Frequency)=0;

Spectrum_Reference(isnan(Spectrum))=0;
Spectrum_Reference(Frequency<Min_Frequency)=0;

%% First Subtract the Reference

Spectrum=Spectrum-Spectrum_Reference;

%% Spatial Filtering 

Spectrum((N_f+1):N_t)=0;

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Time=Time(1:N_f);
Position=C*Time;
Position_micron=Position*1E6;

Signal=fft(Spectrum);
Signal=Signal(1:N_f); 

Signal(1:1600)=0;

[maxvalue maxindex]=max(abs(Signal),[],1);
needshift=-maxindex;

Signal=circshift(Signal,needshift);

%% Data separation (Spatial Filtering)

Signal_1=Signal;        %DC+self interference
Signal_2=Signal;        %Upper
Signal_3=Signal;        %Lower

pixel_1=1600;               % the end of 1
pixel_2=308;               % the end of 2
pixel_3=950;               % the end of 3

Signal_1((pixel_1+1):end)=0;


Signal_1((pixel_1+1):end)=0;

Signal_2((pixel_2+1):(length(Signal_1)-pixel_1))=0;

Signal_3((pixel_3+1):end)=0;
Signal_3(1:pixel_2)=0;

%% Again FD

Spectrum_1=(ifft(Signal_1));       %can take Real, since the imaginary part should be relatively small, and it came from error
Spectrum_2=(ifft(Signal_2));
Spectrum_3=(ifft(Signal_3));

%Spectrum_1=Spectrum_1(1:N_f);
%Spectrum_2=Spectrum_2(1:N_f);
%Spectrum_3=Spectrum_3(1:N_f);
%% Theory - Sample Model (n1 - n - n2)

n1=1;

% Assumed n2 = BK7
C1 = 1.03961212; 
C2 = 0.00600069867; 
C3 = 0.231792344; 
C4 = 0.0200179144; 
C5 = 1.01046945; 
C6 = 103.560653;


Frequency_total=C/abs(2*(Position(2)-Position(1)));

Frequency=[0:Frequency_total/(length(Signal)-1):Frequency_total];%/2是因為一來一回
Frequency=Frequency';

Frequency_Center_Index=find(Frequency>Frequency_Center,1,'first');

Frequency_Considered_Min_Index=find(Frequency>Frequency_Considered_Min,1,'first');
Frequency_Considered_Max_Index=find(Frequency>Frequency_Considered_Max,1,'first');

Wavelength_micron=(C./Frequency)*1E6;
n_bk7=(C1*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C2)+C3*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C4)+C5*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C6)+1).^0.5;

n_bk7=abs(n_bk7);

n_bk7(isnan(n_bk7))=0;
n2=n_bk7;

%% To generate the original spectrum, assuming reference is BK7, too

r_BK7=((1-n_bk7)./(n_bk7+1));
r_BK7_new=r_BK7;

n_bk7_New=n_bk7;
%%

NewStyle=abs(Spectrum_3./Spectrum_2);
NewStyle_phase=angle(NewStyle);
Frequency_NewStyle=Frequency;
plot(Frequency_NewStyle,NewStyle);



%% OLD


thickness_temp=2.68E-6;
n_should=1.6;

Starting_pixel_f_considered=250;
Ending_pixel_f_considered=350;

reversed=0;

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
Data=importdata('111010_Green (2500micros 325mA 100ave)r2.txt');

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
%Signal_Carrier=Signal_Carrier(1:round(length(Signal_Carrier)*ROI_ratio));
%Signal_Envelope=Signal_Envelope(1:round(length(Signal_Envelope)*ROI_ratio));

%% Shifting


[maxvalue maxindex]=max(Signal_Envelope,[],1);
needshift=-maxindex;

Signal=circshift(Signal,needshift);
Signal_Carrier=circshift(Signal_Carrier,needshift);
Signal_Envelope=circshift(Signal_Envelope,needshift);

%% Data separation


Signal_Upper=Signal_Carrier;
Signal_Lower=Signal_Carrier;
Signal_Upper_envelope=Signal_Envelope;
Signal_Lower_envelope=Signal_Envelope;

pixel_1=310;               % the end of upper (the larger peak)
pixel_2=900;               % the end of lower (the larger peak)
pixel_3=7300;               % the start of upper (the larger peak)

Signal_Upper(pixel_1:pixel_3,1)=0;

Signal_Upper_envelope(pixel_1:pixel_3,1)=0;

Signal_Lower(1:pixel_1,1)=0;
Signal_Lower(pixel_2:end)=0;

Signal_Lower_envelope(1:pixel_1,1)=0;
Signal_Lower_envelope(pixel_2:end)=0;

[maxvalue max_U]=max(Signal_Upper_envelope(1:round(length(Signal_Upper_envelope)/2)),[],1);
[maxvalue max_L]=max(Signal_Lower_envelope(1:round(length(Signal_Lower_envelope)/2)),[],1);

OPD=abs(Position(max_U)-Position(max_L));

%% Again FD
Starting_pixel_f=100;
Ending_pixel_f=500;

Signal_Upper_f=fft(Signal_Upper);
Signal_Lower_f=fft(Signal_Lower);

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

%% Theory


c=3E8;
Frequency_Final_total=c/abs(2*(Position(2)-Position(1)));

Frequency_Final=[0:Frequency_Final_total/(length(Signal)-1):Frequency_Final_total];%/2是因為一來一回
Frequency_Final=Frequency_Final';
Wavelength_Final=c./Frequency_Final;

center_wavelength=0.56; %micron

center_frequency=c/(center_wavelength*1E-6);

center_frequency_index=find(Frequency_Final-center_frequency>0,1,'first');

plot(Wavelength_Final,Signal_Upper_f_amp);


    n1=1;

    C1 = 1.03961212; 
    C2 = 0.00600069867; 
    C3 = 0.231792344; 
    C4 = 0.0200179144; 
    C5 = 1.01046945; 
    C6 = 103.560653;


    Wavelength_micron=Wavelength_Final*1E6;

    wavelength_index_400nm=find(Wavelength_micron<0.4,1,'first');

    wavelength_index_700nm=find(Wavelength_micron<0.7,1,'first');

    n_bk7=(C1*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C2)+C3*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C4)+C5*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C6)+1).^0.5;


    n2=n_bk7;
    
    n_bk7_Old=n_bk7;
%%
OldStyle=abs(A_exp);
OldStyle_phase=angle(OldStyle);
Frequency_OldStyle=Frequency_Final;
plot(Frequency_NewStyle,NewStyle,Frequency_OldStyle,OldStyle);
plot(Frequency_NewStyle,NewStyle_phase,Frequency_OldStyle,OldStyle_phase);

plot(Frequency_NewStyle,n_bk7_New,Frequency_OldStyle,n_bk7_Old);