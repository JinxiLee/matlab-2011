clear all

%% Note

% Variables:
% Wavelength-dep: n, k
% Wavelength-indep: d0, d, eff_selfinter, eff_refsam < 這些好像省不了
range_specified=1;

%% Setting

Load_Previous_Result=0;

n_should=1.5;
Wavelength_Center=1300;

Number_of_loop=1;

Wavelength_Considered_Min=1001;          %nm
Wavelength_Considered_Max=1999;


Max_Wavelength=2000;             %nm
Min_Wavelength=1000;             %nm
N_f=8192;
N_t=8192*8;

Number_of_variable=3;           %Wavelength indep. variables

%% Data Loading

cd('D:\TuanShu\');
%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Data=importdata('111118_Green 5-3.jws.txt');
Data_Reference=importdata('111118_Green 5-3.jws.txt');

Wavelength=Data(end:-1:1,1);           %nm
Spectrum_Old=Data(end:-1:1,2)/100;
Spectrum_Reference_Old=Data_Reference(end:-1:1,2);

C=3E8;

Frequency_Old=C./(Wavelength*1E-9);
Spectrum_Frequency=(Spectrum_Old);        %/max(Spectrum_Old.*((Wavelength*1E-9).^2)/C);  N more normalization
Spectrum_Reference_Frequency=(Spectrum_Reference_Old);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);

Max_Frequency=C/(Min_Wavelength*1E-9);             %Hz
Min_Frequency=C/(Max_Wavelength*1E-9);             %Hz

Frequency_Center=C/(Wavelength_Center*1E-9);

Frequency_Considered_Min=C/(Wavelength_Considered_Max*1E-9);             %Hz
Frequency_Considered_Max=C/(Wavelength_Considered_Min*1E-9);             %Hz

Frequency=0:Max_Frequency/(N_f-1):Max_Frequency;
Frequency=Frequency';

Frequency_Center_Index=find(Frequency>Frequency_Center,1,'first');

Frequency_Considered_Min_Index=find(Frequency>Frequency_Considered_Min,1,'first');
Frequency_Considered_Max_Index=find(Frequency>Frequency_Considered_Max,1,'first');

Spectrum=interp1(Frequency_Old,Spectrum_Frequency,Frequency);
Spectrum_Reference=interp1(Frequency_Old,Spectrum_Reference_Frequency,Frequency);

Spectrum(isnan(Spectrum))=0;
Spectrum(Frequency<Min_Frequency)=0;

Spectrum_Reference(isnan(Spectrum))=0;
Spectrum_Reference(Frequency<Min_Frequency)=0;

Spectrum((N_f+1):N_t)=0;
Spectrum_Reference((N_f+1):N_t)=0;

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=C*Time;
Position_micron=Position*1E6;

Spectrum_s=Spectrum;
Spectrum_s(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index)=Spectrum_s(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index)-mean(Spectrum(Frequency_Considered_Min_Index:Frequency_Considered_Max_Index));
Spectrum_s(1:Frequency_Considered_Min_Index-1)=0;
Spectrum_s((Frequency_Considered_Max_Index+1):end)=0;


Signal=fft(Spectrum_s);
Signal_Reference=fft(Spectrum_Reference);
%% Again FD
Spectrum=Spectrum(1:N_f);     %2 to compensate the amplitude reduction from Hilbert transform
%% Divided by Spectrum_Reference

[maxvalue maxindex]=max(abs(Signal));
Thickness_0=Position(maxindex);
Thickness_temp=(Thickness_0-(Thickness_0*0.1)):Thickness_0*00.1:(Thickness_0+(Thickness_0*0.1));

%Spectrum=Spectrum-Spectrum_Reference;
%% Theory - Sample Model (n1 - n - n2)

n1=1;

% Assumed n2 = BK7
C1 = 1.03961212; 
C2 = 0.00600069867; 
C3 = 0.231792344; 
C4 = 0.0200179144; 
C5 = 1.01046945; 
C6 = 103.560653;


Wavelength_micron=(C./Frequency)*1E6;
n_bk7=(C1*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C2)+C3*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C4)+C5*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C6)+1).^0.5;

n_bk7=abs(n_bk7);

n_bk7(isnan(n_bk7))=0;
n2=n_bk7;

%% To generate the original spectrum, assuming reference is BK7, too

r_BK7=((1-n_bk7)./(n_bk7+1));
t_AN100=(2*(n_bk7)./(n_bk7+1));

%% Fitting

n_final(1:length(Frequency),1:length(n_should))=1.5;
k_final(1:length(Frequency),1:length(n_should))=0;
n(1:length(Frequency),1)=1.5;
k(1:length(Frequency),1)=0;

Spectrum_check(1:length(Spectrum),1:length(n_should))=0;
Spectrum_Ratio_check(1:length(Spectrum),1:length(n_should))=0;
T_check(1:length(Spectrum),1:length(n_should))=0;
value_temp=99999999999999;
for p=1:length(Thickness_temp)                                                        % p: Loop 1, for different solutions
    
    % n Template Generation
    n_o=(n_should-0.2):0.0005:(n_should+0.2);
    k_o=0:0.1:0.1;
    k_o=k_o';

    n_empty=n_o;
    n_empty(:)=1;
    k_empty=k_o;
    k_empty(:)=1;
    n_temp=k_empty*n_o;
    k_temp=k_o*n_empty;
    
% initial guesses:

    index_2_old=0;
    %d_check(1:length(Spectrum))=0;
    value_total=0;
    for j=1:(Frequency_Considered_Max_Index-Frequency_Considered_Min_Index+1)
            
        if j <= ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)
            index_now=Frequency_Center_Index+j-1;
        elseif j > ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)           %direction reverse!!
            index_now=-j+((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)+Frequency_Center_Index;
        end
            
            % About the conditions
        r1_r=((n_temp+i*k_temp)-n1)./(n_temp+n1+i*k_temp);
        t1=2*(n1)./(n_temp+n1+i*k_temp);
        t1_r=2*(n_temp+i*k_temp)./(n_temp+n1+i*k_temp);
        t2=2*(n_temp+i*k_temp)./(n_temp+n2(index_now)+i*k_temp);
        r2=((n_temp+i*k_temp)-n2(index_now))./((n_temp+i*k_temp)+n2(index_now));   
        d=exp(i*2*pi.*Frequency(index_now)/C.*(n_temp+i*k_temp).*Thickness_temp(p));   %注意! -1*n!
        T_Temp=real((t_AN100(index_now).^2).*(t1.*t2.*d)./(1+r1_r.*r2.*(d.^2))).^2;
        
       Merit=(real(T_Temp-Spectrum(index_now)).^2);

        [value index_1]=min(Merit);
        [value index_2]=min(value);
        index_2_temp=index_2;

        T_check(index_now)=T_Temp(index_1(index_2_temp), index_2);
        
        
    end

    if value_total < value_temp
        Thickness_final=Thickness_temp(p);
        T_check_Final=T_check;
        value_temp=value_total;
    end
    
end
plot(Frequency,Spectrum,Frequency,T_check);
%% Checking
