clear all

%% Note

% Variables:
% Wavelength-dep: n, k
% Wavelength-indep: d0, d, eff_selfinter, eff_refsam < 這些好像省不了
range_specified=1;

%% Setting

Number_of_Loop=20;

Load_Previous_Result=0;

n_should=1.74;
delta_n=0.08;
Wavelength_Center=540;

Wavelength_Considered_Min=530;          %nm
Wavelength_Considered_Max=550;


Max_Wavelength=800;             %nm
Min_Wavelength=300;             %nm
N_f=8192;
N_t=8192*8;

Number_of_variable=3;           %Wavelength indep. variables

%% Data Loading

cd('D:\TEST\');
%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Data=importdata(sprintf('TEST%i.txt',350));
Data_Reference=importdata(sprintf('TEST%i.txt',50));
Data_Spectroscopy=importdata('111118_Green 5-1.jws.txt');

Wavelength=Data(:,1);           %nm
Spectrum_Old=Data(:,2);
Spectrum_Reference_Old=Data_Reference(:,2);

Spectroscopy_Old=Data_Spectroscopy(end:-1:1,2)/100;
Wavelength_Spectroscopy=Data_Spectroscopy(end:-1:1,1);           %nm

C=3E8;

Frequency_Old=C./(Wavelength*1E-9);
Frequency_Spectroscopy_Old=C./(Wavelength_Spectroscopy*1E-9);
Spectrum_Frequency=(Spectrum_Old.*((Wavelength*1E-9).^2)/C);        %/max(Spectrum_Old.*((Wavelength*1E-9).^2)/C);  N more normalization
Spectrum_Reference_Frequency=(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);    %/max(Spectrum_Reference_Old.*((Wavelength*1E-9).^2)/C);
Spectroscopy_Frequency=Spectroscopy_Old;                        %NOTE!

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

Frequency_Spectroscopy_Considered_Min_Index=find(Frequency_Spectroscopy_Old>Frequency_Considered_Min,1,'first');
Frequency_Spectroscopy_Considered_Max_Index=find(Frequency_Spectroscopy_Old>Frequency_Considered_Max,1,'first');


Spectrum=interp1(Frequency_Old,Spectrum_Frequency,Frequency);
Spectrum_Reference=interp1(Frequency_Old,Spectrum_Reference_Frequency,Frequency);
Spectroscopy=interp1(Frequency_Spectroscopy_Old,Spectroscopy_Frequency,Frequency);

Spectrum(isnan(Spectrum))=0;
Spectrum(Frequency<Min_Frequency)=0;

Spectrum_Reference(isnan(Spectrum))=0;
Spectrum_Reference(Frequency<Min_Frequency)=0;


Spectroscopy(isnan(Spectroscopy))=0;
Spectroscopy(Frequency<Min_Frequency)=0;


%% Spatial Filtering 

Spectrum((N_f+1):N_t)=0;
Spectrum_Reference((N_f+1):N_t)=0;

Time_total=1/(Max_Frequency/(N_f-1));
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2是因為一來一回
Time=Time';
Position=C*Time;
Position_micron=Position*1E6;

Signal=fft(Spectrum);
Signal_Reference=fft(Spectrum_Reference);

%% Data separation (Spatial Filtering)

Signal_Upper=Signal;
Signal_Lower=Signal;

%Signal_1=Signal;        %DC+self interference
%Signal_2=Signal;        %Upper
%Signal_3=Signal;        %Lower

pixel_1=600;               % the end of 1
pixel_2=1600;%2000               % the end of 2
pixel_3=1098;               %for saperation
%pixel_3=1600;               % the end of 3

%Signal_1((pixel_1+1):(length(Signal_1)-pixel_1))=0;

%Signal((pixel_2+1):(length(Signal)-pixel_2))=0;
Signal(1:pixel_1)=0;
Signal((pixel_2+1):(length(Signal)-pixel_2))=0;
Signal((length(Signal)-pixel_1+1):end)=0;

Signal((pixel_2+1):end)=0;
%Signal_Upper(1:pixel_1)=0;
%Signal_Upper((pixel_3+1):(length(Signal)-pixel_3))=0;
%Signal_Upper((length(Signal)-pixel_1+1):end)=0;

%Signal_Lower(1:pixel_3)=0;
%Signal_Lower((pixel_2+1):(length(Signal)-pixel_2))=0;
%Signal_Lower((length(Signal)-pixel_3+1):end)=0;

Signal_Reference(1:pixel_1)=0;
Signal_Reference((pixel_2+1):end)=0;
%Signal(((length(Signal)-pixel_1)+1):end)=0;

%Signal_3((pixel_3+1):(length(Signal_3)-pixel_3))=0;
%Signal_3(1:pixel_2)=0;
%Signal_3(((length(Signal_3)-pixel_2)+1):end)=0;

%plot(Position,Signal_1,Position,Signal_2,Position,Signal_3);
plot(Position,abs(Signal),Position,abs(Signal_Reference));
%% Again FD

%Spectrum_1=real(ifft(Signal_1));       %can take Real, since the imaginary part should be relatively small, and it came from error
%Spectrum_2=real(ifft(Signal_2));
%Spectrum_3=real(ifft(Signal_3));

Spectrum=(ifft(Signal));
Spectrum=2*(Spectrum(1:N_f));     %2 to compensate the amplitude reduction from Hilbert transform


%Spectrum_Upper=(ifft(Signal_Upper));
%Spectrum_Upper=real(Spectrum_Upper(1:N_f));     %2 to compensate the amplitude reduction from Hilbert transform

%Spectrum_Lower=(ifft(Signal_Lower));
%Spectrum_Lower=real(Spectrum_Lower(1:N_f));     %2 to compensate the amplitude reduction from Hilbert transform

Spectrum_Reference=(ifft(Signal_Reference));
Spectrum_Reference=2*(Spectrum_Reference(1:N_f));     %2 to compensate the amplitude reduction from Hilbert transform

%Spectrum_1=Spectrum_1(1:N_f);
%Spectrum_2=Spectrum_2(1:N_f);
%Spectrum_3=Spectrum_3(1:N_f);

%% Divided by Spectrum_Reference

[maxvalue maxindex]=max(abs(Signal));
Thickness_0=Position(maxindex);
    
    %Thickness:the peak position of the second interface
[maxvalue maxindex]=max(abs(Signal_Reference));
Thickness=(Position(maxindex)-Thickness_0);

%Thickness=2.468E-6;
Spectrum=Spectrum./Spectrum_Reference;

%Spectrum_Ratio=Spectrum_Lower./Spectrum_Upper;
%Spectrum=Spectrum-Spectrum_Reference;
%Spectrum=Spectrum-Spectrum_Sample;

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
Spectrum_Original=Spectrum_Reference./(r_BK7).^2;

%% Fitting

n_final(1:length(Frequency),1:length(n_should))=1.5;
k_final(1:length(Frequency),1:length(n_should))=0;
Spectroscopy_Amplitude_Check(1:length(Spectrum),1:length(n_should))=0;
T_check(1:length(Spectrum),1:length(n_should))=0;

for q=1:length(n_should)
    
value_temp=10000000000;
n(1:length(Frequency),1)=1.5;
k(1:length(Frequency),1)=0;

%T_max_check(1:length(Spectrum),1:length(n_should))=0;

value_total_1(1:Number_of_Loop)=0;
value_total_2(1:Number_of_Loop)=0;

    n_o=(n_should(q)-delta_n):0.00001:(n_should(q)+delta_n);
    k_o=-0.05:0.00001:0.1;
    %k_o=k_o';
    delta_n_number=round(delta_n/(n_o(2)-n_o(1)))/2;
    %n_empty=n_o;
    %n_empty(:)=1;
    %k_empty=k_o;
    %k_empty(:)=1;
    %n_temp=k_empty*n_o;
    %k_temp=k_o*n_empty;
    n_temp=n_o;
    k_temp=k_o;
    
    n_should_index=find(n_o-n_should(q)>0,1,'first');
    
for p=1:Number_of_Loop
    
    % n Template Generation

    % For OCT    
% For Spectroscopy    
    Spectrum_Amplitude_Check(1:length(Spectrum))=0;
    for j=1:(Frequency_Considered_Max_Index-Frequency_Considered_Min_Index+1)
            
        if j <= ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)
            index_now=Frequency_Center_Index+j-1;
        elseif j > ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)           %direction reverse!!
            index_now=-j+((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)+Frequency_Center_Index;
        end
            
            % About the conditions
        
        r1=(n1-(n(index_now)+i*k_temp))./(n(index_now)+n1+i*k_temp);
        r1_r=((n(index_now)+i*k_temp)-n1)./(n(index_now)+n1+i*k_temp);    
        t1=2*(n1)./(n(index_now)+n1+i*k_temp);
        t1_r=2*(n(index_now)+i*k_temp)./(n(index_now)+n1+i*k_temp);
        t2=2*(n(index_now)+i*k_temp)./(n(index_now)+n2(index_now)+i*k_temp);
        r2=((n(index_now)+i*k_temp)-n2(index_now))./((n(index_now)+i*k_temp)+n2(index_now));   
        d=exp(i*2*pi.*Frequency(index_now)/C.*(n(index_now)+i*k_temp).*Thickness);   %注意! -1*n!
        d0=exp(i*2*pi.*Frequency(index_now)/C.*Thickness_0);   %注意! -1*n!
        dref=exp(i*2*pi.*Frequency(index_now)/C.*(Thickness_0+Thickness));   %注意! -1*n!
        
        Spectrum_Amplitude_Temp=abs(((r1.*(d0.^2))+t1.*t1_r.*r2.*((d0.*d).^2))./(r_BK7(index_now).*(dref.^2)));                       %神說是cosine
        %Spectroscopy_Temp=abs((t_AN100(index_now).^2).*(t1.*t2.*d)./(1+r1_r.*r2.*(d.^2))).^2;
        %T_Temp=abs((t_AN100(index_now).^2).*(t1.*t2.*d)./(1+r1_r.*r2.*(d.^2))).^2;
        Merit=((Spectrum_Amplitude_Temp-abs(Spectrum(index_now))).^2); %(abs(abs(Spectrum_Temp)-abs(Spectrum(index_now))).^2)+
          

        [value index]=min(Merit);
        value_total_1(p)=value_total_1(p)+value;
        Spectrum_Amplitude_Check(index_now)=Spectrum_Amplitude_Temp(index);
        k(index_now)=k_temp(index);
            
    end
    

    index_old=0;
    Spectrum_Phase_Check(1:length(Spectrum))=0;
    d0_Check(1:length(Spectrum))=0;
    dref_Check(1:length(Spectrum))=0;
    for j=1:(Frequency_Considered_Max_Index-Frequency_Considered_Min_Index+1)
            
        if j <= ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)
            index_now=Frequency_Center_Index+j-1;
        elseif j > ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)           %direction reverse!!
            index_now=-j+((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)+Frequency_Center_Index;
        end
            
        r1=(n1-(n_temp+i*k(index_now)))./(n_temp+n1+i*k(index_now));
        r1_r=((n_temp+i*k(index_now))-n1)./(n_temp+n1+i*k(index_now));    
        t1=2*(n1)./(n_temp+n1+i*k(index_now));
        t1_r=2*(n_temp+i*k(index_now))./(n_temp+n1+i*k(index_now));
        t2=2*(n_temp+i*k(index_now))./(n_temp+n2(index_now)+i*k(index_now));
        r2=((n_temp+i*k(index_now))-n2(index_now))./((n_temp+i*k(index_now))+n2(index_now));   
        d=exp(i*2*pi.*Frequency(index_now)/C.*(n_temp+i*k(index_now)).*Thickness);   %注意! -1*n!
        d0=exp(i*2*pi.*Frequency(index_now)/C.*Thickness_0);   %注意! -1*n!
        dref=exp(i*2*pi.*Frequency(index_now)/C.*(Thickness_0+Thickness));   %注意! -1*n!
        
        Spectrum_Phase_Temp=angle(((r1.*(d0.^2))+t1.*t1_r.*r2.*((d0.*d).^2))./(r_BK7(index_now).*(dref.^2)));                       %神說是cosine
        
        Merit=(Spectrum_Phase_Temp-angle(Spectrum(index_now))).^2; %(abs(abs(Spectrum_Temp)-abs(Spectrum(index_now))).^2)+
          
        if range_specified == 1
            if (index_now == Frequency_Center_Index) || (index_now == (Frequency_Center_Index -1))
                range_upper=min(n_should_index+delta_n_number,length(Merit));
                range_lower=max(n_should_index-delta_n_number,1);
            elseif (index_now > Frequency_Center_Index) || (index_now < (Frequency_Center_Index -1))
                range_upper=min(index_old+delta_n_number,length(Merit));
                range_lower=max(index_old-delta_n_number,1);
            end
            [value index]=min(Merit(range_lower:range_upper));
            index=range_lower+index-1;
            index_old=index;
        else
            [value index]=min(Merit);
        end
        
        value_total_2(p)=value_total_2(p)+value;
        Spectrum_Phase_Check(index_now)=Spectrum_Phase_Temp(index);
        dref_Check(index_now)=dref;
        d0_Check(index_now)=d0;
        n(index_now)=n_temp(index);    %index1,index2=k,n  in situ saving all the solutions
            
    end
    
    
end

    n_final(:,q)=n;
    k_final(:,q)=k;

end

%% Checking

plot(Wavelength_micron,abs(Spectrum),Wavelength_micron,Spectrum_Amplitude_Check);

plot(Wavelength_micron,angle(Spectrum),Wavelength_micron,Spectrum_Phase_Check);
%plot(Wavelength_micron,Spectrum_Check,Wavelength_micron,Spectrum);

%plot(Wavelength_micron,Spectroscopy_Check,Wavelength_micron,Spectroscopy);
plot(Wavelength_micron,n_final);
plot(Wavelength_micron,k_final);
plot(1:Number_of_Loop,value_total_1,1:Number_of_Loop,value_total_2);
dlmwrite('Thickness_0.txt',Thickness_0,'delimiter','\t','newline','pc','precision','%.12f');
dlmwrite('Thickness.txt',Thickness,'delimiter','\t','newline','pc','precision','%.12f');
