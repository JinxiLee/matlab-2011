clear all

%% Setting


thickness_temp=3.084022018600621E-6;
n_should=1.6:0.1:2.2;

Starting_pixel_f_considered=257;
Ending_pixel_f_considered=309;

reversed=0;

range_specified=1;

Max_Wavelength=800;             %nm
Min_Wavelength=300;             %nm
N_f=8192;
N_t=8192*16;
ROI_ratio=1/16;                  %only consider the first ROI_ratio data in TD

DC_cutoff=800;                 %in TD

%% Data Loading

cd('D:\');
%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
Data=importdata('111031_New Green (400mA 3ms 100ave) near c 3.txt');

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
pixel_2=1300;               % the end of lower (the larger peak)
pixel_3=6400;               % the start of upper (the larger peak)

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

Signal_Upper_f=ifft(Signal_Upper);
Signal_Lower_f=ifft(Signal_Lower);

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


A_exp_diff(1:length(A_exp))=0;
A_exp_diff(2:end)=diff(A_exp);


A_exp(isnan(A_exp))=0;

A_exp_abs(isnan(A_exp_abs))=0;
 
A_exp_phase(isnan(A_exp_phase))=0;

A_exp_diff(isnan(A_exp_diff))=0;

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


%% Fitting

n(1:length(Frequency_Final),1:length(n_should))=1.75;        %df is the same as frequency, but only take a portion of it (1000*df)
k(1:length(Frequency_Final),1:length(n_should))=0;
A(1:length(Frequency_Final),1:length(n_should))=0;
value_temp=10000000000;
for q=1:length(n_should)

    n_o=(n_should(q)-0.2):0.001:(n_should(q)+0.2);
    k_o=-0.1:0.001:0.1;
    k_o=k_o';

    n_empty=n_o;
    n_empty(:)=1;
    k_empty=k_o;
    k_empty(:)=1;
    n_temp=k_empty*n_o;
    k_temp=k_o*n_empty;

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

    %n_should=OPD/thickness_temp(q);
    %n_should=OPD*1.5/thickness_temp(q);
    
    n_should_index=find(n_o-n_should(q)>0,1,'first');
    index_2_old=0;
    index_2_center=0;
    
    for p=1:(Ending_pixel_f_considered-Starting_pixel_f_considered+1)            %start from center_frequency_index, to Ending_pixel_f, back to center_frequency_index, finally Starting_pixel_f_considered
        if p <= ((Ending_pixel_f_considered-center_frequency_index)+1)
          index_now=center_frequency_index+p-1;
        elseif p > ((Ending_pixel_f_considered-center_frequency_index)+1)           %direction reverse!!
          index_now=-p+((Ending_pixel_f_considered-center_frequency_index)+1)+center_frequency_index;
        end
        
        t1=2*(n1)./(n_temp+n1+i*k_temp);
        t1_r=2*(n_temp+i*k_temp)./(n_temp+n1+i*k_temp);
        r1=(n1-(n_temp+i*k_temp))./(n_temp+n1+i*k_temp);
        r2=((n_temp+i*k_temp)-n2(index_now))./((n_temp+i*k_temp)+n2(index_now));    
        d=exp(i*2*pi.*Frequency_Final(index_now)/c.*(n_temp+i*k_temp).*thickness_temp);   %注意! -1*n!
        A_temp=(d.^2).*t1.*t1_r.*r2./r1;
        A_temp_abs=abs(A_temp);
        A_temp_phase=angle(A_temp);
        DD=((A_temp_abs-A_exp_abs(index_now)).^2)+((A_temp_phase-A_exp_phase(index_now)).^2);
        %if (index_2_old > 11) && (index_2_old < size(DD,2)-10)
        %    [value index_1]=min(DD(:,index_2_old-10:index_2_old+10));
        %    [value index_2]=min(value);
        %    index_2_temp=index_2;
        %    index_2=index_2_old-10+index_2-1;
        %    index_2_old=index_2;
        %else
        if range_specified == 1
            if (index_now == center_frequency_index) || (index_now == (center_frequency_index -1))
                range_upper=min(n_should_index+40,size(DD,2));
                range_lower=max(n_should_index-40,1);
            elseif (index_now > center_frequency_index) || (index_now < (center_frequency_index -1))
                range_upper=min(index_2_old+40,size(DD,2));
                range_lower=max(index_2_old-40,1);
            end
            [value index_1]=min(DD(:,range_lower:range_upper));
            [value index_2]=min(value);
            index_2_temp=index_2;
            index_2=range_lower+index_2-1;
            index_2_old=index_2;
        else
            [value index_1]=min(DD);
            [value index_2]=min(value);
            index_2_temp=index_2;
        end
            %index_2_old=index_2;
            %index_2_temp=index_2;
        %end
        n(index_now,q)=n_temp(index_1(index_2_temp), index_2);    %index1,index2=k,n
        k(index_now,q)=k_temp(index_1(index_2_temp), index_2);
        A(index_now,q)=A_temp(index_1(index_2_temp), index_2);    %index1,index2=k,n
        %value_total=value_total+abs(value)^2;
    end
    
    A_diff(1:length(A(:,q)))=0;
    A_diff(2:end)=diff(A(:,q));
    value_total=sum((A_diff-A_exp_diff).^2);
    
    if value_total < value_temp
        finalindex=q;                % the same value for all wavelength!!
        value_temp=value_total;
    end
end

n_final=n(:,finalindex);
k_final=k(:,finalindex);
n_should_final=n_should(finalindex);                % the same value for all wavelength!!


t1_check=2*(n1)./(n_final+n1+i*k_final);
t1_r_check=2*(n_final+i*k_final)./(n_final+n1+i*k_final);
r1_check=(n1-(n_final+i*k_final))./(n_final+n1+i*k_final);
r2_check=((n_final+i*k_final)-n2)./((n_final+i*k_final)+n2);    
d_check=exp(i*2*pi.*Frequency_Final/c.*(n_final+i*k_final).*thickness_temp);   %注意! -1*n!
A_check=(d_check.^2).*t1_check.*t1_r_check.*r2_check./r1_check;

%plot(Position_micron,Signal_Upper,Position_micron,Signal_Lower);

%plot(Wavelength_micron(Wavelength_micron<1),n_final(Wavelength_micron<1),Wavelength_micron(Wavelength_micron<1),k_final(Wavelength_micron<1));
%plot(Wavelength_micron(Wavelength_micron<1),Signal_Upper_f_amp(Wavelength_micron<1),Wavelength_micron(Wavelength_micron<1),Signal_Lower_f_amp(Wavelength_micron<1));
plot(Wavelength_micron(Wavelength_micron<1),n(Wavelength_micron<1,:),Wavelength_micron(Wavelength_micron<1),k(Wavelength_micron<1,:));

%plot(Wavelength_micron(Wavelength_micron<1),A_check(Wavelength_micron<1,:),Wavelength_micron(Wavelength_micron<1),A_exp(Wavelength_micron<1,:));
