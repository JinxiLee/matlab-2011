clear all

%% Setting

Load_Previous_Result=1;

n_should=1.6;
Wavelength_Center=550;

Number_of_loop=100;

Wavelength_Considered_Min=540;          %nm
Wavelength_Considered_Max=560;


Max_Wavelength=800;             %nm
Min_Wavelength=300;             %nm
N_f=4096;
N_t=4096*2;

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

Frequency_Center_Index=find(Frequency>Frequency_Center,1,'first');

Frequency_Considered_Min_Index=find(Frequency>Frequency_Considered_Min,1,'first');
Frequency_Considered_Max_Index=find(Frequency>Frequency_Considered_Max,1,'first');

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
Time=[0:Time_total/(N_t-1):Time_total]/2;%/2琌ㄓ
Time=Time';
Position=C*Time;
Position_micron=Position*1E6;

Signal=fft(Spectrum);

%% Data separation (Spatial Filtering)

Signal_1=Signal;        %DC+self interference
Signal_2=Signal;        %Upper
Signal_3=Signal;        %Lower

pixel_1=150;               % the end of 1
pixel_2=323;               % the end of 2
pixel_3=420;               % the end of 3

Signal_1((pixel_1+1):(length(Signal_1)-pixel_1))=0;


Signal_1((pixel_1+1):(length(Signal_1)-pixel_1))=0;

Signal_2((pixel_2+1):(length(Signal_2)-pixel_2))=0;
Signal_2(1:pixel_1)=0;
Signal_2((length(Signal_2)-pixel_1+1):end)=0;

Signal_3((pixel_3+1):(length(Signal_3)-pixel_3))=0;
Signal_3(1:pixel_2)=0;
Signal_3((length(Signal_3)-pixel_2+1):end)=0;
plot(Position,real(Signal_1),Position,real(Signal_2),Position,real(Signal_3));

%% Set the negative space to zero - Noe There is phase
Signal_1((round(length(Signal_1)/2)+1):end)=0;
Signal_2((round(length(Signal_2)/2)+1):end)=0;
Signal_3((round(length(Signal_3)/2)+1):end)=0;
%% Again FD

Spectrum_1=(ifft(Signal_1));       %can take Real, since the imaginary part should be relatively small, and it came from error
Spectrum_2=(ifft(Signal_2));
Spectrum_3=(ifft(Signal_3));

Spectrum_1=Spectrum_1(1:N_f);
Spectrum_2=Spectrum_2(1:N_f);
Spectrum_3=Spectrum_3(1:N_f);
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
Spectrum_Original=Spectrum_Reference./(r_BK7).^2;

%% Overall Variables




%% Fitting

value_temp=10000000000;
n_final(1:length(Frequency),1:length(n_should))=1.5;
k_final(1:length(Frequency),1:length(n_should))=0;
n(1:length(Frequency),1)=1.5;
k(1:length(Frequency),1)=0;


for p=1:length(n_should)                                                        % p: Loop 1, for different solutions
    
    % n Template Generation
    n_o=(n_should(p)-0.3):0.001:(n_should(p)+0.3);
    k_o=0:0.001:0.2;
    k_o=k_o';

    n_empty=n_o;
    n_empty(:)=1;
    k_empty=k_o;
    k_empty(:)=1;
    n_temp=k_empty*n_o;
    k_temp=k_o*n_empty;
    
    n_should_index=find(n_o-n_should(p)>0,1,'first');
    
% initial guesses:

    Current_Scanning_Parameter=1;   %1: Thickness_0, 2: Thickness, 3: Efficiency_1, 4: Efficiency_2, 5: Efficiency_3
    Direction=1;                    %1: for positive, -1: for negative

    Thickness_Pitch=0.1E-6;
    Efficiency_Pitch=0.05;

    if Load_Previous_Result == 0

    %Thickness_0: the peak position of the first interface
        [maxvalue maxindex]=max(abs(Signal_2));
        Thickness_0=Position(maxindex);
        Thickness_0=Thickness_0-20*Thickness_Pitch;
    
    %Thickness:the peak position of the second interface
        [maxvalue maxindex]=max(abs(Signal_3));
        Thickness=(Position(maxindex)-Thickness_0)/n_should(p);
        Thickness=Thickness-20*Thickness_Pitch;
    
    %Efficiency: 1
        %Efficiency_1=0.8;
        Efficiency_2=0.5;
        Efficiency_3=0.5;
        %Efficiency_1=Efficiency_1-10*Efficiency_Pitch;
        %Efficiency_2=Efficiency_2-10*Efficiency_Pitch;
        %Efficiency_3=Efficiency_3-10*Efficiency_Pitch;
        
    elseif Load_Previous_Result == 1
        Thickness_0=importdata('Thickness_0.txt');
        Thickness=importdata('Thickness.txt');
        Efficiency_1=importdata('Efficiency_1.txt');
        Efficiency_2=importdata('Efficiency_2.txt');
        Efficiency_3=importdata('Efficiency_3.txt');
    end
    
    q=1;
	Error_Difference=-1;    % Assumed conv
    Error_Total(1:Number_of_loop)=0;
    Error_Total_1(1:Number_of_loop)=0;
    Error_Total_2(1:Number_of_loop)=0;
    Error_Total_3(1:Number_of_loop)=0;
    Direction_Record(1:Number_of_loop)=0;
    Scanning_Parameter_Record(1:Number_of_loop)=0;
    while (q <= Number_of_loop)                                                    % q: Loop 2, for Efficiency_1, Efficiency_2, Efficiency_3, Thickness, Thickness_0
        index_2_old=0;
        Error_Total_1_Temp=0;       
        Error_Total_2_Temp=0;
        Error_Total_3_Temp=0;
        Error_Total_2_2_Temp=0;
        Error_Total_3_3_Temp=0;
        for j=1:(Frequency_Considered_Max_Index-Frequency_Considered_Min_Index+1)
            
            if j <= ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)
                index_now=Frequency_Center_Index+j-1;
            elseif j > ((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)           %direction reverse!!
                index_now=-j+((Frequency_Considered_Max_Index-Frequency_Center_Index)+1)+Frequency_Center_Index;
            end
            
            % About the conditions
            % Operand_1
            %Spectrum_1_Temp=(abs((n1-(n_temp+i*k_temp))./(n_temp+n1+i*k_temp)).^2).*(Spectrum_Original(index_now))+(abs(((n_temp+i*k_temp)-n2(index_now))./((n_temp+i*k_temp)+n2(index_now)).*(exp(i*2*pi.*Frequency(index_now)/C.*(-1*n_temp+i*k_temp).*Thickness).^2)).^2).*(Spectrum_Original(index_now))+2*Efficiency_1*real(Spectrum_Original(index_now).*(n1-(n_temp+i*k_temp))./(n_temp+n1+i*k_temp).*(n1)./(n_temp+n1+i*k_temp).*(n_temp+i*k_temp)./(n_temp+n1+i*k_temp).*((n_temp+i*k_temp)-n2(index_now))./((n_temp+i*k_temp)+n2(index_now)).*(exp(i*2*pi.*Frequency(index_now)/C.*(-1*n_temp+i*k_temp).*Thickness).^2));
            % Operand_2
            Spectrum_2_Temp_Env=abs(2*Efficiency_2*(Spectrum_Original(index_now).*r_BK7(index_now).*(n1-(n_temp+i*k_temp))./(n_temp+n1+i*k_temp).*exp(i*2*pi.*Frequency(index_now)/C.*Thickness_0).^2));                       %弧琌cosine
            Spectrum_2_Temp_Phase=angle(2*Efficiency_2*(Spectrum_Original(index_now).*r_BK7(index_now).*(n1-(n_temp+i*k_temp))./(n_temp+n1+i*k_temp).*exp(i*2*pi.*Frequency(index_now)/C.*Thickness_0).^2));                       %弧琌cosine
            % Operand_3
            Spectrum_3_Temp_Env=abs(2*Efficiency_2*(Spectrum_Original(index_now).*r_BK7(index_now).*(n1)./(n_temp+n1+i*k_temp).*(n_temp+i*k_temp)./(n_temp+n1+i*k_temp).*((n_temp+i*k_temp)-n2(index_now))./((n_temp+i*k_temp)+n2(index_now)).*(exp(i*2*pi.*Frequency(index_now)/C.*Thickness_0).*exp(i*2*pi.*Frequency(index_now)/C.*(-1*n_temp+i*k_temp).*Thickness)).^2));
            Spectrum_3_Temp_Phase=angle(2*Efficiency_2*(Spectrum_Original(index_now).*r_BK7(index_now).*(n1)./(n_temp+n1+i*k_temp).*(n_temp+i*k_temp)./(n_temp+n1+i*k_temp).*((n_temp+i*k_temp)-n2(index_now))./((n_temp+i*k_temp)+n2(index_now)).*(exp(i*2*pi.*Frequency(index_now)/C.*Thickness_0).*exp(i*2*pi.*Frequency(index_now)/C.*(-1*n_temp+i*k_temp).*Thickness)).^2));
            % Merit Function
            %Merit=((Spectrum_1_Temp-Spectrum_1(index_now)).^2)+((Spectrum_2_Temp-Spectrum_2(index_now)).^2)+((Spectrum_3_Temp-Spectrum_3(index_now)).^2);
            Merit=(abs(Spectrum_2_Temp_Env-abs(Spectrum_2(index_now))).^2)+(abs(Spectrum_3_Temp_Env-abs(Spectrum_3(index_now))).^2)+(abs(Spectrum_2_Temp_Phase-angle(Spectrum_2(index_now))).^2)+(abs(Spectrum_3_Temp_Phase-angle(Spectrum_3(index_now))).^2);
            % About the Solution Selection
        
            if (index_now == Frequency_Center_Index) || (index_now == (Frequency_Center_Index -1))
                range_upper=min(n_should_index+40,size(Merit,2));
                range_lower=max(n_should_index-40,1);
            elseif (index_now > Frequency_Center_Index) || (index_now < (Frequency_Center_Index -1))
                range_upper=min(index_2_old+40,size(Merit,2));
                range_lower=max(index_2_old-40,1);
            end
            [value index_1]=min(Merit(:,range_lower:range_upper));
            [value index_2]=min(value);
            index_2_temp=index_2;
            index_2=range_lower+index_2-1;
            index_2_old=index_2;
            

            
            n(index_now)=n_temp(index_1(index_2_temp), index_2);    %index1,index2=k,n  in situ saving all the solutions
            k(index_now)=k_temp(index_1(index_2_temp), index_2);
            
            % Now I am curious about the contribution of each operands
            %Error_Total_1_Temp= Error_Total_1_Temp+((Spectrum_1_Temp(index_1(index_2_temp), index_2)-Spectrum_1(index_now)).^2);  
            Error_Total_2_Temp= Error_Total_2_Temp+((Spectrum_2_Temp_Env(index_1(index_2_temp), index_2)-abs(Spectrum_2(index_now))).^2);
            Error_Total_3_Temp= Error_Total_3_Temp+((Spectrum_3_Temp_Env(index_1(index_2_temp), index_2)-abs(Spectrum_3(index_now))).^2);
            Error_Total_2_2_Temp= Error_Total_2_2_Temp+((Spectrum_2_Temp_Phase(index_1(index_2_temp), index_2)-angle(Spectrum_2(index_now))).^2);
            Error_Total_3_2_Temp= Error_Total_3_3_Temp+((Spectrum_3_Temp_Phase(index_1(index_2_temp), index_2)-angle(Spectrum_3(index_now))).^2);
            %Error_Total_Temp=Error_Total_1_Temp+Error_Total_2_Temp+Error_Total_3_Temp;
            Error_Total_Temp=Error_Total_2_Temp+Error_Total_3_Temp+Error_Total_2_2_Temp+Error_Total_3_2_Temp;
        end
        % Algorithm: 场常タ, scanㄤい, div碞传scan, try–把计常琌div,
        % 碞场常传は, 繷try
        
        %Error_Total_1(q)=Error_Total_1_Temp;
        Error_Total_2(q)=Error_Total_2_Temp;
        Error_Total_3(q)=Error_Total_3_Temp;
        Error_Total(q)=Error_Total_Temp;
        
        Scanning_Parameter_Record(q)=Current_Scanning_Parameter;
        Direction_Record(q)=Direction;
        
        if q == 1
            Error_Difference=-1;
        else
            Error_Difference=Error_Total(q)-Error_Total(q-1);
        end    
        
        if Error_Difference > 0                    % Start to diverge
            
            if Current_Scanning_Parameter == 1              % Reset privious change
                Thickness_0=Thickness_0-Direction*Thickness_Pitch;
            elseif Current_Scanning_Parameter == 2
                Thickness=Thickness-Direction*Thickness_Pitch;
            %elseif Current_Scanning_Parameter == 3
            %    Efficiency_1=Efficiency_1-Direction*Efficiency_Pitch;
            elseif Current_Scanning_Parameter == 3
                Efficiency_2=Efficiency_2-Direction*Efficiency_Pitch;
            %elseif Current_Scanning_Parameter == 4
            %    Efficiency_3=Efficiency_3-Direction*Efficiency_Pitch;
            end
                
            if Current_Scanning_Parameter < 3       %5: for last parameter Efficiency_3
                Current_Scanning_Parameter=Current_Scanning_Parameter+1;
            else
                Direction=Direction*(-1);
                Current_Scanning_Parameter=1;
            end
        end
        
        if Current_Scanning_Parameter == 1              % Reset privious change
            Thickness_0=Thickness_0+Direction*Thickness_Pitch;
        elseif Current_Scanning_Parameter == 2
            Thickness=Thickness+Direction*Thickness_Pitch;
        %elseif Current_Scanning_Parameter == 3
        %    Efficiency_1=Efficiency_1+Direction*Efficiency_Pitch;
        elseif Current_Scanning_Parameter == 3
            Efficiency_2=Efficiency_2+Direction*Efficiency_Pitch;
        %elseif Current_Scanning_Parameter == 4
        %    Efficiency_3=Efficiency_3+Direction*Efficiency_Pitch;
        end  
        q=q+1;
    end
    
    n_final(:,p)=n;
    k_final(:,p)=k;
    
end
plot(1:Number_of_loop,Error_Total_1,1:Number_of_loop,Error_Total_2,1:Number_of_loop,Error_Total_3);

plot(1:Number_of_loop,Scanning_Parameter_Record,1:Number_of_loop,Direction_Record);

%% Checking


t1_check=2*n1./(n_final+n1+i*k_final);
t1_r_check=2*(n_final+i*k_final)./(n_final+n1+i*k_final);           %_1 means reverse direction
r1_check=(n1-(n_final+i*k_final))./(n_final+n1+i*k_final);
r2_check=((n_final+i*k_final)-n2)./(n_final+n2+i*k_final);

d_check=exp(i*2*pi.*Frequency/C.*(-1*n_final+i*k_final).*Thickness);   %猔種! -1*n!
d0_check=exp(i*2*pi.*Frequency/C.*Thickness_0);   %猔種! -1*n!


%Spectrum_1_check=(abs(r1_check).^2).*(Spectrum_Original)+(abs(r2_check.*(d_check.^2)).^2).*(Spectrum_Original)+2*Efficiency_1*real(Spectrum_Original.*r1_check.*t1_check.*t1_r_check.*r2_check.*(d_check.^2));
Spectrum_2_check=2*Efficiency_2*(Spectrum_Original.*r_BK7.*r1_check.*d0_check.^2);                       %弧琌cosine
Spectrum_3_check=2*Efficiency_3*(Spectrum_Original.*r_BK7.*t1_check.*t1_r_check.*r2_check.*(d0_check.*d_check).^2);

plot(Frequency,Spectrum_2,Frequency,Spectrum_2_check);

dlmwrite('Thickness_0.txt',Thickness_0,'delimiter','\t','newline','pc','precision','%.12f');
dlmwrite('Thickness.txt',Thickness,'delimiter','\t','newline','pc','precision','%.12f');
dlmwrite('Efficiency_1.txt',Efficiency_1,'delimiter','\t','newline','pc','precision','%.12f');
dlmwrite('Efficiency_2.txt',Efficiency_2,'delimiter','\t','newline','pc','precision','%.12f');
dlmwrite('Efficiency_3.txt',Efficiency_3,'delimiter','\t','newline','pc','precision','%.12f');
