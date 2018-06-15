clear all;
close all;
clc;


n_o=1.3:0.001:1.8;

k_o=0:0.001:0.1;
k_o=k_o';

thickness_temp=2.68E-6;

n_should=1.58:1.59;


A=0:0.01:2; %a varying scaling factor

total_OPD=153.4;

range_specified=1;

%% Data pre-loading for Reference selection


cd('D:\110915\2Hz - Green 2\');     %4 (2Hz)
Signal=importdata('1');

Signal=Signal(:,1);
Signal(4500:end)=Signal(4500);

[signal_maxvalue signal_maxindex]=max(Signal,[],1);


cd('D:\110915\2Hz - after Green 2\');     %4 (2Hz)

ref_sig_index_dif=10000;

req_ref_kk=1;

for kk=1:50
    
    
    Reference=importdata(sprintf('%i',kk));
    
    Reference(4500:end)=Reference(4500);
    
    [signal_maxvalue ref_maxindex]=max(Reference,[],1);
    
    if abs(signal_maxindex-ref_maxindex)
        ref_sig_index_dif=abs(signal_maxindex-ref_maxindex);
        req_ref_kk=kk;
    end
    
end


Reference=importdata(sprintf('%i',req_ref_kk));

Reference(4500:end)=Reference(4500);

[maxvalue maxindex]=max(Reference,[],1);
needshift=-maxindex;

Reference=circshift(Reference,needshift);

pixel_1=131;               % the end of upper (the larger peak)
pixel_2=350;               % the end of lower (the larger peak)
pixel_3=9870;               % the start of upper (the larger peak)


Starting_pixel_f=600;
Ending_pixel_f=1800;


%Starting_pixel_f_considered=1000;
%Ending_pixel_f_considered=1300;

Starting_pixel_f_considered=600;
Ending_pixel_f_considered=1800;


Reference_temp_f=fft(Reference);

Reference_temp_f(1:Starting_pixel_f)=0;
Reference_temp_f((length(Reference_temp_f)-Starting_pixel_f+1):length(Reference_temp_f),:)=0;

Reference=real(ifft(Reference_temp_f));

Reference(pixel_2:pixel_3,:)=0;

center_wavelength=0.56; %micron

c=3E8;

center_frequency=c/(center_wavelength*1E-6);

Reference_f=fft(Reference,[],1);
Reference_f(1:Starting_pixel_f,:)=0;
Reference_f(Ending_pixel_f:end,:)=0;


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

%plot(wavelength(500:end),A_exp_abs(500:end));


center_frequency_index=find(frequency-center_frequency>0,1,'first');

%dlmwrite('D:\frequency.txt',frequency,'delimiter','\t','newline','pc');
%dlmwrite('D:\Signal_Upper_f_amp.txt',Signal_Upper_f_amp,'delimiter','\t','newline','pc');
%dlmwrite('D:\Signal_Lower_f_amp.txt',Signal_Lower_f_amp,'delimiter','\t','newline','pc');


%%
%%      Theoritical Calculation and fitting
%%




%% Variables

n_empty=n_o;
n_empty(:)=1;
k_empty=k_o;
k_empty(:)=1;
n_temp=k_empty*n_o;
k_temp=k_o*n_empty;

n1=1;
n2(1:length(frequency),1)=1.5;        %n2 is Bk7 here

%n2 - 1 = C1λ2/(λ2-C2) + C3λ2/(λ2-C4) + C5λ2/(λ2-C6)


C1 = 1.03961212; 
C2 = 0.00600069867; 
C3 = 0.231792344; 
C4 = 0.0200179144; 
C5 = 1.01046945; 
C6 = 103.560653;


wavelength_micron=wavelength*1E6;

wavelength_index_400nm=find(wavelength_micron<0.4,1,'first');

wavelength_index_700nm=find(wavelength_micron<0.7,1,'first');

n_bk7=(C1*(wavelength_micron.^2)./((wavelength_micron.^2)-C2)+C3*(wavelength_micron.^2)./((wavelength_micron.^2)-C4)+C5*(wavelength_micron.^2)./((wavelength_micron.^2)-C6)+1).^0.5;


n2(Starting_pixel_f_considered:Ending_pixel_f_considered,1)=n_bk7(Starting_pixel_f_considered:Ending_pixel_f_considered);


%% Defination

r1=(n1-(n_temp+i*k_temp))./(n_temp+n1+i*k_temp);

n(1:length(frequency),1:length(thickness_temp))=1.75;        %df is the same as frequency, but only take a portion of it (1000*df)
k(1:length(frequency),1:length(thickness_temp))=0;



for jj=1:1

%% 在n上似乎有ambiguity (~0.1)


%%
%%      The First part: A
%%


%% Data Loading


cd('D:\110915\2Hz - Green 2\');     %4 (2Hz)
Signal=importdata(sprintf('%i',jj));

Signal=Signal(:,1);
Signal(4500:end)=Signal(4500);
position=[0:total_OPD/(length(Signal)-1):total_OPD]';  

%% Pre-shifting (以後shift還是shift到0點附近, 雖然這樣定範圍比較麻煩, 但phase比較單純, 不會跟想像中異號)
%% (結果改了之後還是一樣 不知道為什麼) (可能跟MatLabs的fft定義有關, 接下來直接把本程式的d中的phase加個負號)

[maxvalue maxindex]=max(Signal,[],1);
needshift=-maxindex;

Signal=circshift(Signal,needshift);

%% DC filtering (to avoid gap bwtween set zero and raw-data DC)

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

%OPD=abs(Upper_max_index-Lower_max_index)*total_OPD/(length(Signal)-1)*1E-6;

%% To frequency domain

Signal_Upper_f=fft(Signal_Upper,[],1);
Signal_Lower_f=fft(Signal_Lower,[],1);

%Signal_Upper_f(1:Starting_pixel_f,:)=0;
%Signal_Upper_f((length(Signal_Upper_f)-Starting_pixel_f+1):length(Signal_Upper_f),:)=0;
%Signal_Upper_f(Ending_pixel_f:(length(Signal_Upper_f)-Ending_pixel_f),:)=0;
Signal_Upper_f(1:Starting_pixel_f,:)=0;
Signal_Upper_f(Ending_pixel_f:end,:)=0;


%Signal_Lower_f(1:Starting_pixel_f,:)=0;
%Signal_Lower_f((length(Signal_Lower_f)-Starting_pixel_f+1):length(Signal_Lower_f),:)=0;
%Signal_Lower_f(Ending_pixel_f:(length(Signal_Lower_f)-Ending_pixel_f),:)=0;
Signal_Lower_f(1:Starting_pixel_f,:)=0;
Signal_Lower_f(Ending_pixel_f:end,:)=0;

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
value_temp=1000000000000;
n_should_finalindex=1;
A_finalindex=0;
for q=1:length(n_should)
    n_should_index=find(n_o-n_should(q)>0,1,'first');
    value_total=0;
    %n_should=OPD/thickness_temp(q);
    %n_should=OPD*1.5/thickness_temp(q);
    
    index_2_old=0;
    index_2_center=0;
    
    for p=1:(Ending_pixel_f_considered-Starting_pixel_f_considered+1)            %start from center_frequency_index, to Ending_pixel_f, back to center_frequency_index, finally Starting_pixel_f_considered
        if p <= ((Ending_pixel_f_considered-center_frequency_index)+1)
          index_now=center_frequency_index+p-1;
        elseif p > ((Ending_pixel_f_considered-center_frequency_index)+1)           %direction reverse!!
          index_now=-p+((Ending_pixel_f_considered-center_frequency_index)+1)+center_frequency_index;
        end
        
        r2=((n_temp+i*k_temp)-n2(index_now))./(n_temp+n2(index_now)+i*k_temp);    
        d=exp(i*2*pi.*frequency(index_now)/c.*(-1*n_temp+i*k_temp).*thickness_temp);   %注意! -1*n!
        A_temp=(d.^2).*r2./r1;
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
          if index_now == center_frequency_index          
                    [value index_1]=min(DD(:,n_should_index-20:n_should_index+20));

                [value index_2]=min(value);
                index_2_temp=index_2;
                index_2=n_should_index-20+index_2-1;
                index_2_old=index_2;
                index_2_center=index_2;
            elseif index_now > center_frequency_index
                if (index_2_old >= 21) && (index_2_old < size(DD,2)-20)
                    [value index_1]=min(DD(:,index_2_old-20:index_2_old+20));
                    [value index_2]=min(value);
                  index_2_temp=index_2;
                    index_2=index_2_old-20+index_2-1;
                   index_2_old=index_2;
                else
                    [value index_1]=min(DD);
                    [value index_2]=min(value);
                    index_2_temp=index_2;
                    index_2_old=index_2;
                end
            elseif index_now == (center_frequency_index -1)
                if (index_2_center >= 21) && (index_2_center < size(DD,2)-20)
                    [value index_1]=min(DD(:,index_2_center-20:index_2_center+20));
                    [value index_2]=min(value);
                    index_2_temp=index_2;
                    index_2=index_2_center-20+index_2-1;
                    index_2_old=index_2;
                else
                    [value index_1]=min(DD);
                    [value index_2]=min(value);
                    index_2_temp=index_2;
                    index_2_old=index_2;
                end
            elseif index_now < (center_frequency_index -1)
                if (index_2_old >= 21) && (index_2_old < size(DD,2)-20)
                    [value index_1]=min(DD(:,index_2_old-20:index_2_old+20));
                    [value index_2]=min(value);
                  index_2_temp=index_2;
                    index_2=index_2_old-20+index_2-1;
                   index_2_old=index_2;
                else
                    [value index_1]=min(DD);
                    [value index_2]=min(value);
                    index_2_temp=index_2;
                    index_2_old=index_2;
                end
          end
        else
                    [value index_1]=min(DD);
                    [value index_2]=min(value);
                    index_2_temp=index_2;
                    index_2_old=index_2;
        end
        n(index_now,q)=n_temp(index_1(index_2_temp), index_2);    %index1,index2=k,n
        k(index_now,q)=k_temp(index_1(index_2_temp), index_2);
    end
    
    %% to compare with the original data (either in TD or SD)
    
    r1_calculated=(n1-(n(:,q)+i*k(:,q)))./(n(:,q)+n1+i*k(:,q));
    r2_calculated=((n(:,q)+i*k(:,q))-n2)./(n(:,q)+n2+i*k(:,q));
    d_calculated=exp(i*2*pi.*frequency/c.*(-1*n(:,q)+i*k(:,q)).*thickness_temp); 
    
    Signal_calculated=real(ifft(Reference_f.*r1_calculated+Reference_f.*(d_calculated.^2).*r2_calculated));
    
    for rr=1:length(A)
    
        value_total=sum(abs(Signal-A(rr)*Signal_calculated).^2);
    
        if value_total < value_temp
            n_should_finalindex=q;                % the same value for all wavelength!!
            A_finalindex=rr;
            value_temp=value_total;
        end
        
    end
end

n_final=n(:,n_should_finalindex);
k_final=k(:,n_should_finalindex);
n_should_final=n_should(n_should_finalindex);                % the same value for all wavelength!!


%% check


%t1_check=2*n1./(n_final+n1+i*k_final);
%t1_1_check=2*(n_final+i*k_final)./(n_final+n1+i*k_final);           %_1 means reverse direction
%t2_check=2*(n_final+i*k_final)./(n_final+n2+i*k_final);
%t2_1_check=2*n2./(n_final+n2+i*k_final);
%r1_check=(n1-(n_final+i*k_final))./(n_final+n1+i*k_final);
%r1_1_check=((n_final+i*k_final)-n1)./(n_final+n1+i*k_final);
%r2_check=((n_final+i*k_final)-n2)./(n_final+n2+i*k_final);
%r2_1_check=(n2-(n_final+i*k_final))./(n_final+n2+i*k_final);

%d_check=exp(i*2*pi.*frequency/c.*(-1*n_final+i*k_final).*thickness_final); 


%A_check=(d_check.^2).*r2_check./r1_check;

%A_check_abs=abs(A_check);

%A_check_phase=angle(A_check);

%plot(wavelength(500:end),Signal_Upper_f_amp(500:end),wavelength(500:end),Signal_Lower_f_amp(500:end));
%plot(frequency,Signal_Upper_f_amp,frequency,Signal_Lower_f_amp);
%plot(frequency,n_false,frequency,n_final,frequency,k_false,frequency,k_final);

r1_check=(n1-(n_final+i*k_final))./(n_final+n1+i*k_final);
r2_check=((n_final+i*k_final)-n2)./(n_final+n2+i*k_final);
d_check=exp(i*2*pi.*frequency/c.*(-1*n_final+i*k_final).*thickness_temp); 
    
Signal_check=real(ifft(Reference_f.*r1_check+Reference_f.*(d_check.^2).*r2_check));
    

%plot(frequency,A_check_abs,frequency,A_check_phase,frequency,A_false_abs,frequency,A_false_phase);

%plot(frequency,A_check_abs,frequency,A_check_phase,frequency,A_exp_abs,frequency,A_exp_phase);
%plot(wavelength(wavelength_index_700nm:wavelength_index_400nm),Signal_Upper_f_amp(wavelength_index_700nm:wavelength_index_400nm),wavelength(wavelength_index_700nm:wavelength_index_400nm),Signal_Lower_f_amp(wavelength_index_700nm:wavelength_index_400nm));

%plot(wavelength(wavelength_index_700nm:wavelength_index_400nm),n_final(wavelength_index_700nm:wavelength_index_400nm),wavelength(wavelength_index_700nm:wavelength_index_400nm),k_final(wavelength_index_700nm:wavelength_index_400nm));


plot(position,Signal,position,Signal_check);

dlmwrite(sprintf('n_final_%i.txt',jj),n_final,'delimiter','\t','newline','pc');
dlmwrite(sprintf('k_final_%i.txt',jj),k_final,'delimiter','\t','newline','pc');

end
