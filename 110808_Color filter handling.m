clear all;
close all;
clc;

total_OPD=115.8;        %micron
axial_resolution=1.5;   %micron


ratiop=[0:0.1:1];   
ratios=[0:0.1:1];


%ratiop=[1.2:0.01:1.3];   
%ratios=[1.7:0.01:1.8];
 
QQ=0.9:0.01:1.1;

start_position=61;
end_position=67;


Is_o=importdata('D:\110808\110808_20Hz mirror A-scan');
%Is_o=importdata('D:\110808\110808_2Hz mirror A-scan');
Is_o=Is_o(:,1:10);
%Is_o=importdata('D:\110426_Original PSOCT data\PALC_DC_high
%Z_1V_1_s.txt');

Is=sum(Is_o,2)/size(Is_o,2);

position=[0:total_OPD/(length(Is)-1):total_OPD]';  

%plot(position,Is_o);


%Is_background=[Is(1):((Is(length(Is))-Is(1))/(length(Is)-1)):Is(end)]';

%Is=Is-Is_background;

%% Background substration
%Starting_pixel=1000;

Starting_pixel=2370; %GLCD
Ending_pixel=4000;


%Starting_pixel=2000; 
%Ending_pixel=3000;

Is_o(1:Starting_pixel,:)=0;

Is_o(Ending_pixel:end,:)=0;

%% Zero padding and filtering

Starting_pixel_f=400;

Ending_pixel_f=1200;


Is_o_ftemp=fft(Is_o,[],1);

Is_o_ftemp(1:Starting_pixel_f,:)=0;

Is_o_ftemp((length(Is_o_ftemp)-Starting_pixel_f+1):length(Is_o_ftemp),:)=0;

Is_o_ftemp(Ending_pixel_f:(length(Is_o_ftemp)-Ending_pixel_f),:)=0;


Is_o_ftemp1=Is_o_ftemp(1:round(length(Is_o_ftemp)/2),:);
Is_o_ftemp2=Is_o_ftemp((round(length(Is_o_ftemp)/2)+1):end,:);
Is_o_ftemp=[Is_o_ftemp1;zeros(length(Is_o_ftemp)*7,size(Is_o_ftemp,2)); Is_o_ftemp2];
Is_o_new=real(ifft(Is_o_ftemp,[],1));
Is_o_ftemp_for_hilbert=2*Is_o_ftemp;
Is_o_ftemp_no_phase=abs(Is_o_ftemp);
Is_o_ftemp_for_hilbert((round(size(Is_o_ftemp_for_hilbert,1)/2)+1):end,:)=0;
Is_o_new_envelope=abs(ifft(Is_o_ftemp_for_hilbert,[],1));
%Is_o_new=ifftshift(real(ifft(abs(Is_o_ftemp))));
dlmwrite('D:\Is_o_ftemp_no_phase.txt',Is_o_ftemp_no_phase,'delimiter','\t','newline','pc');
position=[0:total_OPD/(length(Is_o_new)-1):total_OPD]';  

plot(position,Is_o_new,position,Is_o_new_envelope);

%% lambda genaration

c=3E8;

dx=(total_OPD*(1E-6))/(length(Is)-1);    %m

dt=2*dx/c;    %this is for real OPD

f_total=4/dt;

frequency=1:f_total/(length(Is_o_new)-1):f_total;

wavelength=c./frequency'*1E6;

plot(wavelength(500:end),abs(Is_o_ftemp(500:end,:)));

%plot(frequency,abs(Is_o_ftemp));

%% use known spectrum

spectrumData=dlmread('D:\Ce.txt');

wavelength_old=spectrumData(:,1)/1000;
spectrumPower_old=spectrumData(:,2);

spectrumPower_new=interp1(wavelength_old,spectrumPower_old,wavelength);
spectrumPower_new(isnan(spectrumPower_new))=0;
spectrumPower_new=spectrumPower_new/max(spectrumPower_new);

dlmwrite('D:\spectrumPower_new.txt',spectrumPower_new,'delimiter','\t','newline','pc');
%plot(position,Is_o_new,position,Is_o_new_envelope);

dlmwrite('D:\wavelength.txt',wavelength,'delimiter','\t','newline','pc');


%% to find max and align

[maxvalue maxindex]=max(Is_o_new,[],1);
meanmaxindex=mean(maxindex);
needshift=round((round(size(Is_o_new,1))/2-maxindex));


for j=1:size(Is_o_new,2)
    Is_o_new(:,j)=circshift(Is_o_new(:,j),needshift(j));
end



%% to find zeros

Starting_pixel_n=38000;

Ending_pixel_n=43000;


for k=1:size(Is_o_new,2)
tag_Is=0;
tag_Ip=0;
x=1;
y=1;
for j=1:(Ending_pixel_n-Starting_pixel_n)
    if tag_Is ==1
        if Is_o_new(Starting_pixel_n-1+j,k)<0
            Zero_Is(x,k)=Starting_pixel_n-1+j;
            x=x+1;
            tag_Is=-1;
        end
    elseif tag_Is == -1
        if Is_o_new(Starting_pixel_n-1+j,k)>0
            Zero_Is(x,k)=Starting_pixel_n-1+j;
            x=x+1;
            tag_Is=1;
        end
    end
        
     
    if tag_Is ==0
        if Is_o_new(Starting_pixel_n-1+j,k)>0
            tag_Is=1;
        else
            tag_Is=-1;
        end
    end
    
end
end

Zero_Is_diff=diff(Zero_Is,1);

Zero_Is_diff(abs(Zero_Is_diff)>1000)=0;

Zero_Is(Zero_Is==0)=max(max(Zero_Is));

plot(position(Zero_Is(1:(size(Zero_Is,1)-1))),1./Zero_Is_diff*total_OPD/(length(Is_o_new)-1));

%% 先只用Zeros來做看看, first to generate the new position array
% to find first and last zero in the range
position_diff=(length(Zero_Is)-1)*0.14/(Zero_Is(end)-Zero_Is(1));
position=0:position_diff:((length(position)-1)*position_diff);
Starting_Zero_Is_pos=position(Zero_Is(1));
Ending_Zero_Is_pos=position(Zero_Is(end));

%% set the diff to 0.135402469462589 micron (acquired from PSF calculation)

set_Zero_diff_Is_pos=(Ending_Zero_Is_pos-Starting_Zero_Is_pos)/(length(Zero_Is)-1);
position_Is_unified=position;

for j=1:(Zero_Is(end)-Zero_Is(1)-1)
        left_Zero_index=Zero_Is(find(Zero_Is<(Zero_Is(1)+j),1,'last'));
        right_Zero_index=Zero_Is(find(Zero_Is>=(Zero_Is(1)+j),1,'first'));
        if (isempty(left_Zero_index)==0)&&(isempty(right_Zero_index)==0)
            small_diff=set_Zero_diff_Is_pos/(right_Zero_index-left_Zero_index);
            position_Is_unified(Zero_Is(1)+j)=position_Is_unified(Zero_Is(1)+j-1)+small_diff;
        end
end

Is_interp=interp1(position,Is,position);



%for j=1:(length(Zero_Ip)-1)
%        Ceter_Ip(j)=(Zero_Ip(j)+Zero_Ip(j+1))/2;
%end

%position=[0:total_OPD/(length(Ip)-1):total_OPD]';  < Old position



%max_window=max(max(Zero_Is_diff),max(Zero_Ip_diff));

%max_window=3;

%window(1:max_window)=1;
%window=window';

%Is_windowed=convn(Is,window,'same');

%Ip_windowed=convn(Ip,window,'same');


[Is_peakvalue Is_peakindex]=max(Is_interp(1:saparation_index));

position_shift=position-position(Is_peakindex);



%% to find first peak (for Is and Ip)

%space: position (micron)
%power: Is_carrier

%% first peak generation

N_f=2^10;
N_t=N_f*4;

spectrumData=dlmread('D:\Ce.txt');

wavelength_old=spectrumData(:,1);
spectrumPower_old=spectrumData(:,2);

spectrumPower_old=spectrumPower_old-mean(spectrumPower_old(1:1000));

c=3E8;                     %m/sec
freq=c./(wavelength_old*1E-9);     %Hz orignal freq array
d_f=max(freq)/(N_f-1);
fx=0:d_f:max(freq);        %freq after interpolation
d_t=1/(d_f*(N_t-1));
%d_t=1/(d_f*2*N_f);
time=[-0.5*(N_t-1)*d_t:d_t:0.5*N_t*d_t]'/2;%/2是因為一來一回
%time=[-(N_f-1)*d_t:d_t:N_f*d_t]/2;  %/2是因為一來一回
space_FD=c*time; % 暫定只用array前1/100 (大約也有100 micron左右吧)
space_FD=(space_FD-space_FD(1))*1E6;                     % shift to zero & m to micron
%plot(lambda,inter,lambda2,inter);

%S0=inter;
%plot(lambda,S0/max(S0),lambda,inter/max(inter),lambda,ref/max(ref));
%S0=S0/max(S0);        %不先normal應該也沒差?
spectrumPower=interp1(freq,spectrumPower_old.*wavelength_old.^2,fx);
spectrumPower(isnan(spectrumPower))=0;


CS=fft(spectrumPower,N_t)';     %with minus time

CS=fftshift(CS);
%CS=real(fft(S_padded));     %with minus time
CS_normal=CS/max(abs(CS));
CS_envelope=abs(CS);    

%CS_carrier=real(CS_normal);

[inter_peakvalue_FD inter_peakindex_FD]=max(CS_envelope);
inter_position_FD=space_FD(inter_peakindex_FD);
FWHM_FD_right=space_FD(find(CS_envelope>0.5*inter_peakvalue_FD, 1, 'last'));
FWHM_FD_left=space_FD(find(CS_envelope>0.5*inter_peakvalue_FD, 1, 'first'));
FWHM_inter_FD=FWHM_FD_right-FWHM_FD_left;
FWHM_FD=FWHM_inter_FD;

space_FD=space_FD-inter_position_FD;

%% false carrier 

False_Carrier=cos(pi*space_FD/0.14);

CS_carrier=CS_envelope.*False_Carrier/max(CS_envelope);


PSF_zero_index=find(CS_carrier(inter_peakindex_FD:end)<0,1,'first');


CS_interp=interp1(space_FD,CS_carrier,position_shift);

for j=1:length(CS_interp)
    if isnan(CS_interp(j))
        CS_interp(j)=0;
    end
end
%% first peak substraction 

CS_real=real(CS_interp);

tempmaxs=max(Is_interp);
goodShifts=1;
goodRatios=0;

for j=1:length(ratios)
    for k=1:20   
        Is_sub=Is_interp-circshift(CS_real,[0 10-k])/max(CS_real)*Is_peakvalue*ratios(j);  %this 0.9 can be critical...

        if max(Is_sub(1:saparation_index))<tempmaxs
            tempmaxs=max(Is_sub(1:saparation_index));
            goodRatios=ratios(j);
            goodShifts=k;
        end
    end
end

%% to use good ratio

        Is_sub=Is_interp;%-circshift(CS_real,[0 10-goodShifts])/max(CS_real)*Is_peakvalue*goodRatios;  %this 0.9 can be critical...
        Subs=circshift(CS_real,[0 5-goodShifts])/max(CS_real)*Is_peakvalue*goodRatios;
        %plot(position,Is_interp,position,circshift(CS_real,[0 5-goodShifts])/max(CS_real)*Is_peakvalue*goodRatios);
   

%% filtering and manual hilbert

Is_f_o=(fft(Is_sub));
Is_f=Is_f_o;

start_index_of_spectrum=200;

end_index_of_spectrum=550;

Is_f(end_index_of_spectrum:end)=0;


Is_f(1:start_index_of_spectrum)=0;



Is_hilberted=abs(ifft(Is_f));

Is_hilberted_real=real(ifft(Is_f));

Is_carrier=ifft(Is_f);

angs=angle(ifft(Is_f));

position_new=[0:total_OPD/(length(Is_hilberted)-1):total_OPD]';  

%% flip sub

%Is_sub=Is_hilberted;
%Ip_sub=Ip_hilberted;

%if (2*Is_peakindex+1)<=length(Is)
%Is_hilberted((Is_peakindex+1):(2*Is_peakindex+1))=Is_hilberted((Is_peakindex+1):(2*Is_peakindex+1))-flipdim(Is_hilberted(1:(Is_peakindex)),1);

%Ip_hilberted((Is_peakindex+1):(2*Is_peakindex+1))=Ip_hilberted((Is_peakindex+1):(2*Is_peakindex+1))-flipdim(Ip_hilberted(1:(Is_peakindex)),1);
%else
%Is_sub((Is_peakindex+1):end)=Is_hilberted((Is_peakindex+1):end)-flipdim(Is_hilberted((Is_peakindex+1-(length(Is)-(Is_peakindex))):(Is_peakindex)),1);

%Ip_sub((Is_peakindex+1):end)=Ip_hilberted((Is_peakindex+1):end)-flipdim(Ip_hilberted((Is_peakindex+1-(length(Is)-(Is_peakindex))):(Is_peakindex)),1);
%end

%% phase (and to filter


Is_3=Is_interp*1E3;

%% Phase handling

for j=1:length(angs)
    if j<length(angs)
    if angs(j+1) > angs(j) + 0.5*2*pi
        angs(j+1:end)=angs(j+1:end)-2*pi;
    elseif angs(j+1) < angs(j) - 0.5*2*pi       
        angs(j+1:end)=angs(j+1:end)+2*pi;
    end
    
end

%plot(position_new,Is_hilberted,position_new,Ip_hilberted,position_new,phase);

%plot(position_new,phase);

%plot(position_new(11200:11600),Is_interp(11200:11600),position_new(11200:11600),Subs(11200:11600),position_new(11200:11600),Is_sub(11200:11600));
%plot(position_new(11200:11680),Is_hilberted_real(11200:11680),position_new(11200:11680),Is_hilberted(11200:11680),position_new(11200:11680),Ip_hilberted_real(11200:11680),position_new(11200:11680),Ip_hilberted(11200:11680));
%plot(position_new(11200:11680),phase(11200:11680));

%dlmwrite('D:\position.txt',position,'delimiter','\t','newline','pc');
%dlmwrite('D:\position_new.txt',position_new,'delimiter','\t','newline','pc');
%dlmwrite('D:\Is.txt',Is,'delimiter','\t','newline','pc');
%dlmwrite('D:\Ip.txt',Ip,'delimiter','\t','newline','pc');
%dlmwrite('D:\Is_hilberted.txt',Is_hilberted,'delimiter','\t','newline','pc');
%dlmwrite('D:\Ip_hilberted.txt',Ip_hilberted,'delimiter','\t','newline','pc');
%dlmwrite('D:\Is_hilberted_real.txt',Is_hilberted_real,'delimiter','\t','newline','pc');
%dlmwrite('D:\Ip_hilberted_real.txt',Ip_hilberted_real,'delimiter','\t','newline','pc');
%dlmwrite('D:\Angs.txt',angs,'delimiter','\t','newline','pc');
%dlmwrite('D:\Angp.txt',angp,'delimiter','\t','newline','pc');
%dlmwrite('D:\Ang_diff.txt',ang_diff,'delimiter','\t','newline','pc');

%plot(position,real(Is_carrier),position,real(Is_carrier_sub),position,CS_real); figure(gcf)


%plot(position,Is_carrier_sub,position,Ip_carrier_sub); figure(gcf)
