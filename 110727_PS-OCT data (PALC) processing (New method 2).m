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


Is_o=importdata('D:\PALC\average100\AC1V30Hz_20Hz_116um_100_1\Cut Data\s_Lee.txt');
%Is_o=importdata('D:\110426_Original PSOCT data\PALC_DC_high Z_1V_1_s.txt');

Is_o(1:500)=0;
Is_o(4500:end)=0;

Is=sum(Is_o,2)/size(Is_o,2);

%Is_background=[Is(1):((Is(length(Is))-Is(1))/(length(Is)-1)):Is(end)]';

%Is=Is-Is_background;

Ip_o=importdata('D:\PALC\average100\AC1V30Hz_20Hz_116um_100_1\Cut Data\p_Lee.txt');
%Ip_o=importdata('D:\110426_Original PSOCT data\PALC_DC_high Z_1V_1_p.txt');


Ip_o(1:500)=0;
Ip_o(4500:end)=0;

Ip=sum(Ip_o,2)/size(Ip_o,2);


position=[0:total_OPD/(length(Ip)-1):total_OPD]';  
%% background substraction
Starting_pixel=2500;

Ending_pixel=3300;

saparation_index=2900;

%% to find zeros
tag_Is=0;
tag_Ip=0;
x=1;
y=1;

for j=1:(saparation_index-Starting_pixel)
    if tag_Is ==1
        if Is(Starting_pixel-1+j)<0
            Zero_Is(x)=Starting_pixel-1+j;
            x=x+1;
            tag_Is=-1;
        end
    elseif tag_Is == -1
        if Is(Starting_pixel-1+j)>0
            Zero_Is(x)=Starting_pixel-1+j;
            x=x+1;
            tag_Is=1;
        end
    end
        
    if tag_Ip ==1
        if Ip(Starting_pixel-1+j)<0
            Zero_Ip(y)=Starting_pixel-1+j;
            y=y+1;
            tag_Ip=-1;
        end
    elseif tag_Ip == -1
        if Ip(Starting_pixel-1+j)>0
            Zero_Ip(y)=Starting_pixel-1+j;
            y=y+1;
            tag_Ip=1;
        end
    end
     
    if tag_Is ==0
        if Is(Starting_pixel-1+j)>0
            tag_Is=1;
        else
            tag_Is=-1;
        end
    end
    
    if tag_Ip ==0
        if Ip(Starting_pixel-1+j)>0
            tag_Ip=1;
        else
            tag_Ip=-1;
        end
    end
end

Zero_Is_diff=diff(Zero_Is);
Zero_Ip_diff=diff(Zero_Ip);

%% 先只用Zeros來做看看, first to generate the new position array
% to find first and last zero in the range
position_diff=(length(Zero_Is)-1)*0.13/(Zero_Is(end)-Zero_Is(1));
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

Is_interp=interp1(position_Is_unified,Is,position);

Starting_Zero_Ip_pos=position(Zero_Ip(1));
Ending_Zero_Ip_pos=position(Zero_Ip(end));

set_Zero_diff_Ip_pos=(Ending_Zero_Ip_pos-Starting_Zero_Ip_pos)/(length(Zero_Ip)-1);
position_Ip_unified=position;

for j=1:(Zero_Ip(end)-Zero_Ip(1)-1)
        left_Zero_index=Zero_Ip(find(Zero_Ip<(Zero_Ip(1)+j),1,'last'));
        right_Zero_index=Zero_Ip(find(Zero_Ip>=(Zero_Ip(1)+j),1,'first'));
        if (isempty(left_Zero_index)==0)&&(isempty(right_Zero_index)==0)
            small_diff=set_Zero_diff_Ip_pos/(right_Zero_index-left_Zero_index);
            position_Ip_unified(Zero_Ip(1)+j)=position_Ip_unified(Zero_Ip(1)+j-1)+small_diff;
        end
end

Ip_interp=interp1(position_Ip_unified,Ip,position);


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
[Ip_peakvalue Ip_peakindex]=max(real(Ip_interp(1:Is_peakindex+100)));

Ip_interp=circshift(Ip_interp,Is_peakindex-Ip_peakindex);    %s-pol shift by 200 pixel to minimize total phase retardation 275 for mid min 12 for minmin, 5 for poeriod min

position_shift=position-position(Is_peakindex);



%% to find first peak (for Is and Ip)

%space: position (micron)
%power: Is_carrier

%% first peak generation

N_f=2^12;
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

CS_carrier=real(CS_normal);

[inter_peakvalue_FD inter_peakindex_FD]=max(CS_envelope);
inter_position_FD=space_FD(inter_peakindex_FD);
FWHM_FD_right=space_FD(find(CS_envelope>0.5*inter_peakvalue_FD, 1, 'last'));
FWHM_FD_left=space_FD(find(CS_envelope>0.5*inter_peakvalue_FD, 1, 'first'));
FWHM_inter_FD=FWHM_FD_right-FWHM_FD_left;
FWHM_FD=FWHM_inter_FD;

space_FD=space_FD-inter_position_FD;

PSF_zero_index=find(CS_carrier(inter_peakindex_FD:end)<0,1,'first');


CS_interp=interp1(space_FD,CS_carrier,position_shift);


%% first peak substraction 

CS_real=real(CS_interp);

tempmaxs=max(Is_interp);

for j=1:length(ratios)

Is_sub=Is_interp-CS_real/max(CS_real)*Is_peakvalue*ratios(j);  %this 0.9 can be critical...

if max(Is_sub(1:saparation_index))<tempmaxs
    tempmaxs=max(Is_sub(1:saparation_index));
    goodRatios=ratios(j);
end

end

for j=1:length(ratiop)

Ip_sub=(Ip_interp-CS_real/max(CS_real)*Ip_peakvalue*ratiop(j));  %this 0.9 can be critical...

if max(Ip_sub(1:saparation_index))<tempmaxs
    tempmaxs=max(Ip_sub(1:saparation_index));
    goodRatiop=ratiop(j);
end

end

Is_sub(1:Starting_pixel)=0;

Is_sub(Ending_pixel:end)=0;

Ip_sub(1:Starting_pixel)=0;

Ip_sub(Ending_pixel:end)=0;
%% to use good ratio

%% filtering and manual hilbert

Is_f_o=(fft(Is_sub));
Is_f=Is_f_o;

Ip_f_o=(fft(Ip_sub));
Ip_f=Ip_f_o;

start_index_of_spectrum=200;

end_index_of_spectrum=550;

Is_f(start_index_of_spectrum:(length(Is_f)-start_index_of_spectrum))=0;

Ip_f(start_index_of_spectrum:(length(Ip_f)-start_index_of_spectrum))=0;



Is_hilberted=abs(ifft(Is_f));

Is_hilberted_real=real(ifft(Is_f));

Is_carrier=ifft(Is_f);

angs=angle(ifft(Is_f));

Ip_hilberted=abs(ifft(Ip_f));

Ip_hilberted_real=real(ifft(Ip_f));

Ip_carrier=ifft(Ip_f);

angp=angle(ifft(Ip_f));

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

phase=atan(Ip_hilberted./Is_hilberted)*180/pi;

Is_3=Is_interp*1E3;

Ip_3=Ip_interp*1E3;

%% Phase handling

for j=1:length(angs)
    if j<length(angs)
    if angs(j+1) > angs(j) + 0.5*2*pi
        angs(j+1:end)=angs(j+1:end)-2*pi;
    elseif angs(j+1) < angs(j) - 0.5*2*pi       
        angs(j+1:end)=angs(j+1:end)+2*pi;
    end
    
    if angp(j+1) > angp(j) + 0.5*2*pi
        angp(j+1:end)=angp(j+1:end)-2*pi;
    elseif angp(j+1) < angp(j) - 0.5*2*pi       
        angp(j+1:end)=angp(j+1:end)+2*pi;
    end
    end
end

ang_diff=angs-angp;
%ang_diff(length(ang_diff):length(ang_diff))=0;
ang_diff(length(ang_diff)-30:length(ang_diff))=0;

%plot(position_new,angs,position_new,angp,position_new,ang_diff,position_new,Ip_new*10000000); figure(gcf)


Itotal=(Is_hilberted.^2+Ip_hilberted.^2).^0.5;

plot(position_new,Is_hilberted,position_new,Ip_hilberted,position_new,Itotal);

plot(position_new,phase);

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
