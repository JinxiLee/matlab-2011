clear all;
close all;
clc;

total_OPD=115.8;        %micron
axial_resolution=1.5;   %micron


Is_o=importdata('D:\110426_Original PSOCT data\PALC_DC_high Z_1V_1_s.txt');

Is=sum(Is_o,2)/size(Is_o,2);

%Is_background=[Is(1):((Is(length(Is))-Is(1))/(length(Is)-1)):Is(end)]';

%Is=Is-Is_background;

Ip_o=importdata('D:\110426_Original PSOCT data\PALC_DC_high Z_1V_1_p.txt');

Ip=sum(Ip_o,2)/size(Ip_o,2);

%% background substraction
Starting_pixel=2220;

Ending_pixel=2981;

Is_background_slope=(Is(Ending_pixel)-Is(Starting_pixel))/(Ending_pixel-Starting_pixel-1);

Is_background=[Is(Starting_pixel)-Is_background_slope*Starting_pixel:Is_background_slope:Is(Ending_pixel)+Is_background_slope*(length(Is)-Ending_pixel)]';

Is=Is-Is_background;

Ip_background_slope=(Ip(Ending_pixel)-Ip(Starting_pixel))/(Ending_pixel-Starting_pixel-1);

Ip_background=[Ip(Starting_pixel)-Ip_background_slope*Starting_pixel:Ip_background_slope:Ip(Ending_pixel)+Ip_background_slope*(length(Ip)-Ending_pixel)]';

Ip=Ip-Ip_background;

Is(1:Starting_pixel)=0;

Is(Ending_pixel:end)=0;


Ip(1:Starting_pixel)=0;

Ip(Ending_pixel:end)=0;

%% filtering and manual hilbert

Is_f=(fft(Is));

Ip_f=(fft(Ip));

start_index_of_spectrum=200;

end_index_of_spectrum=550;

Is_f(1:start_index_of_spectrum)=0;

Ip_f(1:start_index_of_spectrum)=0;

Is_f(end_index_of_spectrum:end)=0;
Ip_f(end_index_of_spectrum:end)=0;

zeros(1:3*length(Is_f))=0;

Is_f1=Is_f(1:length(Is_f)/2);
Is_f2=Is_f(length(Is_f)/2+1:length(Is_f));


Ip_f1=Ip_f(1:length(Ip_f)/2);
Ip_f2=Ip_f(length(Ip_f)/2+1:length(Ip_f));

Is_f=[Is_f1;zeros';Is_f2];
Ip_f=[Ip_f1;zeros';Ip_f2];

%Is_f(end_index_of_spectrum:length(Is_f)-end_index_of_spectrum)=0;

%Ip_f(end_index_of_spectrum:length(Ip_f)-end_index_of_spectrum)=0;



%Is_f(length(Is_f)-start_index_of_spectrum:length(Is_f))=0;

%Ip_f(length(Is_f)-start_index_of_spectrum:length(Ip_f))=0;

Is_o=Is;

Ip_o=Ip;

Is=abs(ifft(Is_f));

Is_carrier=ifft(Is_f);

angs=angle(ifft(Is_f));


%Is_real=real(ifft(Is_f));
%Is_imag=imag(ifft(Is_f));


Ip=abs(ifft(Ip_f));

Ip_carrier=ifft(Ip_f);

angp=angle(ifft(Ip_f));

%Ip_real=real(ifft(Ip_f));

%Ip_imag=imag(ifft(Ip_f));

Is_f_3=real(Is_f)*1E3;

Ip_f_3=real(Ip_f)*1E3;

%% hilbert transform

%Is_h=abs(hilbert(Is));

%Ip_h=abs(hilbert(Ip));



position=[0:total_OPD/(length(Ip)-1):total_OPD]';  

%plot(1:length(Is),Is_f,1:length(Ip),Ip_f);

%% phase (and to filter

phase=atan(Ip./Is)*180/pi;

%plot(1:length(Is),Is,1:length(Ip),Ip);
%Ip_background=[Ip(1):(Ip(length(Is))-Ip(1))/(length(Ip)-1):Ip(end)]';

%Ip=Ip-Ip_background;

%FFTIs=fft(Is);

%FFTIp=fft(Ip);

%FFTIs(1:240)=0;

%FFTIs(4800-240:4800)=0;

%FFTIp(1:240)=0;

%FFTIp(4800-240:4800)=0;

%Is_filtered=abs(ifft(FFTIs));

%Ip_filtered=abs(ifft(FFTIp));

%clear Is_o Ip_o;

%phase=(Ip-Ip(1))./(Is-Is(1));


pixel_size=total_OPD/length(Ip);   %micron

spatial_PSF=gaussmf(position,[axial_resolution 0]);

spectrum=fft(spatial_PSF);

%plot(1:length(Is),phase);

%plot(position,Is,position,Ip);

plot(position,phase);

Is_3=Is*1E3;

Ip_3=Ip*1E3;

[first_max_s first_max_index_s]=max(Is(1:saparation_index));
first_max_position_s=position(first_max_index_s);

[first_max_p first_max_index_p]=max(Ip(1:saparation_index));
first_max_position_p=position(first_max_index_p);


%% to find first peak (for Is_carrier and Ip_carrier)

%space: position (micron)
%power: Is_carrier


[Is_peakvalue Is_peakindex]=max(real(Is_carrier));
Is_peakposition=position(Is_peakindex);
FWHM_Ispeak_right=position(find(Is>0.5*Is_peakvalue, 1, 'last'));
FWHM_Ispeak_left=position(find(Is>0.5*Is_peakvalue, 1, 'first'));
FWHM_Ispeak=FWHM_Ispeak_right-FWHM_Ispeak_left;

position_Is_new=position-Is_peakposition;      %for interpolation of known PSF substraction


[Ip_peakvalue Ip_peakindex]=max(real(Ip_carrier(1:Is_peakindex+500)));


Ip_carrier=circshift(Ip_carrier,Is_peakindex-Ip_peakindex);    %s-pol shift by 200 pixel to minimize total phase retardation 275 for mid min 12 for minmin, 5 for poeriod min

%% first peak generation

N_f=2^12;
N_t=N_f*4;

spectrumData=dlmread('D:\Ce.txt');

wavelength_old=spectrumData(:,1);
spectrumPower_old=spectrumData(:,2);

spectrumPower_old=spectrumPower_old-mean(spectrumPower_old(1:1000));

plot(wavelength_old,spectrumPower_old);



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
spectrumPower=interp1(freq,spectrumPower_old,fx);
spectrumPower(isnan(spectrumPower))=0;


CS=fft(spectrumPower,N_t)';     %with minus time

CS=fftshift(CS);
%CS=real(fft(S_padded));     %with minus time
CS_normal=CS/max(abs(CS));
CS_envelope=abs(CS);             

[inter_peakvalue_FD inter_peakindex_FD]=max(CS_envelope);
inter_position_FD=space_FD(inter_peakindex_FD);
FWHM_FD_right=space_FD(find(CS_envelope>0.5*inter_peakvalue_FD, 1, 'last'));
FWHM_FD_left=space_FD(find(CS_envelope>0.5*inter_peakvalue_FD, 1, 'first'));
FWHM_inter_FD=FWHM_FD_right-FWHM_FD_left;
FWHM_FD=FWHM_inter_FD;

space_FD=space_FD-inter_position_FD;

plot(space_FD,real(CS));


CS_interp=interp1(space_FD,CS,position_Is_new);


%% first peak substraction 

CS_real=real(CS_interp);

Is_carrier_sub=Is_carrier-CS_interp/max(real(CS_interp))*Is_peakvalue*0.9;  %this 0.9 can be critical...

Ip_carrier_sub=Ip_carrier-CS_interp/max(real(CS_interp))*Ip_peakvalue*0.9;


%plot(position,Is_filtered,position,Ip_filtered);

plot(position,real(Is_carrier),position,real(Ip_carrier),position,real(Is_carrier_sub),position,real(Ip_carrier_sub)); figure(gcf)
%plot(position,log(abs(FFTIs)),position,log(abs(FFTIp)));



%% PHase handling

for j=1:length(angs)
    if j<length(angs)
    if angs(j+1) > angs(j) + 0.95*2*pi
        angs(j+1:end)=angs(j+1:end)-2*pi;
    elseif angs(j+1) < angs(j) - 0.9*2*pi       
        angs(j+1:end)=angs(j+1:end)+2*pi;
    end
    
    if angp(j+1) > angp(j) + 0.95*2*pi
        angp(j+1:end)=angp(j+1:end)-2*pi;
    elseif angp(j+1) < angp(j) - 0.9*2*pi       
        angp(j+1:end)=angp(j+1:end)+2*pi;
    end
    end
end

angs=circshift(angs,-22);    %s-pol shift by 200 pixel to minimize total phase retardation 275 for mid min 12 for minmin, 5 for poeriod min
ang_diff=angs-angp;
%ang_diff(length(ang_diff):length(ang_diff))=0;
ang_diff(length(ang_diff)-30:length(ang_diff))=0;


%plot(position,angs,position,angp,position,ang_diff,position,Ip*10000000); figure(gcf)

%plot(position,real(Is_carrier),position,real(Is_carrier_sub),position,CS_real); figure(gcf)


%plot(position,Is_carrier_sub,position,Ip_carrier_sub); figure(gcf)
