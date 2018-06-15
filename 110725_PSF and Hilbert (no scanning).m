clear all;
close all;
clc;

N_f=2^12;
N_t=N_f*8;

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
space_FD=c*time; % 暫定只用array前1/100 (大約也有100 micron左右吧)
space_FD=(space_FD-space_FD(1))*1E6;                     % shift to zero & m to micron

spectrumPower=interp1(freq,spectrumPower_old.*wavelength_old.^2,fx);
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

CS_real=real(CS);


disatance=0.85;    %micron

Sig1=circshift(CS_real,round(disatance*1E-6/d_t/c));
Sig2=circshift(CS_real,-round(disatance*1E-6/d_t/c));
Sig_hilbert_1=abs(hilbert(Sig1))+abs(hilbert(Sig2));
Sig=Sig1+Sig2;
Sig_hilbert_2=abs(hilbert(Sig));
[inter_peakvalue_left_1 inter_peakindex_left_1]=max(Sig_hilbert_1(1:inter_peakindex_FD));
[inter_peakvalue_right_1 inter_peakindex_right_1]=max(Sig_hilbert_1((inter_peakindex_FD+1):end));
[inter_peakvalue_left_2 inter_peakindex_left_2]=max(Sig_hilbert_2(1:inter_peakindex_FD));
[inter_peakvalue_right_2 inter_peakindex_right_2]=max(Sig_hilbert_2((inter_peakindex_FD+1):end));
distance_calc_1=(inter_peakindex_right_1+inter_peakindex_FD-inter_peakindex_left_1-1)*d_t*c/(1E-6)/2;
distance_calc_2=(inter_peakindex_right_2+inter_peakindex_FD-inter_peakindex_left_2-1)*d_t*c/(1E-6)/2;



plot(space_FD,Sig,space_FD,Sig_hilbert_2);

disatance=disatance';
distance_calc_2=distance_calc_2';
distance_calc_1=distance_calc_1';

plot(space_FD,Sig_hilbert_2,space_FD,Sig);
