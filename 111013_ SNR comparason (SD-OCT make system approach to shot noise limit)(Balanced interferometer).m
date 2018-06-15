clear all;

averaging=700;
wavelength=0.56;    %(micron, 暫時用center wavelength, 實際上應該是要用spectrum)
linerate=70000;
N=4096;                             %number of pixel
Plog=[-6:0.01:0]';
P=10.^(Plog);              %A: 3.7mW, B: 4.6mW, C and D: 16mW
T=1/linerate;                % sec, for a A-scan
BW=0.1;             %(micron)
Ps=0.25*P;
Pr=0.25*P;

qe=0.6;                             %(CCD quantum efficiency (electron to photon))
      %for a A-scan


h=6.626E-34;        %plank constant
c=3E8;              %light speed
d_freq=(c/((wavelength-0.5*BW)*1E-6))-(c/((wavelength+0.5*BW)*1E-6));

%% for each pixel 

Q_DC=qe*(Pr+Ps)/(h*(c/(wavelength*10^-6)))/N*T;
Q_sig=2*qe*sqrt(Pr.*Ps)/(h*(c/(wavelength*10^-6)))/N*T*averaging;
dark_variance=(14^2)*averaging;
shot_variance=Q_DC*averaging;
excess_variance=(Q_DC.^2)/T/d_freq*averaging;

%% after sumup and fft

% sumup

dark_variance_total=dark_variance*N;
shot_variance_total=shot_variance*N;
excess_variance_total=excess_variance*N;
Q_sig_total=Q_sig*N/2;                  %/2 for fft

SNR=(Q_sig_total.^2)./(dark_variance_total+shot_variance_total+excess_variance_total);

SNR_dark=(Q_sig_total.^2)./dark_variance_total;
SNR_shot=(Q_sig_total.^2)./shot_variance_total;
SNR_excess=(Q_sig_total.^2)./excess_variance_total;

SNRlog_dark=10*log10(SNR_dark);
SNRlog_shot=10*log10(SNR_shot);
SNRlog_excess=10*log10(SNR_excess);


SNRlog=10*log10(SNR);

plot(Plog,SNRlog,Plog,SNRlog_dark,Plog,SNRlog_shot,Plog,SNRlog_excess);



