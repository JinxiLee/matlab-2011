clear all;

%P=0.001;            %(log(Watt), incident power)
Plog=[-10:0.1:0]';
wavelength=0.76;    %(micron, 暫時用center wavelength, 實際上應該是要用spectrum)
%Rslog=[-10:0.1:0];    %(log10(reflectance )of sample arm)

%Rs=10.^Rslog;
P=10.^Plog;
%Ps=0.5*Rs*P;
%Pr=0.5*P;           %(power back from ref arm, assumed rr=1)
Ps=0.5*P;
Pr=Ps;

qe=0.5;                             %(CCD quantum efficiency (electron to photon))
N=4096;                             %number of pixel

ti=1;                               %for a A-scan

h=6.626E-34;        %plank constant
c=3E8;              %light speed
d_freq=(c/(0.65E-6))-(c/(0.85E-6));

SigSD=(2*qe*((Pr.*Ps).^0.5)/h/(c/(wavelength*10^-6)))*ti/2;        %for one pixel after FFT, N is timed back asuumed FFT
shut_SD2=(qe*(Pr+Ps)/N/(h*(c/(wavelength*10^-6))))*ti*N/N;                   %noise is not timed back after FFT
dark_SD2=((1.5/4096*117500)^2)*N/N;                                               %full well:117500, dark noise=1.5LSB
excess_SD2=((qe/(h*(c/(wavelength*10^-6))))^2)*ti*((Pr+Ps).^2)/(N^2)/d_freq*N/N;

SigTD=(2*qe*((Pr.*Ps).^0.5)/h/(c/(wavelength*10^-6)))*ti/N;        %for one pixel after FFT, N is timed back asuumed FFT
shut_TD2=(qe*(Pr+Ps)/(h*(c/(wavelength*10^-6))))*ti/N;                   %noise is not timed back after FFT
dark_TD2=((1.5/4096*117500)^2);                                               %full well:117500, dark noise=1.5LSB
excess_TD2=((qe/(h*(c/(wavelength*10^-6))))^2)*ti/N*((Pr+Ps).^2)/d_freq;

%SNRshut=10*log10((Sig.^2)./shut_2);
%SNRexcess=10*log10((Sig.^2)./excess_2);
%SNRdark=10*log10((Sig.^2)./dark_2);

SNRshut_SD=10*log10((SigSD.^2)./(shut_SD2));
SNRshut_TD=10*log10((SigTD.^2)./(shut_TD2));
SNRtotal_SD=10*log10((SigSD.^2)./(shut_SD2+excess_SD2+dark_SD2));
SNRtotal_TD=10*log10((SigTD.^2)./(shut_TD2+excess_TD2+dark_TD2));

plot(Plog,SNRtotal_SD,Plog,SNRtotal_TD);



