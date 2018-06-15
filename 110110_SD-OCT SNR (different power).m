clear all;

%P=0.001;            %(log(Watt), incident power)
Plog=[-8:0.1:-3]';  %dB
wavelength=0.76;    %(micron, 暫時用center wavelength, 實際上應該是要用spectrum)
%Rslog=[-10:0.1:0];    %(log10(reflectance )of sample arm)

%Rs=10.^Rslog;
P=(10.^Plog);  %W
%Ps=0.5*Rs*P;
%Pr=0.5*P;           %(power back from ref arm, assumed rr=1)
Ps=0.5*P;
Pr=Ps;

qe=0.5;                             %(CCD quantum efficiency (electron to photon))
N=4096;                             %number of pixel
%framrate=1;                        %(Hz)
%Nlateralpixel=1000;                  %(for B-scan)
%duty=1;                  
ti=1/1000;        %(second, integration time corresponding to the B-scan rate)

h=6.626E-34;        %plank constant
c=3E8;              %light speed
d_freq=(c/(0.65E-6))-(c/(0.85E-6));

Sig=(2*qe*((Pr.*Ps).^0.5)/h/(c/(wavelength*10^-6)))*ti/N*N^0.5;        %for one pixel after FFT, N is timed back asuumed FFT
shut_2=(qe*(Pr+Ps)/(h*(c/(wavelength*10^-6))))*ti/N;                   %noise is not timed back after FFT
dark_2=(1.5/4096*117500)^2;                                               %full well:117500, dark noise=1.5LSB
excess_2=((qe/(h*(c/(wavelength*10^-6))))^2)*ti*((Pr+Ps).^2)/N/d_freq;
SNRshut=10*log10((Sig.^2)./shut_2);
SNRexcess=10*log10((Sig.^2)./excess_2);
SNRdark=10*log10((Sig.^2)./dark_2);
SNRtotal=10*log10((Sig.^2)./(shut_2+excess_2+dark_2));
plot(Plog,SNRshut,Plog,SNRexcess,Plog,SNRdark,Plog,SNRtotal);
ti_ms=ti*1000;



