clear all;
center_lambda=0.760;  %micron
bandwidth=0.180;      %micron
FWHM_pixel=600;     %(pixel)
lambda_res=bandwidth/FWHM_pixel;
lens_lambda_res=0.5*lambda_res;
x=[0:0.01:1500]';       %micron
R=10*log10((sinc(2*(lambda_res/center_lambda^2).*x).^2).*exp(-4/log(2)*((lens_lambda_res/center_lambda^2).*x).^2));
index_3dB=find(R<-3,1,'first');
index_10dB=find(R<-10,1,'first');
index_0=find(x>(2*(lambda_res/center_lambda^2))^-1,1,'first');
x_0=x(index_0);
x_3dB=x(index_3dB);
x_10dB=x(index_10dB);
x=x/1000;
plot(x,R);