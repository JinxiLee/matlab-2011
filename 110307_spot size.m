clear all;

x=1:1/100:20;
center=1.6;   %micron
HM1=10;     %micron
HM2=9.8;     %micron

spot=2*gaussmf(x,[center/2/(2*log(2))^0.5 10])+gaussmf(x,[HM1/2/(2*log(2))^0.5 10])+gaussmf(x,[HM2/2/(2*log(2))^0.5 10]);
spot=spot/max(spot);
FWHM=x(find(spot>0.5,1,'last'))-x(find(spot>0.5,1,'first'));