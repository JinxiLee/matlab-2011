clear all;

ne0=1.607;
no=1.494;

theta=0:pi/100:pi/2;

ne=ne0*no./((ne0^2)*(sin(theta).^2)+(no^2)*(cos(theta).^2)).^0.5;

plot(theta/pi*180,ne);