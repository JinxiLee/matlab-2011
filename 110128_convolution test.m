clear all;

x=0:pi/1000:20*pi;
Q1=sin(x);
Q2(1:length(x))=0;
Q2(x<10*pi)=1/1000;
Q2(x<9*pi)=0;
Q3=convn(Q1,Q2);
plot(Q3);
