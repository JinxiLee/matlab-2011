clear all;



n_o=1.70:0.001:1.80;
k_o=0:0.1:2;
k_o=k_o';

n_empty=n_o;
n_empty(:)=1;
k_empty=k_o;
k_empty(:)=1;
n=k_empty*n_o;
k=k_o*n_empty;




n1=1;
n2=1.5;
%n=1.75;
%k=0;
l=1;           %micron

lambda=0.76;       %micron


t1=2*(n+i*k)./(n+n1+i*k);           %_1 means reverse direction
t1_1=2*n1./(n+n1+i*k);
t2=2*n2./(n+n2+i*k);
t2_1=2*(n+i*k)./(n+n2+i*k);
r1=(n1-(n+i*k))./(n+n1+i*k);
r1_1=((n+i*k)-n1)./(n+n1+i*k);
r2=((n+i*k)-n2)./(n+n2+i*k);
r2_1=(n2-(n+i*k))./(n+n2+i*k);


d=exp(i*2*pi./lambda.*(n+i*k).*l);

T_field=t1.*t2.*d./(1-r1_1.*r2.*(d.^2));

T=n1/n2*abs(T_field).^2;

Tmax=max(abs(T));

QQ=n1/n2*t1.*t2;

QQmax=max(abs(QQ));

%R=abs(r1+t1.*t1_1.*d.*r2./(1-(d.^2).*r1_1.*r2)).^2;        也不能這樣做,因為膜層有吸收......................
%T=1-R;

imagesc(T,'xdata',n_o,'ydata',k_o);
%plot(lambda,T);