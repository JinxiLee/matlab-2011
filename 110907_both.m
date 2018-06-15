clear all;



n_o=1.50:0.001:2;
k_o=0:0.0001:0.05;
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


%% T
t1=2*n1./(n+n1+i*k);
t1_1=2*(n+i*k)./(n+n1+i*k);           %_1 means reverse direction
t2=2*(n+i*k)./(n+n2+i*k);
t2_1=2*n2./(n+n2+i*k);
r1=(n1-(n+i*k))./(n+n1+i*k);
r1_1=((n+i*k)-n1)./(n+n1+i*k);
r2=((n+i*k)-n2)./(n+n2+i*k);
r2_1=(n2-(n+i*k))./(n+n2+i*k);


d=exp(i*2*pi./lambda.*(n+i*k).*l);

T_field=t1.*t2.*d./(1-r1_1.*r2.*(d.^2));

T=n2/n1*abs(T_field).^2;

Tmax=max(abs(T));

QQ=n1/n2*t1.*t2;

QQmax=max(abs(QQ));

%R=abs(r1+t1.*t1_1.*d.*r2./(1-(d.^2).*r1_1.*r2)).^2;        也不能這樣做,因為膜層有吸收......................
%T=1-R;


%% A



%R_upper=( ( ( ((n-n1).^2)  +  k.^2       )  ./ ( ((n+n1).^2)  +  k.^2       )  ).^0.5 ) .* exp(i*atan(2*k*n1./((n.^2)+(k.^2)-n1^2)));

R_upper=r1;

%R_lower=( ( ( ((n-n2).^2)  +  k.^2       )  ./ ( ((n+n2).^2)  +  k.^2       )  ).^0.5 ) .* exp(i*atan(-2*k*n2*1./((n.^2)+(k.^2)-n2^2)));

R_lower=r2;


E_upper=r1;


%E_upper_ver2=R_upper_ver2;

E_lower=(d.^2).*r2;

%E_lower_ver2=exp(2*2*pi*i/lambda.*(n+i*k).*l).*R_lower_ver2;

A=E_lower./E_upper;


A=(d.^2).*r2./r1;
%A_ver2=E_lower_ver2./E_upper_ver2;








imagesc(T,'xdata',n_o,'ydata',k_o);
imagesc(abs(A),'xdata',n_o,'ydata',k_o);
imagesc(angle(A),'xdata',n_o,'ydata',k_o);
%plot(lambda,T);