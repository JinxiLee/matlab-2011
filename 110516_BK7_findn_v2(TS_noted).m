clear all; close all; clc
% sample: BK7, 6mm
temp=dlmread('Cr.txt');
lambda=temp(:,1);
lambda=lambda*0.001;
S0=temp(:,2);
d=6000;
slope=0.05;
N_fx1=4096;
N_fx2=ceil(slope*d/min(lambda)*5);
if N_fx1>=N_fx2
    N_fx=N_fx1;
end
if N_fx1<N_fx2
    N_fx=N_fx2;
end
lambda=interp1(1.1:(1.7-1.1)/1000:1.7,lambda,1.1:(1.7-1.1)/(N_fx-1):1.7);
S0=interp1(1.1:(1.7-1.1)/1000:1.7,S0,1.1:(1.7-1.1)/(N_fx-1):1.7);
dx_max=0.01;
ilambda=1./lambda;
dfx=max(ilambda)/(N_fx-1);
fx=0:dfx:max(ilambda);
N_x=2^(ceil(log2(N_fx/dx_max/max(fx))));
lambda0=lambda(S0==max(S0));

C1 = 1.03961212; C2 = 0.00600069867; C3 = 0.231792344; C4 = 0.0200179144; C5 = 1.01046945; C6 = 103.560653;
n_BK7 = sqrt(1+ C1.*lambda.^2./(lambda.^2-C2) + C3.*lambda.^2./(lambda.^2-C4) + C5*lambda.^2./(lambda.^2-C6));
n0_BK7= sqrt(1+ C1.*lambda0.^2./(lambda0.^2-C2) + C3.*lambda0.^2./(lambda0.^2-C4) + C5*lambda0.^2./(lambda0.^2-C6));

read_data='experiment'; % experiment % simulation (固定data) % simulation_1 (可調整)
process='same'; % cut % move(round) % same % noise % shift (zero) 
save='n'; % y/n
transform='fft'; % fft/ ft 
z0=-100;
for time=1:1
switch lower(read_data)
    case 'experiment'
        beam='focusedd'; %parallel/focusedd
        number=24;
        s_number=sprintf('%d',number);
        if beam=='focusedd'
            filename=(strcat('D:\xlemon\101222\Z Scan\Rename\',sprintf('%d.txt',number)));
        end
        if beam=='parallel'
            filename=(strcat('D:\xlemon\101224\Z Scan\Rename\',sprintf('%d.txt',number)));
        end            
        temp1= dlmread(filename);
        CS1=temp1(:,2);
        xx=temp1(:,1);
    case 'simulation'
       CS1=importdata('CS_BK7_115000.txt'); % N_x=115000
       xx=dlmread('xx_115000.txt');
    case 'simulation1'
%         z0=200;
        S=S0.*exp(i*4*pi./lambda*d.*(n_BK7-n0_BK7)).*exp(i*2*pi./lambda.*z0);
        S1=interp1(ilambda,S.*lambda.^2,fx);
        S1(isnan(S1))=0;
        N_x=4096*8*2;
        CS1=real(fftshift(fft(S1,N_x)));
        CS1=CS1/max(CS1);
        xx=1/2*(-1/dfx*(1/2):1/dfx/N_x:1/dfx*(1/2-1/N_x));
end

N_x1=size(CS1,1);
N_x2=size(CS1,2);
if N_x1>N_x2
    N_x=N_x1;
else N_x=N_x2;
end

switch lower(process)
    case 'cut'    %substitute with zero
        CS2=zeros(1,N_x);
        idxx=find(CS1==max(CS1));
        range=5000;
        CS2(idxx-range:idxx+range)=CS1(idxx-range:idxx+range);
        CS1=CS2;
    case 'move'
        CS2=CS1;
        ff=find(CS1==max(CS1));            % 把data shift到正確位置
        fff=find(xx>max(xx)/2+100.35,1);
%         fff=find(xx>max(xx)/2+(time*10),1);
        CS1(fff-ff+1:fff-ff+N_x)=CS2;
        CS1(1:fff-ff)=CS2(N_x-fff+ff+1:N_x);
        CS1=CS1(1:N_x);
    case 'same'
    case 'noise'
        noise=0.02;
        CS1=CS1+(rand(1,N_x)-0.5)*noise;
    case 'shift'
        CS2=zeros(N_x,1);
        ff=find(CS1==max(CS1));           % 把data shift到正確位置
        fff=find(xx>max(xx)/2+100.35,1);  % correct
%         fff=find(xx>max(xx)/2+time*10,1);
        CS2(fff-N_x/20:fff+N_x/20)=CS1(ff-N_x/20:ff+N_x/20);
        CS1=CS2;
        CS1=CS1(1:115000);
    case 'padding'
        s_padding=sprintf('%d',padding);
end
CS=CS1/max(CS1);
SNR=10*log10(1/var(CS(1:11500)));
dx=(xx(end)-xx(1))/(N_x-1);
N_fx=N_x;
df=1/2/N_x/dx;
fx=0:df:df*(N_fx-1);

switch lower(transform)
    case 'fft'
        S01=ifft(ifftshift(CS1));
    case 'ft' % 寫不出來
        syms x
        for nfx=1:N_fx;
            f=fx(nfx);
            Fun=@(x)CS1.*exp(i*f*(x));
            S01(nfx)=quad(Fun,-10,10);
        end
end
S1=interp1(fx,S01,ilambda);
S1=S1./(lambda.^2);

h_CS=abs(hilbert(CS1));
h_CS=h_CS/max(h_CS);
FWHM(time)=xx(find(h_CS>0.5,1,'last'))-xx(find(h_CS>0.5,1,'first'));
peak(time)=xx(find(h_CS==max(h_CS)));
FWHM0=3.0414;

figure(1);plot(xx,CS1/max(CS1),'k',xx,h_CS,'b--')
xlabel('Optical Path Difference (\mum)')
ylabel('Intensity (a.u.)')
legend('Carrier','Envelope')
switch lower(read_data)
    case 'simulation1'
        title(strcat('Simulation Shift ',sprintf('%d',z0/2),'\mum : BK7 6mm Dispersive Signal'))
    case 'simulation'
        title(strcat(read_data,' (',process,') ',': BK7 6mm Dispersive Signal'))
    case 'experiment'
        title(strcat(read_data,' (',process,') ',': BK7 6mm Dispersive Signal'))
end
set(gcf,'color',[1 1 1])
F_CS1=getframe(gcf);

figure(2);plot(xx,CS1/max(CS1),'k',xx,h_CS,'b--')
xlabel('Optical Path Difference (\mum)')
ylabel('Intensity (a.u.)')
legend('Carrier','Envelope')
switch lower(read_data)
    case 'simulation1'
        title(strcat('Simulation Shift ',sprintf('%d',z0/2),'\mum : BK7 6mm Dispersive Signal'))
    case 'simulation'
        title(strcat(read_data,' (',process,') ',': BK7 6mm Dispersive Signal'))
    case 'experiment'
        title(strcat(read_data,' (',process,') ',': BK7 6mm Dispersive Signal','(',s_number,')'))
end
switch lower(read_data)
    case 'simulation'
        axis([80 120 -1 1])
    case 'simulation1'
        axis([80+z0/2 120+z0/2 -1 1])
    case 'experiment'
        axis([peak(time)-20 peak(time)+20 -1 1])
end
set(gcf,'color',[1 1 1])
F_CS1_enlarge=getframe(gcf);

figure(3);plot(fx,S01/max(S01))
xlabel('frequency')
axis([0.5 1 -1 1])
set(gcf,'color',[1 1 1])
F_frequency=getframe(gcf);

figure(4);plot(lambda,S0/max(S0),'k--',lambda,S1/max(S1),'b')
legend('S0','S1')
xlabel('Wavelength(\mum)')
ylabel('Nomalized intensity (a.u.)')
set(gcf,'color',[1 1 1])
F_spectrum=getframe(gcf);

t1=log(S1./S0);
t1=imag(t1);
si1=size(t1,1);
si2=size(t1,2);
if si1>si2  % 確認找到的是正確的，不是1
    si=si1;
else si=si2; 
end
index=find(S0==max(S0));  % 切成>lambd0和<lambda0兩部份
for p=1:index                % 先找<lambda0的部份
    t2(p)=t1(index+1-p);
end
m1=zeros(1,si);
for s=1:index-1;
    if t2(s+1)<t2(s)-4
        m1(s+1:index)=m1(s+1:index)+1;
    end
end
for r=1:index
    m2(r)=m1(index+1-r);
end
for pp=index+1:si              % 算<lambda0的
    if t1(pp)>t1(pp-1)+4.2
        m1(pp:si)=m1(pp:si)-1;
    end
end
m2(index+1:si)=m1(index+1:si);
if si1>si2
    m2=m2';
end

t3=t1+2*pi.*m2;

n1=n0_BK7+(t3./4/pi/d.*lambda);
n0_1=n1(find(lambda==lambda0));
n01=n0_BK7+(t1./4/pi/d.*lambda);
zz=polyfit(lambda,n1,3);
z=polyfit(lambda,n1,2);
zzz=polyfit(lambda,n1,1);
n2=zz(1)*lambda.^3+zz(2)*lambda.^2+zz(3)*lambda.^1+zz(4);
n3=z(1)*lambda.^2+z(2)*lambda.^1+z(3);
n4=zzz(1)*lambda+z(2);
figure(5);plot(lambda,n_BK7,'k',lambda,n01,'m-.',lambda,n1,'r:',lambda,n3,'b--')
xlabel('Wavelength (\mum)')
ylabel('Refractive index')
legend('BK7','Back BK7','Reconstructed BK7','Polyfit 2')
set(gcf,'color',[1 1 1])
F_n=getframe(gcf);

figure(6);plot(lambda,n_BK7,'k',lambda,n3,'b--')
xlabel('Wavelength (\mum)')
ylabel('Refractive index')
legend('BK7','Derived BK7')
set(gcf,'color',[1 1 1])
F_n1=getframe(gcf);

% 用新建的n 來看訊號寬度是否與模擬相同
na=sin(atan(0.4));
thida=asin(na)*180/pi;
ddd=cot(asin(na));
N_fx=4096;
N_x=4096*32;
temp1=zeros(100,N_x);
ilambda=1./lambda;
dfx=max(ilambda)/(N_fx-1);
fx=0:dfx:max(ilambda);
% focused
if beam=='focusedd'
xx=1:N_x;
g=gaussmf(xx, [15500 ,length(xx)/2]);
gg=g(N_x/2:N_x/399:end);
ggd=0:0.005:199*0.005;
% figure(9);plot(ggd,gg)
% xlabel('off center (\mum)')
% ylabel('factor')
for pp=1:200
dd=0.005*(pp-1);
d2=6000./cos(atan(dd/ddd));
S21=S0.*exp(i*4*pi./lambda*d2.*(n1-n0_1));
ilambda=1./lambda;
dfx=max(ilambda)/(N_fx-1);
fx=0:dfx:max(ilambda);
T1=interp1(ilambda,S21.*lambda.^2,fx);
T1(isnan(T1))=0;
CS21=real(fftshift(fft(T1,N_x)));
x=1/2*(-1/dfx*(1/2):1/dfx/N_x:1/dfx*(1/2-1/N_x));
CS_t=CS21;
CSS_t=CS_t;
CS_t=CS_t/max(CS_t);
CS_t_envelope=abs(hilbert(CS_t));
% temp1(pp,:)=CS_t_envelope*(sinc((pp-1)*0.005)); %sinc(0:0.005:1).^2
temp1(pp,:)=CS_t_envelope*gg(pp);
ff(pp)=x(find(CS_t_envelope>max(CS_t_envelope)/2,1,'last'))-x(find(CS_t_envelope>max(CS_t_envelope)/2,1,'first'));
end
CS1_t_envelope=mean(temp1)/max(mean(temp1));
FWHM21(time)=x(find(CS1_t_envelope>max(CS1_t_envelope)/2,1,'last'))-x(find(CS1_t_envelope>max(CS1_t_envelope)/2,1,'first'));
end
% parallel
if beam=='parallel'
S=S0.*exp(i*4*pi./lambda*d.*(n1-n0_1));
S2=interp1(ilambda,S.*lambda.^2,fx);
S2(isnan(S2))=0;
CS2=real(fftshift(fft(S2,N_x)));
CS2=CS2/max(abs(CS2));
h_CS2=abs(hilbert(CS2));
h_CS2=h_CS2/max(h_CS2);
x=1/2*(-1/dfx*(1/2):1/dfx/N_x:1/dfx*(1/2-1/N_x));
FWHM2(time)=x(find(h_CS2>0.5,1,'last'))-x(find(h_CS2>0.5,1,'first'));
peak2(time)=x(find(h_CS2==max(h_CS2)));
% diff_FWHM(time)=FWHM2(time)-FWHM0; 
% diff_peak(time)=peak2(time)-peak(time);
% if diff_peak(time)<1
%     peak_position=peak;
% end    
figure(7);plot(x,CS2,'k',x,h_CS2,'b--')
title('Forward Signal')
legend('Carrier','Envelope')
xlabel('Optical Path Difference (\mum)')
ylabel('Intensity (a.u.)')
axis([peak2(time)-20 peak2(time)+20 -1 1])
set(gcf,'color',[1 1 1])
F_CS2=getframe(gcf);
end
% 使拓寬的變窄
S3=S1.*exp(-i*4*pi./lambda.*(n1-n0_1)*d);
S3=interp1(ilambda,S3.*lambda.^2,fx);
S3(isnan(S3))=0;
CS3=real(fftshift(fft(S3,N_x)));
CS3=CS3/max(abs(CS3));
h_CS3=abs(hilbert(CS3));
h_CS3=h_CS3/max(h_CS3);
x=1/2*(-1/dfx*(1/2):1/dfx/N_x:1/dfx*(1/2-1/N_x));
FWHM3(time)=x(find(h_CS3>0.5,1,'last'))-x(find(h_CS3>0.5,1,'first'));
peak3(time)=x(find(h_CS3==max(h_CS3)));
if peak3(time)<0.01
    peak_position=peak(time);
end
figure(8);plot(x,CS3,'k',x,h_CS3,'b--')
title('Compensated Signal')
legend('Carrier','Envelope')
xlabel('Optical Path Difference (\mum)')
ylabel('Intensity (a.u.)')
axis([-7 7 -1 1])
set(gcf,'color',[1 1 1])
F_CS3=getframe(gcf);

switch lower(save)
    case 'y'           
        switch lower(process)
            case 'cut'
                if read_data=='experiment'
                    imwrite(F_CS1.cdata,['D:\xlemon\result\BK7\',read_data,'\cut_',s_number,'.tif']);
                    imwrite(F_CS1_enlarge.cdata,['D:\xlemon\result\BK7\',read_data,'\cut_1_',s_number,'.tif']);
                    imwrite(F_frequency.cdata,['D:\xlemon\result\BK7\',read_data,'\frequency_cut_',s_number,'.tif']);
                    imwrite(F_spectrum.cdata,['D:\xlemon\result\BK7\',read_data,'\spectrum_cut_',s_number,'.tif']);
                    imwrite(F_n.cdata,['D:\xlemon\result\BK7\',read_data,'\n_cut_',s_number,'.tif']);
                    imwrite(F_n1.cdata,['D:\xlemon\result\BK7\',read_data,'\n1_cut_',s_number,'.tif']);
                    imwrite(F_CS2.cdata,['D:\xlemon\result\BK7\',read_data,'\CS2_cut_',s_number,'.tif']);
                    imwrite(F_CS3.cdata,['D:\xlemon\result\BK7\',read_data,'\CS3_cut_',s_number,'.tif']);
                end
                imwrite(F_CS1.cdata,['D:\xlemon\result\BK7\',read_data,'\cut.tif']);
                imwrite(F_CS1_enlarge.cdata,['D:\xlemon\result\BK7\',read_data,'\cut_1.tif']);
                imwrite(F_frequency.cdata,['D:\xlemon\result\BK7\',read_data,'\frequency_cut.tif']);
                imwrite(F_spectrum.cdata,['D:\xlemon\result\BK7\',read_data,'\spectrum_cut.tif']);
                imwrite(F_n.cdata,['D:\xlemon\result\BK7\',read_data,'\n_cut.tif']);
            case 'move'
                imwrite(F_CS1.cdata,['D:\xlemon\result\BK7\',read_data,'\move.tif']);
                imwrite(F_CS1_enlarge.cdata,['D:\xlemon\result\BK7\',read_data,'\move_1.tif']);
                imwrite(F_frequency.cdata,['D:\xlemon\result\BK7\',read_data,'\frequency_move.tif']);
                imwrite(F_spectrum.cdata,['D:\xlemon\result\BK7\',read_data,'\spectrum_move.tif']);
                imwrite(F_n.cdata,['D:\xlemon\result\BK7',read_data,'\n_move.tif']);
            case 'noise'
                imwrite(F_CS1.cdata,['D:\xlemon\result\BK7\',read_data,'\noise_',sprintf('%.1f',noise),'.tif']);
                imwrite(F_CS1_enlarge.cdata,['D:\xlemon\result\BK7\',read_data,'\noise_',sprintf('%.1f',noise),'_1.tif']);
                imwrite(F_frequency.cdata,['D:\xlemon\result\BK7\',read_data,'\frequency_noise_',sprintf('%.1f',noise),'.tif']);
                imwrite(F_spectrum.cdata,['D:\xlemon\result\BK7\',read_data,'\spectrum_noise_',sprintf('%.1f',noise),'.tif']);
                imwrite(F_n.cdata,['D:\xlemon\result\BK7',read_data,'\n_noise_',sprintf('%.1f',noise),'.tif']);
            case 'shift'
                imwrite(F_CS1.cdata,['D:\xlemon\result\BK7\',read_data,'\shift_',s_number,'.tif']);
                imwrite(F_CS1_enlarge.cdata,['D:\xlemon\result\BK7\',read_data,'\shift_1_',s_number,'.tif']);
                imwrite(F_frequency.cdata,['D:\xlemon\result\BK7\',read_data,'\frequency_shift_',s_number,'.tif']);
                imwrite(F_spectrum.cdata,['D:\xlemon\result\BK7\',read_data,'\spectrum_shift_',s_number,'.tif']);
                imwrite(F_n.cdata,['D:\xlemon\result\BK7\',read_data,'\n_shift_',s_number,'.tif']);
                imwrite(F_n1.cdata,['D:\xlemon\result\BK7\',read_data,'\n1_shift_',s_number,'.tif']);
                imwrite(F_CS2.cdata,['D:\xlemon\result\BK7\',read_data,'\CS2_shift_',s_number,'.tif']);
                imwrite(F_CS3.cdata,['D:\xlemon\result\BK7\',read_data,'\CS3_shift_',s_number,'.tif']);
            case 'same' % simulaton1 + same = shift
                switch lower(read_data)
                    case 'simulation1'
                        imwrite(F_CS1.cdata,['D:\xlemon\result\BK7\',read_data,'\shift_',sprintf('%d',z0/2),'.tif']);
                        imwrite(F_CS1_enlarge.cdata,['D:\xlemon\result\BK7\',read_data,'\shift_1_',sprintf('%d',z0/2),'.tif']);
                        imwrite(F_frequency.cdata,['D:\xlemon\result\BK7\',read_data,'\frequency_shift_',sprintf('%d',z0/2),'.tif']);
                        imwrite(F_spectrum.cdata,['D:\xlemon\result\BK7\',read_data,'\spectrum_shift_',sprintf('%d',z0/2),'.tif']);
                        imwrite(F_n.cdata,['D:\xlemon\result\BK7\',read_data,'\n_shift_',sprintf('%d',z0/2),'.tif']);
                        imwrite(F_n1.cdata,['D:\xlemon\result\BK7\',read_data,'\n1_shift_',sprintf('%d',z0/2),'.tif']);
                        imwrite(F_CS2.cdata,['D:\xlemon\result\BK7\',read_data,'\CS2_shift_',sprintf('%d',z0/2),'.tif']);
                        imwrite(F_CS3.cdata,['D:\xlemon\result\BK7\',read_data,'\CS3_shift_',sprintf('%d',z0/2),'.tif']);
                    case 'experiment'
                        imwrite(F_CS1.cdata,['D:\xlemon\result\BK7\',read_data,'\_',s_number,'.tif']);
                        imwrite(F_CS1_enlarge.cdata,['D:\xlemon\result\BK7\',read_data,'\_1_',s_number,'.tif']);
                        imwrite(F_frequency.cdata,['D:\xlemon\result\BK7\',read_data,'\frequency_',s_number,'.tif']);
                        imwrite(F_spectrum.cdata,['D:\xlemon\result\BK7\',read_data,'\spectrum_',s_number,'.tif']);
                        imwrite(F_n.cdata,['D:\xlemon\result\BK7\',read_data,'\n_',s_number,'.tif']);
                        imwrite(F_n1.cdata,['D:\xlemon\result\BK7\',read_data,'\n1_',s_number,'.tif']);
                        imwrite(F_CS2.cdata,['D:\xlemon\result\BK7\',read_data,'\CS2_',s_number,'.tif']);
                        imwrite(F_CS3.cdata,['D:\xlemon\result\BK7\',read_data,'\CS3_',s_number,'.tif']);
%                         dlmwrite('D:\xlemon\result\BK7\',read_data,'\FWHM3_',s_number,'.txt',FWHM3);
%                         dlmwrite('D:\xlemon\result\BK7\',read_data,'\FWHM2_',s_number,'.txt',FWHM2);
%                         dlmwrite('D:\xlemon\result\BK7\',read_data,'\FWHM_',s_number,'.txt',FWHM);
%                         dlmwrite('D:\xlemon\result\BK7\',read_data,'\peak3_',s_number,'.txt',peak3);
%                         dlmwrite('D:\xlemon\result\BK7\',read_data,'\peak2_',s_number,'.txt',peak);
%                         dlmwrite('D:\xlemon\result\BK7\',read_data,'\peak_',s_number,'.txt',peak);
                end
            case 'padding'
                imwrite(F_CS1.cdata,['D:\xlemon\result\BK7\',read_data,'\CS_',s_padding,'.tif']);
                imwrite(F_CS1_enlarge.cdata,['D:\xlemon\result\BK7\',read_data,'\CS_1_',s_padding,'.tif']);
                imwrite(F_frequency.cdata,['D:\xlemon\result\BK7\',read_data,'\frequency_',s_padding,'.tif']);
                imwrite(F_spectrum.cdata,['D:\xlemon\result\BK7\',read_data,'\spectrum_',s_padding,'.tif']);
                imwrite(F_n.cdata,['D:\xlemon\result\BK7\',read_data,'\n_',s_padding,'.tif']); 
                imwrite(F_n1.cdata,['D:\xlemon\result\BK7\',read_data,'\n1_',s_padding,'.tif']);
                imwrite(F_CS2.cdata,['D:\xlemon\result\BK7\',read_data,'\CS2_',s_padding,'.tif']);
                imwrite(F_CS3.cdata,['D:\xlemon\result\BK7\',read_data,'\CS3_',s_padding,'.tif']);
                FF=[FWHM FWHM2 FWHM3];
                filename=strcat('D:\xlemon\result\BK7\',read_data,'\padding_FWHM_',sprintf('%d.txt',padding));
                dlmwrite(filename,FF);
        end
    case 'n',
end
z0=z0+20;
end

