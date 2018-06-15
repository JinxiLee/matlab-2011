clear all;
close all;
clc;

%% PSF generation

N_f=2^12;
N_t=N_f*4;

spectrumData=dlmread('D:\Ce.txt');

wavelength_old=spectrumData(:,1);
spectrumPower_old=spectrumData(:,2);

spectrumPower_old=spectrumPower_old-mean(spectrumPower_old(1:1000));

c=3E8;                     %m/sec
freq=c./(wavelength_old*1E-9);     %Hz orignal freq array
d_f=max(freq)/(N_f-1);
fx=0:d_f:max(freq);        %freq after interpolation
d_t=1/(d_f*(N_t-1));
%d_t=1/(d_f*2*N_f);
time=[-0.5*(N_t-1)*d_t:d_t:0.5*N_t*d_t]'/2;%/2是因為一來一回
%time=[-(N_f-1)*d_t:d_t:N_f*d_t]/2;  %/2是因為一來一回
space_FD=c*time; % 暫定只用array前1/100 (大約也有100 micron左右吧)
space_FD=(space_FD-space_FD(1))*1E6;                     % shift to zero & m to micron
%plot(lambda,inter,lambda2,inter);

%S0=inter;
%plot(lambda,S0/max(S0),lambda,inter/max(inter),lambda,ref/max(ref));
%S0=S0/max(S0);        %不先normal應該也沒差?
spectrumPower=interp1(freq,spectrumPower_old.*wavelength_old.^2,fx);
spectrumPower(isnan(spectrumPower))=0;


CS=fft(spectrumPower,N_t)';     %with minus time

CS=fftshift(CS);
%CS=real(fft(S_padded));     %with minus time
CS_normal=CS/max(abs(CS));
CS_envelope=abs(CS);             

[inter_peakvalue_FD inter_peakindex_FD]=max(CS_envelope);
inter_position_FD=space_FD(inter_peakindex_FD);
%FWHM_FD_right=space_FD(find(CS_envelope>0.5*inter_peakvalue_FD, 1, 'last'));
%FWHM_FD_left=space_FD(find(CS_envelope>0.5*inter_peakvalue_FD, 1, 'first'));
%FWHM_inter_FD=FWHM_FD_right-FWHM_FD_left;
%FWHM_FD=FWHM_inter_FD;

space_FD=space_FD-inter_position_FD;

%PSF=real(CS_normal((round(length(space_FD)/2-250)):(round(length(space_FD)/2+250))));
PSF=real(CS_normal)';


%% Sine generation

period=28;   %micron

cosinefunction(1:length(space_FD))=0;
sinefunction(1:length(space_FD))=0;
for j=1:length(space_FD)
if (space_FD(j)<period/8) && (space_FD(j)>-period/8)
cosinefunction(j)=100*cos(2*pi*(space_FD(j)+period/8)/period);
sinefunction(j)=100*sin(2*pi*(space_FD(j)+period/8)/period);
end
end



%% Manual convolution

% only partial of PSF

[PSF_max_value PSF_max_index]=max(PSF);

partial_index=4096;

PSF_o=PSF;
PSF(1:(PSF_max_index-partial_index-1))=0;
PSF((PSF_max_index+partial_index+1):end)=0;


cosinefunction_convoluted(1:length(cosinefunction))=0;
sinefunction_convoluted(1:length(sinefunction))=0;

for j=1:length(cosinefunction)
cosinefunction_temp=circshift(PSF,[0 -(PSF_max_index-j)]).*cosinefunction;
sinefunction_temp=circshift(PSF,[0 -(PSF_max_index-j)]).*sinefunction;
cosinefunction_convoluted(j)=sum(cosinefunction_temp);
sinefunction_convoluted(j)=sum(sinefunction_temp);
end

%% To integrate PSF

% to find the center of PSF




%test(1:length(space_FD))=0;

%test((round(length(space_FD)/2-10)):(round(length(space_FD)/2+10)))=1;





%PSF(1:100)=0;
%PSF(50)=1;

%clear time CS_envelope freq CS_normal CS_envelope

%result=convn(sinefunction,PSF,'same');
%result_C=conv(cosinefunction,PSF);
%result_S=conv(sinefunction,PSF);
%sub(1:length(result))=0;

%%for j=1:length(sub)
 %   if (j-(7881-252))>0 && (j-(7881-252))<length(PSF)
  %  sub(j)=sub(j)+PSF(j-(7881-252));
   % end
    %if (j-(9007-252))>0 && (j-(9007-252))<length(PSF)
    %sub(j)=sub(j)+PSF(j-(9007-252));
    %end
%end

%sub=sub/max(sub)*max(result);

%result=result-sub;
%result=convn(test,sinefunction,'same');
%plot(space_FD,test,space_FD,sinefunction,space_FD,result);


%result_C=result_C((1+length(PSF)/2):(length(result_C)-length(PSF)/2)+1)-result_C(1);
%result_C=result_C';
%result_S=result_S((1+length(PSF)/2):(length(result_S)-length(PSF)/2)+1)-result_S(1);
%result_S=result_S';

angle=atan(sinefunction_convoluted./cosinefunction_convoluted);
%%
birefringence=diff(angle)*0.56/4/pi./diff(space_FD)';

plot(space_FD,sinefunction,space_FD,cosinefunction,space_FD,sinefunction_convoluted,space_FD,cosinefunction_convoluted);

dlmwrite('D:\sinefunction.txt',sinefunction,'delimiter','\t','newline','pc');
dlmwrite('D:\cosinefunction.txt',cosinefunction,'delimiter','\t','newline','pc');
dlmwrite('D:\space_FD.txt',space_FD,'delimiter','\t','newline','pc');
dlmwrite('D:\sinefunction_convoluted.txt',sinefunction_convoluted,'delimiter','\t','newline','pc');
dlmwrite('D:\cosinefunction_convoluted.txt',cosinefunction_convoluted,'delimiter','\t','newline','pc');
dlmwrite('D:\birefringence.txt',birefringence,'delimiter','\t','newline','pc');
dlmwrite('D:\angle.txt',angle,'delimiter','\t','newline','pc');
%plot(space_FD,result_C,space_FD,result_S);
%plot(space_FD,sinefunction,space_FD,result);
%plot(space_FD,real(CS_normal),space_FD,sinefunction);
