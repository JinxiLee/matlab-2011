clear all


N_f=2^12;
N_t=N_f*32;

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

CS_carrier=real(CS_normal);

[inter_peakvalue_FD inter_peakindex_FD]=max(CS_envelope);
inter_position_FD=space_FD(inter_peakindex_FD);
FWHM_FD_right=space_FD(find(CS_envelope>0.5*inter_peakvalue_FD, 1, 'last'));
FWHM_FD_left=space_FD(find(CS_envelope>0.5*inter_peakvalue_FD, 1, 'first'));
FWHM_inter_FD=FWHM_FD_right-FWHM_FD_left;
FWHM_FD=FWHM_inter_FD;

space_FD=space_FD-inter_position_FD;

tag=0;
x=1;
y=1;

for j=1:length(space_FD)
    if tag ==1
        if CS_carrier(j)<0
            Zero_CS_carrier(x)=j;
            x=x+1;
            tag=-1;
        end
    elseif tag == -1
        if CS_carrier(j)>0
            Zero_CS_carrier(x)=j;
            x=x+1;
            tag=1;
        end
    end
        
     
    if tag ==0
        if CS_carrier(j)>0
            tag=1;
        else
            tag=-1;
        end
    end
    
end

Zero_CS_carrier_diff=diff(Zero_CS_carrier);

index_diff=mean(Zero_CS_carrier_diff((round(length(Zero_CS_carrier_diff)/2):(round(length(Zero_CS_carrier_diff)/2)))));
position_diff=index_diff*(space_FD(2)-space_FD(1));
