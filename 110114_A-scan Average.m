clear all;
close all;
clc;

N_f=4096*4;
N_t=N_f;

%% To calculate SNR with envelope (assume incoherent summation between
%% frames)


data_array_all=importdata('D:\110117_ascan averge1.txt');

clear averaging_factor  data_array

averaging_factor=100;

data_array(1:fix(size(data_array_all,1)/averaging_factor),size(data_array_all,2))=0;

for j=1:fix(size(data_array_all,1)/averaging_factor)
data_array(j,:)=sum(data_array_all(1+(j-1)*averaging_factor:j*averaging_factor,:),1);
end

clear data_array_all

if size(data_array,1)>100
data_array=data_array(1:100,:);
end
%% FD data calculation

%ref=data_array(real_max_index,:);

%% Calibrate by asign 444 and 888 pixel
grating_pitch=1000000/600;      %nm
%center_wavelength=760;       %nm

wavelength_1=444;
wavelength_2=888;

index_1=3500;
index_2=1850;

spe_pixel=(wavelength_2-wavelength_1)/(index_2-index_1);   %spectral res of CCD (nm), can be minus, in fact its minus in current system
pixel=[1:length(data_array)]';  %1~4096
lambda_o=(pixel-index_1)*spe_pixel+wavelength_1;  %only used in peak finding

[dc_max_value dc_max_index]=max(data_array(:,1));
found_peak_wavelength=lambda_o(dc_max_index);
incidence_angle=asin(found_peak_wavelength/2/grating_pitch);     %假設是照此角度入射,rad

theta_1=asin((wavelength_1-0.5*found_peak_wavelength)/grating_pitch);
theta_2=asin((wavelength_2-0.5*found_peak_wavelength)/grating_pitch);

Q=(tan(theta_1-incidence_angle)-tan(theta_2-incidence_angle))/(index_1-index_2);      %roughly flens/pixel size
%%tan!
lambda=grating_pitch*sin((atan((pixel-dc_max_index)*Q)+incidence_angle))+found_peak_wavelength/2;   %atan!

%assumed_BW=137.28;

%short_coef=tan(asin((center_wavelength/2-assumed_BW/2)/grating_pitch)-incidence_angle);         %tan!
%long_coef=tan(asin((center_wavelength/2+assumed_BW/2)/grating_pitch)-incidence_angle);          %tan!



%[peak index_peak]=max(ref);
%index_long=find(ref>0.5*max(ref),1,'last');
%index_short=find(ref>0.5*max(ref),1,'first');

%Q=178/(index_long-index_short);
%lambda=((-(pixel-index_peak)*Q+found_peak_wavelength));   %前面的負號很重要! 影響freq domain是否接近guassian (不過為什麼差一個負號dispersion的broaden效應好像也會變大? 因為carrier in lambda domain有chirp, 但在freq domain應該不太有)

%Q=(long_coef-short_coef)/(index_long-index_short);      %roughly flens/pixel size
%lambda=grating_pitch*sin((asin((index_peak-pixel)*Q)+incidence_angle))+center_wavelength/2;

%Q=(long_coef-short_coef)/(index_long-index_short);      %roughly flens/pixel size
%lambda=grating_pitch*sin((atan((index_peak-pixel)*Q)+incidence_angle))+center_wavelength/2;   %atan!

%使用新的演算法後, 最短wavelength從166micron -> 10micron, 可能造成點數不夠, 所以先去掉不要的data

%ref=ref(lambda>400);
%lambda2=lambda2(lambda>300);
lambda_2=lambda(1400:2500);

%% x-axis
c=3E8;                     %m/sec
freq=c./(lambda_2*1E-9);     %Hz orignal freq array
d_f=max(freq)/(N_f-1);
fx=0:d_f:max(freq);        %freq after interpolation
d_t=1/(d_f*N_t);
%d_t=1/(d_f*2*N_f);
time=[-0.5*(N_t-1)*d_t:d_t:0.5*N_t*d_t]'/2;%/2是因為一來一回
%time=[-(N_f-1)*d_t:d_t:N_f*d_t]/2;  %/2是因為一來一回
length_space_FD=round(length(time)/2);
space_FD=c*time(1:length_space_FD); % 暫定只用array前1/100 (大約也有100 micron左右吧)
space_FD=(space_FD-space_FD(1))*1E6;                     % shift to zero & m to micron
%plot(lambda,inter,lambda2,inter);

clear S CS_envelope

CS_envelope(1:size(data_array,1),1:length_space_FD)=0;
S(1:size(data_array,1),1:length(fx))=0;

for j=1:size(data_array,1)

%% data related
S0=data_array(j,:)-mean(data_array(j,1:100));
%S0=S0-mean(S0(lambda>300)); %background light substration
S0=S0(1400:2500);
%S0=inter;
%plot(lambda,S0/max(S0),lambda,inter/max(inter),lambda,ref/max(ref));
%S0=S0/max(S0);        %不先normal應該也沒差?
S(j,:)=interp1(freq,S0,fx);
S(isnan(S))=0;
%CS=fft(S,N_t)';     %with minus time
%CS=real(fft(S_padded));     %with minus time
%CS_normal=CS/max(abs(CS));
%CS=CS(1:length_space_FD);
%CS_carrier=real(CS);
%CS_envelope(j,:)=abs(CS);                                        %data 1 to record (need decimate? maybe cut is better)
%% 



%plot(space,CS_normal,space,CS_envelope);

%% After the acquire of envelope


%%  DC PSF FWHM
%value_DC_FD=CS_envelope(1);
%FWHM_DC_FD=2*(space_FD(find(CS_envelope<0.5*value_DC_FD, 1, 'first'))-space_FD(1));

%%  Interfered PSF FWHM
%space_min_for_inter_peak_FD=10;  %10fixed
%space_min_for_inter_peak_index_FD=find(space_FD>space_min_for_inter_peak_FD, 1, 'first');
%[inter_peakvalue_FD inter_peakindex_FD]=max(CS_envelope(space_min_for_inter_peak_index_FD:end)); %注意! 此處的inter_peakindex_FD會抓錯, 因為CS_envelope從space_min_for_inter_peak_index_FD開始而非0開始
%inter_position_FD=space_FD(inter_peakindex_FD+space_min_for_inter_peak_index_FD-1);
%FWHM_FD_right=space_FD(find(CS_envelope(space_min_for_inter_peak_index_FD:end)>0.5*inter_peakvalue_FD, 1, 'last'));
%FWHM_FD_left=space_FD(find(CS_envelope(space_min_for_inter_peak_index_FD:end)>0.5*inter_peakvalue_FD, 1, 'first'));
%FWHM_inter_FD=FWHM_FD_right-FWHM_FD_left;
%FWHM_FD=FWHM_inter_FD;

end

CS_envelope_mean=mean(CS_envelope,1);

S_mean=mean(S,1);

CS_envelope_mean_ex(1:size(data_array,1),1:length_space_FD)=0;

S_mean_ex(1:size(data_array,1),1:length(fx))=0;

for j=1:size(data_array,1)
CS_envelope_mean_ex(j,:)=CS_envelope_mean;

S_mean_ex(j,:)=S_mean;
end

RMS_CS=sqrt(sum((CS_envelope-CS_envelope_mean_ex).^2,1)/size(data_array,1));


RMS_S=sqrt(sum((S-S_mean_ex).^2,1)/size(data_array,1));

RMS_total=sqrt(sum(RMS_S.^2));

CS_from_mean=abs(fft(S_mean,16384*2));

peak_CS=max(max(CS_from_mean));

SNR=20*log10(mean(peak_CS)/2/RMS_total);


%SNR=20*log10(4.4E7./RMS_CS);
%RMS_CS_envelope=sqrt(sum((CS_envelope-CS_mean).^2,2)/M);


%interference_efficiency=inter_peakvalue_FD/value_DC_FD*2;
%plot(space_FD,CS_normal,space_FD,CS_envelope);
%dlmwrite('Interferogram.txt',M,'delimiter','\t','newline','pc');
%BW=lambda(find(S0>0.5,1,'last'))-lambda(find(S0>0.5,1,'first'));
%x_Res=x(find(CS_envelope>0.5,1,'last'))-x(find(CS_envelope>0.5,1,'first'))
%;
%plot(space_TD,inter_position_FD,space_TD,interference_efficiency,space_TD,FWHM_FD);
%plot(space_TD,inter_position_FD,space_TD,FWHM_FD);
%interference_efficiency=interference_efficiency';
%inter_position_FD=inter_position_FD';

%imagesc(10*log10(CS_envelope./maxs)); figure(gcf)
%imagesc(10*log10(CS_envelope(2500:3500,16:116)./max(max(CS_envelope(2500:3500,16:116),[],1))),'xdata',[0:5:500],'ydata',[0:0.009375850226843:0.009375850226843*1000]); figure(gcf);
%dlmwrite('power_envelope.txt',power_envelope,'delimiter','\t','newline','pc');
%dlmwrite('power.txt',power,'delimiter','\t','newline','pc');
%dlmwrite('space_TD.txt',space_TD,'delimiter','\t','newline','pc');
%dlmwrite('CS_envelope.txt',CS_envelope,'delimiter','\t','newline','pc');
%dlmwrite('interference_efficiency.txt',interference_efficiency','delimiter','\t','newline','pc');
%dlmwrite('dispersion_expansion_ratio.txt',dispersion_expansion_ratio','delimiter','\t','newline','pc');
%dlmwrite('inter_position_FD.txt',inter_position_FD','delimiter','\t','newline','pc');