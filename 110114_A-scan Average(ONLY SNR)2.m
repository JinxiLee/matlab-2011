clear all;
close all;
clc;

N_f=4096*4;
N_t=N_f;

%% To calculate SNR with envelope (assume incoherent summation between
%% frames)

data_array_all=importdata('D:\110117_ascan averge1.txt');

clear   data_array averaging_factor

data_array_all=data(1500:end,:);


averaging_factor=500;

data_array(1:fix(size(data_array_all,1)/averaging_factor),size(data_array_all,2))=0;

for j=1:fix(size(data_array_all,1)/averaging_factor)
data_array(j,:)=sum(data_array_all(1+(j-1)*averaging_factor:j*averaging_factor,:),1);
end

clear data_array_all

if size(data_array,1)>100
data_array=data_array(1:100,:);
end









clear S CS_envelope S_mean_ex

S=data_array-mean(data_array(1,1:100));
S=S(:,1000:3000);

S_mean=mean(S,1);
S_mean_ex(1:size(S,1),1:size(S,2))=0;

for j=1:size(S,1)
S_mean_ex(j,:)=S_mean;
end

RMS_S=sqrt(sum((S-S_mean_ex).^2,1)/size(S,1));

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