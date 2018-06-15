clear all;
close all;
clc;

total_OPD=32.5325;
axial_resolution=1.5;   %micron

lateral_step=5;    %micron
length_lateral=1000;

window_ref=10;

array=1;

cd('D:\111229\Grating (100micron) step  5micron total 200points_1\cut\');

for jj=1:length(array)

data_Bscan=importdata(sprintf('Grating (100micron) step  5micron total 200points_%i.txt',array(jj)));
ref_Bscan=importdata('mirror 5micron total 200points  ave1.txt');
data_Bscan=data_Bscan(:,1:round(length_lateral/lateral_step));
ref_Bscan=ref_Bscan(:,1:round(length_lateral/lateral_step));
%Is=sum(Is_o,2)/size(Is_o,2);

axial_position=[0:total_OPD/(size(data_Bscan,1)-1):total_OPD]';  

lateral_position=[0:lateral_step:lateral_step*(size(data_Bscan,2)-1)]';  

%% filtering and manual hilbert

data_Bscan_f=fft(data_Bscan,[],1);
ref_Bscan_f=fft(ref_Bscan,[],1);

start_index_of_spectrum=50;

end_index_of_spectrum=170;

data_Bscan_f(1:start_index_of_spectrum,:)=0;
data_Bscan_f(end_index_of_spectrum:end,:)=0;
data_Bscan_env=abs(ifft(data_Bscan_f,[],1));


ref_Bscan_f(1:start_index_of_spectrum,:)=0;
ref_Bscan_f(end_index_of_spectrum:end,:)=0;
ref_Bscan_env=abs(ifft(ref_Bscan_f,[],1));

%plot(axial_position,data_Bscan_env(:,1));

data_Bscan_env(1:50,:)=0;
data_Bscan_env((size(data_Bscan_env,1)-49:end),:)=0;

ref_Bscan_env(1:50,:)=0;
ref_Bscan_env((size(ref_Bscan_env,1)-49:end),:)=0;

%% Finding the inerface
[value_max index_max]=max(data_Bscan_env,[],1);
profile=axial_position(index_max);
plot(lateral_position,profile);


[value_max_ref index_max_ref]=max(ref_Bscan_env,[],1);
profile_ref=axial_position(index_max_ref);
plot(lateral_position,profile_ref);

%% smooth ref
profile_ref=smooth(profile_ref,window_ref);
profile_ref=circshift(profile_ref,-10);
%% Sub ref

profile_calibrated=profile-profile_ref;

plot(lateral_position,profile_calibrated);

%% FWHM calculaion
for j=1:size(data_Bscan_env,2)
    index_left=find(data_Bscan_env(:,j)>value_max(j)/2,1,'first');
    index_right=find(data_Bscan_env(:,j)>value_max(j)/2,1,'last');
    FWHM(j)=axial_position(index_right)-axial_position(index_left);
end

%% Grating

Area_1_left=127;  %index
Area_1_right=134;  %index

Area_2_left=137;  %index
Area_2_right=144;  %index

%% Reference points

Position_1=10;

Position_2=190;

%% To solve Obliquity
angle=atan((profile_calibrated(Position_1)-profile_calibrated(Position_2))/(lateral_position(Position_1)-lateral_position(Position_2)))/pi*180; %-0.48 for small 2, 
for j=1:size(profile_calibrated,1)
    profile_tilted(j)=profile_calibrated(j)-lateral_position(j)*tan(angle*pi/180);
end
    % Height 1
Height_1=mean(profile_tilted(Area_1_left:Area_1_right));
Height_2=mean(profile_tilted(Area_2_left:Area_2_right));
Step_difference(jj)=Height_1-Height_2;


end
plot(lateral_position,profile_calibrated,lateral_position,profile_tilted);

plot(lateral_position,FWHM);
imagesc(10*log10(data_Bscan_env(end:-1:1,:)),'xdata',lateral_position,'ydata',axial_position);
%% Ra, Ry, Rz, Rq calculation

%dlmwrite('D:\position.txt',position,'delimiter','\t','newline','pc');
