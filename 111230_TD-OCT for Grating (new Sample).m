clear all;
close all;
clc;

total_OPD=156;
axial_resolution=1.5;   %micron

lateral_step=1;    %micron
length_lateral=200;

window_ref=10;


start_index_of_spectrum=350;

end_index_of_spectrum=800;

array=1:200;

axial_size=4800;
lateral_size=round(length_lateral/lateral_step);

axial_position=[0:total_OPD/(axial_size-1):total_OPD]';  

lateral_position=[0:lateral_step:lateral_step*(lateral_size-1)]';  
    
cd('D:\grating again (step 1 micron total 200 micron)_1\');
profile_total(1:lateral_size)=0;
profile_original(1:lateral_size,1:length(array))=0;
profile_tilted(1:lateral_size,1:length(array))=0;
for jj=1:length(array)
    
if exist('ref_PROFILE.txt')
    profile_ref=importdata('ref_PROFILE.txt');
else
    
    ref_Bscan=importdata('ref.txt');  
    ref_Bscan=ref_Bscan(:,1:round(length_lateral/lateral_step)); 
    
    
    ref_Bscan_f=fft(ref_Bscan,[],1);
    
    ref_Bscan_f(1:start_index_of_spectrum,:)=0;
    ref_Bscan_f(end_index_of_spectrum:end,:)=0;
    ref_Bscan_env=abs(ifft(ref_Bscan_f,[],1));
    
    ref_Bscan_env(1:50,:)=0;
    ref_Bscan_env(((size(ref_Bscan_env,1)-49):end),:)=0;
    
    [value_max_ref index_max_ref]=max(ref_Bscan_env,[],1);
    profile_ref=axial_position(index_max_ref);

%% smooth ref
    profile_ref=smooth(profile_ref,window_ref);
    dlmwrite('ref_PROFILE.txt',profile_ref,'delimiter','\t','newline','pc');
end
    
    
    
if exist(sprintf('sam_PROFILE_%i.txt',array(jj)))
    profile_original(:,jj)=importdata(sprintf('sam_PROFILE_%i.txt',array(jj)));
else

    data_Bscan=importdata(sprintf('grating again (step 1 micron total 200 micron)_%i.txt',array(jj)));

    data_Bscan=data_Bscan(:,1:round(length_lateral/lateral_step));
%Is=sum(Is_o,2)/size(Is_o,2);


%% filtering and manual hilbert

    data_Bscan_f=fft(data_Bscan,[],1);

    data_Bscan_f(1:start_index_of_spectrum,:)=0;
    data_Bscan_f(end_index_of_spectrum:end,:)=0;
    data_Bscan_env=abs(ifft(data_Bscan_f,[],1));

%plot(axial_position,data_Bscan_env(:,1));

    data_Bscan_env(1:50,:)=0;
    data_Bscan_env((size(data_Bscan_env,1)-49:end),:)=0;

%% Finding the inerface
    [value_max index_max]=max(data_Bscan_env,[],1);
    profile_original(:,jj)=axial_position(index_max);

    dlmwrite(sprintf('sam_PROFILE_%i.txt',array(jj)),profile_original(:,jj),'delimiter','\t','newline','pc');
    
    %% FWHM calculaion
for j=1:size(data_Bscan_env,2)
    index_left=find(data_Bscan_env(:,j)>value_max(j)/2,1,'first');
    index_right=find(data_Bscan_env(:,j)>value_max(j)/2,1,'last');
    FWHM(j)=axial_position(index_right)-axial_position(index_left);
end


    dlmwrite(sprintf('sam_FWHM_%i.txt',array(jj)),profile_original(:,jj),'delimiter','\t','newline','pc');
end
%% Sub ref
profile_calibrated=profile_original(:,jj)-profile_ref;



%% Grating

Area_1_left=30;  %index
Area_1_right=55;  %index

Area_2_left=80;  %index
Area_2_right=105;  %index

%% Reference points

Position_1=1;

Position_2=199;

%% To solve Obliquity
angle=atan((profile_calibrated(Position_1)-profile_calibrated(Position_2))/(lateral_position(Position_1)-lateral_position(Position_2)))/pi*180; %-0.48 for small 2, 
for j=1:size(profile_calibrated,1)
    profile_tilted(j,jj)=profile_calibrated(j)-lateral_position(j)*tan(angle*pi/180);
    profile_total(j)=profile_total(j)+profile_tilted(j,jj);
end
    % Height 1
Height_1=mean(profile_tilted(Area_1_left:Area_1_right,jj));
Height_2=mean(profile_tilted(Area_2_left:Area_2_right,jj));
Step_difference(jj)=Height_1-Height_2;

end

Step_difference_sub=Step_difference-mean(Step_difference);
Step_difference_mean=mean(Step_difference);
reducing_raio=10;


for j=1:round(length(Step_difference_sub)/reducing_raio)
    Step_difference_sub_new(j)=mean(Step_difference_sub(1+((j-1)*reducing_raio)));
end
Error=(sum((Step_difference_sub).^2)/size(Step_difference_sub,1)).^0.5;
Error_new=(sum((Step_difference_sub_new).^2)/length(Step_difference_sub_new)).^0.5;
profile_mean=profile_total/length(array);
plot(lateral_position,profile_tilted);

plot(lateral_position,FWHM);

for j=1:200
    data_Bscan_env_shifted=circshift(data_Bscan_env,1540-round(profile_ref(j)/156*4800));
end
imagesc(10*log10(data_Bscan_env_shifted(end:-1:1,:)),'xdata',lateral_position,'ydata',axial_position-50);
%% Ra, Ry, Rz, Rq calculation

%dlmwrite('D:\position.txt',position,'delimiter','\t','newline','pc');
