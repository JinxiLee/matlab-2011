clear all;
close all;
clc;

total_OPD=156;
axial_resolution=1.5;   %micron

lateral_step=1;    %micron
length_lateral=1000;

window_ref=10;


axial_size=4800;

start_index_TD=1;

end_index_TD=4800;

total_OPD_new=total_OPD/axial_size*(end_index_TD-start_index_TD);

start_index_of_spectrum=100;

end_index_of_spectrum=200;

array=1;

lateral_size=1000;

    
cd('D:\');
profile_total(1:lateral_size)=0;
for jj=1:length(array)
    

    data_Bscan=importdata('1v only need calibration.txt');

    data_Bscan=data_Bscan(:,1:round(length_lateral/lateral_step));
%Is=sum(Is_o,2)/size(Is_o,2);


%% filtering and manual hilbert

    data_Bscan_f=fft(data_Bscan(start_index_TD:end_index_TD,:),[],1);

    data_Bscan_f(1:start_index_of_spectrum,:)=0;
    data_Bscan_f(end_index_of_spectrum:end,:)=0;
    %data_Bscan_f((size(data_Bscan_f,1)+1):10*(size(data_Bscan_f,1)),:)=0;
    data_Bscan_env=abs(ifft(data_Bscan_f,[],1));

axial_position=[0:total_OPD_new/(size(data_Bscan_env,1)-1):total_OPD_new]';  
%plot(axial_position,data_Bscan_env(:,1));

    data_Bscan_env(1:50,:)=0;
    data_Bscan_env((size(data_Bscan_env,1)-49:end),:)=0;

%% Finding the inerface
    [value_max index_max]=max(data_Bscan_env,[],1);
    profile_original=axial_position(index_max);

    dlmwrite(sprintf('All the same place_PROFILE_%i.txt',array(jj)),profile_original,'delimiter','\t','newline','pc');

%% FWHM calculaion
    for j=1:size(data_Bscan_env,2)
        index_left=find(data_Bscan_env(:,j)>value_max(j)/2,1,'first');
        index_right=find(data_Bscan_env(:,j)>value_max(j)/2,1,'last');
        FWHM(j)=axial_position(index_right)-axial_position(index_left);
    end
end
profile_mean=mean(profile_original);
profile_difference=profile_original-profile_mean;
reducing_raio=10;
for j=1:round(length(profile_difference)/reducing_raio)
    profile_difference_new(j)=mean(profile_difference(1+((j-1)*reducing_raio)));
end

Error=(sum((profile_difference).^2)/length(profile_difference)).^0.5;

Error_new=(sum((profile_difference_new).^2)/length(profile_difference_new)).^0.5;

imagesc(10*log10(data_Bscan_env(end:-1:1,:)),'ydata',axial_position);

plot(lateral_position,-profile_difference_new);
%% Ra, Ry, Rz, Rq calculation

%dlmwrite('D:\position.txt',position,'delimiter','\t','newline','pc');
