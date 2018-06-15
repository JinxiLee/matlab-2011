clear all;
close all;
clc;

total_OPD=115.8;        %micron
axial_resolution=1.5;   %micron


data_o=importdata('D:\110426_Original PSOCT data\PALC_0V_1.txt');

data=sum(data_o,2)/size(data_o,2);

clear data_o;

position=[total_OPD/length(data):total_OPD/length(data):total_OPD]';

pixel_size=total_OPD/length(data);   %micron

spatial_PSF=gaussmf(position,[axial_resolution 0]);

spectrum=fft(spatial_PSF);

plot(position,data);

FFTdata=fft(data);

ratio=0.05;

FFTdata(find(spectrum<ratio*max(spectrum),1,'first'):find(spectrum<ratio*max(spectrum),1,'last'))=0;

index_highest_spatial_frequency=find(spectrum<ratio*max(spectrum),1,'first');

filtered_data=ifft(FFTdata);

d_filtered_data=diff(filtered_data);

d_filtered_data(length(filtered_data))=d_filtered_data(length(filtered_data)-1);

%plot(position,d_filtered_data*200,position,filtered_data);

%plot(1:length(position),d_filtered_data);

plot(1:length(position),data);