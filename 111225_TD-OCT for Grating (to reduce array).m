clear all;
close all;
clc;

total_OPD=156;        %micron ask ¾G¤D¹Å
axial_resolution=1.5;   %micron

lateral_step=1;    %micron
length_lateral=1000;

window_ref=10;

array=1:1;

cd('D:\111229\Grating (100micron) step  5micron total 200points_1\');
start_index=1501;
end_index=2500;
start_index_ref=1001;
end_index_ref=2000;
for jj=1:length(array)  
    cd('D:\111229\Grating (100micron) step  5micron total 200points_1\');
    data_Bscan=importdata(sprintf('Grating (100micron) step  5micron total 200points_%i.txt',array(jj))); 
    cd('D:\111229\Grating (100micron) step  5micron total 200points_1\cut\');
    dlmwrite(sprintf('Grating (100micron) step  5micron total 200points_%i.txt',array(jj)),data_Bscan(start_index:end_index,:),'delimiter','\t','newline','pc');
end

cd('D:\111229\Grating (100micron) step  5micron total 200points_1\');
ref_Bscan=importdata('mirror 5micron total 200points  ave1.txt');
cd('D:\111229\Grating (100micron) step  5micron total 200points_1\cut\');
dlmwrite('mirror 5micron total 200points  ave1.txt',ref_Bscan(start_index_ref:end_index_ref,:),'delimiter','\t','newline','pc');
total_OPD_new=total_OPD/size(data_Bscan,1)*(end_index_ref-start_index_ref+1);
