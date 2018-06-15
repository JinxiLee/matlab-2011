clear all;
close all;
clc;

total_OPD=115.8*2;        %micron,  *2 because now the data is for a roundtrip (uncut)
axial_resolution=1.5;   %micron


ratiop=[0:0.1:1];   
ratios=[0:0.1:1];


%ratiop=[1.2:0.01:1.3];   
%ratios=[1.7:0.01:1.8];
 
QQ=0.9:0.01:1.1;

start_position=61;
end_position=67;


array=1:80;

for jj=1:40

cd('D:\110822\10times 3\');
Is_o=importdata(sprintf('%i',array(jj)));
%Is_o=importdata('D:\110426_Original PSOCT data\PALC_DC_high
%Z_1V_1_s.txt');

Is_o=Is_o(:,1);

Is=sum(Is_o,1)/size(Is_o,1);

%plot(position,Is_o);


%Is_background=[Is(1):((Is(length(Is))-Is(1))/(length(Is)-1)):Is(end)]';

%Is=Is-Is_background;

%% Background substration
[max_value max_pixel]=max(Is_o(1:3000));

Starting_pixel=max_pixel-300;   %2Hz
Ending_pixel=max_pixel+300;

%Starting_pixel=600;   %20Hz
%Ending_pixel=1100;

Is_o(1:Starting_pixel,:)=0;

Is_o(Ending_pixel:end,:)=0;

if size(Is_o,1)==10000

Is_all(:,jj)=Is_o;

end


end


position=[0:total_OPD/(size(Is_o,1)-1):total_OPD]';  

for j=1:10000
    position(j,1:38)=position(j);
end


plot(position,Is_all);