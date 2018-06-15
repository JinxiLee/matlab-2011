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


array=1:100;

for jj=1:50

%cd('D:\110827\2Hz\');     %4 (2Hz)
cd('D:\110915\2Hz - after Green 2\');     %4 (2Hz)
Is_o=importdata(sprintf('%i',array(jj)));
%Is_o=importdata('D:\110426_Original PSOCT data\PALC_DC_high
%Z_1V_1_s.txt');

if size(Is_o,1) ~= 10000
    continue;
end

Is_o=Is_o(:,1);

Is=sum(Is_o,1)/size(Is_o,1);

position=[0:total_OPD/(length(Is)-1):total_OPD]';  

%plot(position,Is_o);


%Is_background=[Is(1):((Is(length(Is))-Is(1))/(length(Is)-1)):Is(end)]';

%Is=Is-Is_background;

%% Background substration
[max_value max_pixel]=max(Is_o(1:4500));        %3000 for 20Hz, 4500 for 2Hz

Starting_pixel=max_pixel-300;   %2Hz
Ending_pixel=max_pixel+300;

%Starting_pixel=600;   %20Hz
%Ending_pixel=1100;

Is_o(1:Starting_pixel,:)=0;

Is_o(Ending_pixel:end,:)=0;

%% Zero padding and filtering

Starting_pixel_f=500;

Ending_pixel_f=1200;


Is_o_ftemp=fft(Is_o,[],1);

Is_o_ftemp(1:Starting_pixel_f,:)=0;

Is_o_ftemp((length(Is_o_ftemp)-Starting_pixel_f+1):length(Is_o_ftemp),:)=0;

Is_o_ftemp(Ending_pixel_f:(length(Is_o_ftemp)-Ending_pixel_f),:)=0;


Is_o_ftemp1=Is_o_ftemp(1:round(length(Is_o_ftemp)/2),:);
Is_o_ftemp2=Is_o_ftemp((round(length(Is_o_ftemp)/2)+1):end,:);
Is_o_ftemp=[Is_o_ftemp1;zeros(length(Is_o_ftemp)*7,size(Is_o_ftemp,2)); Is_o_ftemp2];
Is_o_new=real(ifft(Is_o_ftemp,[],1));
Is_o_ftemp_for_hilbert=2*Is_o_ftemp;
Is_o_ftemp_for_hilbert((round(size(Is_o_ftemp_for_hilbert,1)/2)+1):end,:)=0;
Is_o_new_envelope=abs(ifft(Is_o_ftemp_for_hilbert,[],1));
%Is_o_new=ifftshift(real(ifft(abs(Is_o_ftemp))));

position=[0:total_OPD/(length(Is_o_new)-1):total_OPD]';  

%plot()

%plot(position,Is_o_new,position,Is_o_new_envelope);
%% lambda genaration

c=3E8;

dx=(total_OPD*(1E-6))/(size(Is_o_new,1)-1);    %m

dt=2*dx/c;    %this is for real OPD

f_total=1/dt;

frequency=1:f_total/(size(Is_o_new,1)-1):f_total;

wavelength=c./frequency'*1E6;

Is_o_ftemp_abs=abs(Is_o_ftemp);

%plot(wavelength(500:end),abs(Is_o_ftemp(500:end,:)));

%plot(frequency,abs(Is_o_ftemp));

%% to find max and align

[maxvalue maxindex]=max(Is_o_new,[],1);
meanmaxindex=mean(maxindex);
needshift=round((round(size(Is_o_new,1))/2-maxindex));


%for j=1:size(Is_o_new,2)
    %Is_o_new(:,j)=circshift(Is_o_new(:,j),needshift(j));
%end



%% to find zeros
[max_value_n max_pixel_n]=max(Is_o_new(:,1));
Starting_pixel_n=max_pixel_n-1000;

Ending_pixel_n=max_pixel_n+1000;

%Starting_pixel_n=39000;

%Ending_pixel_n=41000;


for k=1:size(Is_o_new,2)
tag_Is=0;
tag_Ip=0;
x=1;
y=1;
for j=1:(Ending_pixel_n-Starting_pixel_n)
    if tag_Is ==1
        if Is_o_new(Starting_pixel_n-1+j,k)<0
            Zero_Is(x,k)=Starting_pixel_n-1+j;
            x=x+1;
            tag_Is=-1;
        end
    elseif tag_Is == -1
        if Is_o_new(Starting_pixel_n-1+j,k)>0
            Zero_Is(x,k)=Starting_pixel_n-1+j;
            x=x+1;
            tag_Is=1;
        end
    end
        
     
    if tag_Is ==0
        if Is_o_new(Starting_pixel_n-1+j,k)>0
            tag_Is=1;
        else
            tag_Is=-1;
        end
    end
    
end
end

Zero_Is_diff=diff(Zero_Is,1);

Zero_Is_diff(abs(Zero_Is_diff)>1000)=0;

Zero_Is(Zero_Is==0)=1;

%plot(position(Zero_Is(1:(size(Zero_Is,1)-1),:)),1./(2*(Zero_Is_diff*total_
%OPD/(length(Is_o_new)-1))*1E-6));

Zero_Diff_position=position(Zero_Is(1:(size(Zero_Is,1)-1),:));

nearest_max_Zero_index=find(Zero_Is>max_pixel_n,1,'first');

considered_number_of_Zeros_half=2;

considered_Zeros=Zero_Is((nearest_max_Zero_index-considered_number_of_Zeros_half):(nearest_max_Zero_index+considered_number_of_Zeros_half));

considered_Zeros_Diff=diff(considered_Zeros);

considered_Zeros_Diff_position=position(considered_Zeros(1:(size(considered_Zeros,1)-1),:));

ave_considered_Zeros_Diff=mean(considered_Zeros_Diff);

ave_considered_Zeros_Diff_position=mean(position(considered_Zeros(1:(size(considered_Zeros,1)-1),:)));

jj_pos(jj)=ave_considered_Zeros_Diff_position;

jj_pos2(((jj-1)*(2*considered_number_of_Zeros_half)+1):jj*(2*considered_number_of_Zeros_half))=considered_Zeros_Diff_position';

jj_diff(jj)=ave_considered_Zeros_Diff;

jj_diff2(((jj-1)*(2*considered_number_of_Zeros_half)+1):jj*(2*considered_number_of_Zeros_half))=considered_Zeros_Diff';

%plot(considered_Zeros_Diff_position,considered_Zeros_Diff);

end
jj_pos=jj_pos';
jj_diff=jj_diff';

plot(jj_pos,jj_diff);


%% to calculate average micron/pixel

known_center_wavelength=0.56;       %micron

known_zeros_difference=known_center_wavelength/4;

calculate_total_length=known_zeros_difference/ave_considered_Zeros_Diff*length(position)/2;