clear all

%% %% 1st part: to find the calibration array

center_wavelength=0.56;     %micron

zero_difference=center_wavelength/4;

array=1:100;

N_total_scan=50;

N_considered_zero=2;        %2 for 5, 3 for 7, 4 for 9 ...

%index_array(1:N_total_scan)=0;      %the position "index" of each array
index_array=1;
%speed_array(1:N_total_scan)=0;      %the ratio between index and position (delta_position(micron)/delta_index(namely, 1))
speed_array=0;
%position_array(1:N_total_scan)=0;      %NOT USING THIS, the position array should be calculated after interpolation

IfContinoue=1;
N_scan=0;

while (IfContinoue>0)

    N_scan=N_scan+1;

    cd('D:\110915\2Hz - after Green 2\');
    IfContinoue=exist(sprintf('%i',array(N_scan)),'file');

    if (IfContinoue == 0)
        continue;
    end


    
%% Reference Data loading
    
    cd('D:\110915\2Hz - after Green 2\');
    Signal=importdata(sprintf('%i',array(N_scan)));
    Signal=Signal(:,1);
%Signal(50000:end)=Signal(50000);
    Signal(4500:end)=Signal(4500);
    
%% DC filtering (to avoid gap bwtween set zero and raw-data DC)

%Starting_pixel_f=100;
%Ending_pixel_f=600;


    Starting_pixel_f=600;
    Ending_pixel_f=1600;

    Signal_temp_f=fft(Signal);

    Signal_temp_f(1:Starting_pixel_f)=0;
    Signal_temp_f(Ending_pixel_f:Starting_pixel_f:(length(Signal_temp_f)-Ending_pixel_f+1))=0;
    Signal_temp_f((length(Signal_temp_f)-Starting_pixel_f+1):length(Signal_temp_f),:)=0;

    Signal=real(ifft(Signal_temp_f));

%% Peak finding

    [maxvalue maxindex]=max(Signal,[],1);
    if index_array(length(index_array)) == maxindex
        continue;
    end
    
    index_array(length(index_array)+1)=maxindex;

%%

% if we know there's N pixel between window, then the speed equal to zero_difference / N
% speed: the position differene of adjacent pixel
% so after the calibration, we will have a power array with non-uniform
% difference, this is where interpolation walks in


%Starting_pixel=maxindex-1000;

%Ending_pixel=maxindex+1000;

    Starting_pixel=maxindex-50;

    Ending_pixel=maxindex+50;

    tag_Is=0;
    tag_Ip=0;
    x=1;
    y=1;

    for j=1:(Ending_pixel-Starting_pixel)
        if tag_Is ==1
            if Signal(Starting_pixel-1+j)<0
                Zeros(x)=Starting_pixel-1+j;
                x=x+1;
                tag_Is=-1;
            end
        elseif tag_Is == -1
            if Signal(Starting_pixel-1+j)>0
                Zeros(x)=Starting_pixel-1+j;
                x=x+1;
                tag_Is=1;
            end
        end
        
     
        if tag_Is ==0
            if Signal(Starting_pixel-1+j)>0
                tag_Is=1;
            else
                tag_Is=-1;
            end
        end
    
    end

    Zeros_diff=diff(Zeros);

    Zeros(Zeros==0)=1;
    Zeros=Zeros(1:(length(Zeros)-1));

    nearest_max_Zero_index=find(Zeros>maxindex,1,'first');

    Zeros_considered=Zeros((nearest_max_Zero_index-N_considered_zero):(nearest_max_Zero_index+N_considered_zero));
    Zeros_diff_considered=Zeros_diff((nearest_max_Zero_index-N_considered_zero):(nearest_max_Zero_index+N_considered_zero));

    ave_Zeros_diff_considered=mean(Zeros_diff_considered);

    speed_array(length(speed_array)+1)=zero_difference/ave_Zeros_diff_considered;

end

%% to generate the new position array

index_array(1)=1;
index_array(length(index_array)+1)=length(Signal);

speed_array(1)=speed_array(2);
speed_array(length(speed_array)+1)=speed_array(length(speed_array));
%plot(index_array,speed_array);

speed_array_interp=interp1(index_array,speed_array,1:length(Signal));

plot(1:length(Signal),speed_array_interp,index_array,speed_array);

for j=1:length(speed_array_interp)
    position_array_interp(j)=sum(speed_array_interp(1:j));
end

%% %% 2nd part: to generate the new PSF data

cd('D:\110915\2Hz - after Green 2\');
Signal=importdata(sprintf('%i',array(N_scan)));