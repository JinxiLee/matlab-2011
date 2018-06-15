clear all;
close all;
clc;

total_OPD=156;        %micron ask ¾G¤D¹Å
axial_resolution=1.5;   %micron

lateral_step=1;    %micron
length_lateral=750;

window_ref=10;

data_Bscan=importdata('D:\111224\111224_LARGE SAMPLE 1000micron ave10 2.txt');
ref_Bscan=importdata('D:\111224\111224_mirror 100micron ave1 2 (for large sample).txt');
data_Bscan=data_Bscan(:,1:round(length_lateral/lateral_step));
ref_Bscan=ref_Bscan(:,1:round(length_lateral/lateral_step));
%Is=sum(Is_o,2)/size(Is_o,2);

axial_position=[0:total_OPD/(size(data_Bscan,1)-1):total_OPD]';  

lateral_position=[0:lateral_step:lateral_step*(size(data_Bscan,2)-1)]';  

%% filtering and manual hilbert

data_Bscan_f=fft(data_Bscan,[],1);
ref_Bscan_f=fft(ref_Bscan,[],1);

start_index_of_spectrum=350;

end_index_of_spectrum=800;

data_Bscan_f(1:start_index_of_spectrum,:)=0;
data_Bscan_f(end_index_of_spectrum:end,:)=0;
data_Bscan_env=abs(ifft(data_Bscan_f,[],1));


ref_Bscan_f(1:start_index_of_spectrum,:)=0;
ref_Bscan_f(end_index_of_spectrum:end,:)=0;
ref_Bscan_env=abs(ifft(ref_Bscan_f,[],1));

plot(axial_position,data_Bscan_env(:,1));

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
%% Sub ref

profile_calibrated=profile-profile_ref;

plot(lateral_position,profile_calibrated);

%% FWHM calculaion
for j=1:size(data_Bscan_env,2)
    index_left=find(data_Bscan_env(:,j)>value_max(j)/2,1,'first');
    index_right=find(data_Bscan_env(:,j)>value_max(j)/2,1,'last');
    FWHM(j)=axial_position(index_right)-axial_position(index_left);
end

angle=0.12; %-0.48 for small 2, 
ratio=1;
for b=1:size(ratio,2)
for q=1:size(angle,2)
current_sign=1;
temp_max=0;
temp_min=0;
number_of_maxs=1;
number_of_mins=1;
maxs=0;
mins=0;
profile_calibrated=profile-profile_ref*ratio(b);
    for j=1:size(profile_calibrated,1)
        profile_tilted(j)=profile_calibrated(j)-lateral_position(j)*tan(angle(q)*pi/180);
    end
    profile_sub(:,q)=profile_tilted-mean(profile_tilted);
    Ra(b,q)=sum(abs(profile_sub(:,q)))/size(profile_sub(:,q),1);
    Rq(b,q)=(sum((profile_sub(:,q)).^2)/size(profile_sub(:,q),1)).^0.5;
    Ry(b,q)=max(profile_sub(:,q))-min(profile_sub(:,q));
    %Rz
    for p=1:size(profile_sub(:,q),1)
        if p == 1
            if profile_sub(1,q) >= 0
                current_sign=1;
            else
                current_sign=-1;
            end
        end
        if current_sign == 1
            if profile_sub(p,q) > temp_max
                temp_max=profile_sub(p,q);
            end
            if profile_sub(p,q) < 0
                maxs(number_of_maxs)=temp_max;
                temp_max=0;
                number_of_maxs=number_of_maxs+1;
                current_sign=-1;
                if profile_sub(p,q) < temp_min
                    temp_min=profile_sub(p,q);
                end
            end
        else
            if profile_sub(p,q) < temp_min
                temp_min=profile_sub(p,q);
            end
            if profile_sub(p,q) > 0
                mins(number_of_mins)=temp_min;
                temp_min=0;
                number_of_mins=number_of_mins+1;
                current_sign=1;
                if profile_sub(p,q) > temp_max
                    temp_max=profile_sub(p,q);
                end
            end
        end
    end
    maxs_sorted=sort(maxs,'descend');
    mins_sorted=sort(mins);
    if size(maxs_sorted,2)>=5 && size(mins_sorted,2)>=5
        Rz(b,q)=(mean(maxs_sorted(1:5))+abs(mean(mins_sorted(1:5))));
    else
        Rz(b,q)=99999999;
    end
end
end

[value_minra index_minra_1]=min(Ra);

[value_minra index_minra_2]=min(value_minra);
index_minra_1=index_minra_1(index_minra_2);

ratio_final=ratio(index_minra_1);

angle_final=angle(index_minra_2);

plot(lateral_position,profile_calibrated,lateral_position,profile_tilted);

plot(lateral_position,profile_sub);

plot(lateral_position,FWHM);
imagesc(10*log10(data_Bscan_env(end:-1:1,:)),'xdata',lateral_position,'ydata',axial_position);
%% Ra, Ry, Rz, Rq calculation

%dlmwrite('D:\position.txt',position,'delimiter','\t','newline','pc');
