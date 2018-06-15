clear all;
close all;
clc;

lateral_position=(1:750)';  


profile_original=importdata('D:\small2alpha.txt');
profile=interp1(profile_original(:,1),profile_original(:,2),(1:750))';
angle=0; %-0.48 for small 2, 
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
profile_calibrated=profile;
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

%% Ra, Ry, Rz, Rq calculation

%dlmwrite('D:\position.txt',position,'delimiter','\t','newline','pc');
