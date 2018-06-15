clear Zero_Is Zero_Ip Zero_Is_diff Zero_Ip_diff Zero_CS

tag_Is=0;
tag_Ip=0;
tag_CS=0;
x=1;
y=1;
z=1;

for j=1:(saparation_index-Starting_pixel)
    if tag_Is ==1
        if Is_interp(Starting_pixel-1+j)<0
            Zero_Is(x)=Starting_pixel-1+j;
            x=x+1;
            tag_Is=-1;
        end
    elseif tag_Is == -1
        if Is_interp(Starting_pixel-1+j)>0
            Zero_Is(x)=Starting_pixel-1+j;
            x=x+1;
            tag_Is=1;
        end
    end
        
    if tag_Ip ==1
        if Ip_interp(Starting_pixel-1+j)<0
            Zero_Ip(y)=Starting_pixel-1+j;
            y=y+1;
            tag_Ip=-1;
        end
    elseif tag_Ip == -1
        if Ip_interp(Starting_pixel-1+j)>0
            Zero_Ip(y)=Starting_pixel-1+j;
            y=y+1;
            tag_Ip=1;
        end
    end
     
      
    if tag_CS ==1
        if CS_interp(Starting_pixel-1+j)<0
            Zero_CS(z)=Starting_pixel-1+j;
            z=z+1;
            tag_CS=-1;
        end
    elseif tag_CS == -1
        if CS_interp(Starting_pixel-1+j)>0
            Zero_CS(z)=Starting_pixel-1+j;
            z=z+1;
            tag_CS=1;
        end
    end
    
    if tag_Is ==0
        if Is_interp(Starting_pixel-1+j)>0
            tag_Is=1;
        else
            tag_Is=-1;
        end
    end
    
    if tag_Ip ==0
        if Ip_interp(Starting_pixel-1+j)>0
            tag_Ip=1;
        else
            tag_Ip=-1;
        end
    end
    
    if tag_CS ==0
        if CS_interp(Starting_pixel-1+j)>0
            tag_CS=1;
        else
            tag_CS=-1;
        end
    end
end


Zero_CS_diff=diff(Zero_CS);
Zero_Is_diff=diff(Zero_Is);
Zero_Ip_diff=diff(Zero_Ip);