clear all

%% Data Loading

cd('D:\TEST\111127 - 3ms 100ave 3mm v0.01 Green - 3\');
%Data=importdata('111010_Green (2500microsec) no word 5 ave 100.txt');
n_final=importdata('n_final.txt');
k_final=importdata('k_final.txt');
Thickness_temp=importdata('Thickness_temp.txt');
Thickness=importdata('Thickness.txt');
Frequency=importdata('Frequency.txt');
Data_Spectroscopy=importdata('111118_Green 5-1.jws.txt');

Spectrum_abs=importdata('Spectrum_abs.txt');

Spectrum_real=importdata('Spectrum_real.txt');
C=3E8;

Spectroscopy_Old=Data_Spectroscopy(end:-1:1,2)/100;
Wavelength_Spectroscopy=Data_Spectroscopy(end:-1:1,1)/1000;     

n2=1;

% Assumed n2 = BK7
C1 = 1.03961212; 
C2 = 0.00600069867; 
C3 = 0.231792344; 
C4 = 0.0200179144; 
C5 = 1.01046945; 
C6 = 103.560653;
for j=1:length(Frequency)
    Frequency(j,1:size(n_final,2))=Frequency(j);
end
for j=1:length(Thickness_temp)
    Thickness(1:size(n_final,1),j)=Thickness_temp(j);
end

Wavelength_micron=(C./Frequency)*1E6;
n_bk7=(C1*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C2)+C3*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C4)+C5*(Wavelength_micron.^2)./((Wavelength_micron.^2)-C6)+1).^0.5;

n_bk7=abs(n_bk7);

n_bk7(isnan(n_bk7))=0;
n1=n_bk7;

r_BK7=((n_bk7-1)./(n_bk7+1));
t_AN100=(2*(1)./(n_bk7+1));
            % About the conditions
        
        r1=(n1-(n_final+i*k_final))./(n_final+n1+i*k_final);
        r1_r=((n_final+i*k_final)-n1)./(n_final+n1+i*k_final);    
        t1=2*(n1)./(n_final+n1+i*k_final);
        t1_r=2*(n_final+i*k_final)./(n_final+n1+i*k_final);
        t2=2*(n_final+i*k_final)./(n_final+n2+i*k_final);
        r2=((n_final+i*k_final)-n2)./((n_final+i*k_final)+n2);   
        d=exp(i*2*pi.*Frequency/C.*(n_final+i*k_final).*Thickness);   %注意! -1*n!
        
        d_n=exp(i*2*pi.*Frequency/C.*(n_final).*Thickness);   %注意! -1*n!
        
        d_k=exp(i*2*pi.*Frequency/C.*(i*k_final).*Thickness);   %注意! -1*n!
        
        Spectrum_Amplitude_Temp=abs((r1+t1.*t1_r.*r2.*(d.^2))./r_BK7);                       %神說是cosine
        %Spectroscopy_Temp=abs((t_AN100(index_now).^2).*(t1.*t2.*d)./(1+r1_r.*r2.*(d.^2))).^2;
        T_Temp=abs((t_AN100).*(t1.*t2.*d)./(1+r1_r.*r2.*(d.^2))).^2;
        
          
        plot(Wavelength_micron,Spectrum_abs,Wavelength_micron,Spectrum_Amplitude_Temp);
        plot(Wavelength_micron,real((d_n.^2)/max(d_n.^2)),Wavelength_micron,(Spectrum_real/max(Spectrum_abs(4095:4880))));
        
        plot(Wavelength_micron,real((d_n.^2)/max(d_n.^2)),Wavelength_micron,real((d_k.^2)/max(d_k.^2)));
        
        plot(Wavelength_micron,T_Temp,Wavelength_Spectroscopy,Spectroscopy_Old);
%%
        K=4;
        plot(Wavelength_micron(Wavelength_micron(:,K)<0.6,K),T_Temp(Wavelength_micron(:,K)<0.6,K),Wavelength_Spectroscopy,Spectroscopy_Old);
        
dlmwrite('T_Temp.txt',T_Temp,'delimiter','\t','newline','pc');