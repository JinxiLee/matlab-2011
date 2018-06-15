clear all;

delta_theta=-1:0.001:1;
output=cos(delta_theta).^2;
%FWHA_2=delta_theta(find(output>10^(-0.1),1,'first'))*180/pi;
FWHA_2=delta_theta(find(output>10^(-0.3),1,'first'))*180/pi;
deltaV_V=FWHA_2/45;

BW_ASE=180;
