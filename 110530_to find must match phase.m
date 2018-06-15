%pos=position(9000:11000);
%ap=angp(9000:11000);
%as=angs(9000:11000);

pos=position(1:4000);
ap=angp(1:4000);
as=angs(1:4000);


for j=1:200
as_shift=circshift(as,-j);
delta_phi=ap-as_shift;
delta_phi=delta_phi(500:length(delta_phi)-500);
delta_variation(j)=max(delta_phi)-min(delta_phi);
end

as_shift=circshift(as,0);
delta_phi=ap-as_shift;
delta_variation(j)=max(delta_phi)-min(delta_phi);
%as_shift=circshift(as,-10);

%delta_phi=ap-as_shift;
plot(1:150,delta_variation(1:150));