clear all;

lateral_res=100;     %micron
space=0:0.1:500;     %micron
spot=gaussmf(space,[lateral_res/2/(2*log(2))^0.5 250]);
plot(space,spot);
integrated_spot(1:length(spot))=0;
for j=1:length(spot)
    integrated_spot(j)=sum(spot(1:j));
end
norm=integrated_spot/max(integrated_spot);
plot(space,norm);
perc=norm(space==lateral_res);