% Define the grid and find the distance of the center of each cell unit
% from events in seq
function [Ggrid,rxy,dA,Xcgrid,Ycgrid,Ngrid]=Bayesian_SetGrid()
global longitude latitude lonMin lonMax latMin latMax deltaGrid
Xgrid = lonMin:deltaGrid:lonMax;
Ygrid = latMin:deltaGrid:latMax;

Xcgrid = Xgrid(1:end-1)+deltaGrid/2;
Ycgrid = Ygrid(1:end-1)+deltaGrid/2;

Ngrid = length(Xcgrid)*length(Ycgrid);

Ggrid = zeros(Ngrid,2);
count=1;
for j=1:length(Ycgrid)
    for k=1:length(Xcgrid)
        Ggrid(count,:) = [Ycgrid(j),Xcgrid(k)];
        count=count+1;
    end
end
rxy = calculate_rxy(latitude,longitude,Ggrid);   % rxy = [length(Seq),length(Ggrid)];
% rxy is distance to each Grids like 370X6e4
dA = prod(topgeo(Ygrid(2)*pi/180,Xgrid(2)*pi/180,Ygrid(1)*pi/180,Xgrid(1)*pi/180));
end