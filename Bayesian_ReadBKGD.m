function [Muxy,Muxy2col,Nbo]=Bayesian_ReadBKGD(Xcgrid,Ycgrid,T_start,Ngrid,dA,UseBKGD)
% Read Background Seismicity

if ~UseBKGD
    Muxy = zeros(length(Ycgrid),length(Xcgrid));
else
    Muxy = load('background seismicity.txt');
end

% Background seismicity becomes a column, mu_xy_Ml_t is the transformed mu_xy_Ml
Muxy2col = zeros(Ngrid,1);
count=1;
for j=1:length(Ycgrid)
    for k=1:length(Xcgrid)
        Muxy2col(count) = Muxy(j,k);
        count=count+1;
    end
end

Nbo = sum(sum(Muxy))*dA*T_start;   
end