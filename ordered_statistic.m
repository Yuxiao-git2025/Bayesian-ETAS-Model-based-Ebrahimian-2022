% Function for Calculating ordered statistics
function [xmean,x50,x16,x84,x02,x98] = ordered_statistic(x)
% <A>
% What we input is sum of every column, then we average the row, get xmean
% <B>
% What we input is initial samples(like 6e4*371), then we average it, get
% xmean, represented each grid points
xmean = mean(x,2);
x50 = median(x,2);
Ns = size(x,2);
x02 = (interp1((1:Ns)',sort(x,2)',normcdf(-2)*Ns,'nearest'))';
x16 = (interp1((1:Ns)',sort(x,2)',normcdf(-1)*Ns,'nearest'))';
x84 = (interp1((1:Ns)',sort(x,2)',normcdf(1)*Ns,'nearest'))';
x98 = (interp1((1:Ns)',sort(x,2)',normcdf(2)*Ns,'nearest'))';   
end




