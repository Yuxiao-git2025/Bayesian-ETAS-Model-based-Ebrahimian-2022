% =========================================================================
%                       Written by Yu.2025/9/13
%                        Refer to <Zechar.2010>
% we usually consider two-sided test to assess the probabilities of
% "at least" and "at most"
% =========================================================================
% If δ1 is very small, the forecast rate is too low (an underprediction);
% If δ2 is very small, the forecast rate is too high (an overprediction);
% =========================================================================
%                            main function
% =========================================================================
function [delta1,delta2]=N_test(forecast, observe)
if length(forecast)==1
    warning('Usually this function will analyse the quantile of <forecast>\n ');
else
    fprintf('Total %d Samples \n',length(forecast));
end
[delta1, delta2]=cal_Quantiles(forecast, observe);
end
% =========================================================================
%                             Sub-function
% =========================================================================
function [delta1, delta2]=cal_Quantiles(foreCounts, obsCounts)
[xs, cdf]=ECDF(foreCounts);
% δ1: P(x ≥ obs_count)
delta1=Overobs(xs, cdf, obsCounts);
% δ2: P(x ≤ obs_count)
delta2=Lessobs(xs, cdf, obsCounts);
end

function [xs, cdf]=ECDF(x)
xs = sort(x(:));
n = numel(xs);
cdf = (1:n)' / n; % cal cdf 
end

function pover=Overobs(xs, cdf, val)
if val > xs(end)
    pover = 0.0;
elseif val < xs(1)
    pover = 1.0;
else
    idx = find(xs >= val-1, 1, 'first');
    pover = 1-cdf(idx);
end
end

function prob = Lessobs(xs, cdf, val)
if val > xs(end)
    prob = 1.0;
elseif val < xs(1)
    prob = 0.0;
else
    idx = find(xs <= val, 1, 'last');
    prob = cdf(idx);
end
end
