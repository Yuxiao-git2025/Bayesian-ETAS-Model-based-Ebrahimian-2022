% Generate Time
%                 using inverse transform sampling
function ts=generate_time(tstart, Maxlam)
u=rand();
IAT=-1/Maxlam*log(1-u);
ts=tstart + IAT;
end

