function [latitude,longitude,M,time,time_MS,tstart,tend]=Bayesian_ReadCatalog(filename)
global T_start T_end
% Read The Catalog
Catalog = load(filename);

latitude  = Catalog(:,1);
longitude = Catalog(:,2);
M         = Catalog(:,3);
% Notice time is real catalog time and time_MS is relative time to initial
% time0
time      = Catalog(:,4);
time_MS   = Catalog(:,5);
% So need to transform forecast time to begin with relative time
tstart = T_start-time(1);
tend   = T_end-time(1); 

end