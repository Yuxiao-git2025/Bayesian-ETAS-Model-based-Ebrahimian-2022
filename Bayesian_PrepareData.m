function [Method,UseBKGD,UseGamma]=Bayesian_PrepareData(filename)
%                  ||  Important Input Parameters  ||
global T_start T_end deltaGrid lonMin lonMax latMin latMax ...
    latitude longitude M time time_MS tstart tend ...
    Ggrid rxy dA Xcgrid Ycgrid Ngrid ...
    Cal_K Appro Mc Mmax vecM Muxy Muxy2col Nbo ...
    PriPDF theta
% Forecasting interval
T_start = 11.6;
T_end   = 12.0;
% T_start = 11.6250;
% T_end   = 12.00;

% Aftershock Zone
% The outregion will not be considered when predict events
% Adjust dGrid so that better fitted like 200X200 Grids
deltaGrid = 0.01;
lonMin = 45;
lonMax = 47;
latMin = 32.50;
latMax = 35.50;

% Read catalog
[latitude,longitude,M,time,time_MS,tstart,tend]=Bayesian_ReadCatalog(filename);

% Set Grids
[Ggrid,rxy,dA,Xcgrid,Ycgrid,Ngrid]=Bayesian_SetGrid();

% Method:Fast  ||  Semi-Fast  ||  Slow
Method='Fast';
[Cal_K,Appro]=Bayesian_SetMethod(Method);

% Magnitude
Mc=3.4;     
Mmax=7.5;    % Any if you want but here empirically it is 7.5
% vecM=[Mc,4,5,6,7];  
% suggest that using a single vecM because calculation is large
vecM=[Mc];  

% Is use Background 
UseBKGD=0;    
% Read BKGD Data
[Muxy,Muxy2col,Nbo]=Bayesian_ReadBKGD(Xcgrid,Ycgrid,T_start,Ngrid,dA,UseBKGD);

% Is ues gamma in ETAS-Model
UseGamma=1;    

% Initial params 
[PriPDF,theta]=Bayesian_SetParams(UseGamma);



end