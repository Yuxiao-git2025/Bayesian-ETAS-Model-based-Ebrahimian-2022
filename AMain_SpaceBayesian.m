%% ========================================================================
% Seismicity forecasting based on a Bayesian spatio-temporal ETAS model
% written by: Hossein Ebrahimian and Fatemeh Jalayer  
%                       Modified by Yu 2025/9/9
clc;
clear;
filename = 'catalog.txt'; 
OutDir='output';
if exist(OutDir,'dir')==0
    mkdir(OutDir);
end
%% ========================================================================
[Method,UseBKGD,UseGamma]=Bayesian_PrepareData(filename);


%%                      <1> MCMC algorithm 
clc;
DoMCMC=1;   
[samples,indexSeq,Muxiyi]=Bayesian_DoMCMC(DoMCMC,OutDir);


%%                      <2> Finding the Robust Estimate 
DoRobust=1;
[Expsamples]=Bayesian_DoRobust(DoRobust,samples,indexSeq,OutDir);

%%                      <3> N-test & S-test
DoNStest=1;
NSsamples=load(fullfile("Savedata",'\Sum1_205.mat'));
[Pmean,Nmean,Nobs,Gridfore,Gridobs,Prob1,Prob2]=...
    Bayesian_DoTest(DoRobust,DoNStest,OutDir,NSsamples.Sum1_205);

% =========================================================================
%%                       Mapping Location 
Bayesian_LocationMapping([],[],'rainbow',5);


%%                      Mapping samples-correlation
Bayesian_CorrMapping(samples);

%%                      Mapping Robust-estimate Grids
Beyesian_GridMapping(Gridfore,'rainbow')


%% At last
% If you use this code, please cite the fillowing two articles:

% (1)Ebrahimian, H., Jalayer, F., Maleki Asayesh, B., Hainzl, S., Zafarani,
% H. Improvements to seismicity forecasting based on a Bayesian
% spatio-temporal ETAS model. Scientific Reports, (2022)
% https://doi.org/10.1038/s41598-022-24080-1.

% (2) Ebrahimian, H. & Jalayer, F. Robust seismicity forecasting based on
% Bayesian parameter estimation for epidemiological spatio-temporal
% aftershock clustering models. Scientific Reports 7, 9803 (2017),
% https://doi.org/10.1038/s41598-017-09962-z.
% -----------------------------------------------------

% Note: This code is derived based on the definitions that are presented in
% (1) above. Therefore, for some input parameters, it is required to see
% the paper accordingly.

        