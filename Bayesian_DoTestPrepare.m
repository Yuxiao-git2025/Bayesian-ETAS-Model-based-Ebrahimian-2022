function [Nmean,N50,N16,N84,N02,N98,Nobs,Pmean,P98, ...
    Nmeanxy,N50xy,N16xy,N84xy,N02xy,N98xy]=Bayesian_DoTestPrepare()
global vecM Xcgrid Ycgrid
    Nmean= zeros(length(vecM),1); N50  = zeros(length(vecM),1);
    N16  = zeros(length(vecM),1); N84  = zeros(length(vecM),1);
    N02  = zeros(length(vecM),1); N98  = zeros(length(vecM),1);
    Nobs = zeros(length(vecM),1);

    Pmean= zeros(length(vecM),1); P98  = zeros(length(vecM),1);

    Nmeanxy= zeros(length(Ycgrid),length(Xcgrid),length(vecM));
    N50xy  = zeros(length(Ycgrid),length(Xcgrid),length(vecM));
    N16xy  = zeros(length(Ycgrid),length(Xcgrid),length(vecM));
    N84xy  = zeros(length(Ycgrid),length(Xcgrid),length(vecM));
    N02xy  = zeros(length(Ycgrid),length(Xcgrid),length(vecM));
    N98xy  = zeros(length(Ycgrid),length(Xcgrid),length(vecM));
end