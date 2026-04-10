%                         Block-Wise MCMC procedure

function [THETA,accept] = sample_posterior_MCMC_updated...
(THETA,seeds,weights,priorPDF,priorPDF_par,weights_prior,likelihoodFunction,data)
%% Extra Calculation (just is needed for ETAS parameter updating)

cal_K = data{1,8};
if cal_K == 1
    seeds(2,:) = [];
    tempK = THETA(2);
    THETA(2) = [];
    priorPDF(2) = [];
    priorPDF_par(2) = [];
    indx_gr_1 = [5,7]-1;
else
    indx_gr_1 = [5,7];
end

%% sampling θ from proposal PDF 
% Refer to SI-2
new_THETA = sampleTheta_kernel(seeds,weights);
% make sure positive
if (any(new_THETA<0) || any(new_THETA(indx_gr_1)<=1.00))
    while (any(new_THETA<0) || any(new_THETA(indx_gr_1)<=1.00)) 
        new_THETA = sampleTheta_kernel(seeds,weights);
    end
end

% After obtaining the parameters, we will go back to calculate the kernel
% function κold and κNew, and then divide them to obtain the proposal 
% distribution ratio
% proposal_ratio = kernelPDF(THETA,seeds,'adaptive')...
% /kernelPDF(new_THETA,seeds,'adaptive');
proposal_ratio = calculateKernel(THETA,seeds,weights)...
    /calculateKernel(new_THETA,seeds,weights);

%% calculate the ratio of priors

if ~strcmp(priorPDF{1,1},'kernel')
    if ~all(strcmp(priorPDF,'normal') | strcmp(priorPDF,'lognormal'))
        ratioPrior = zeros(1,length(THETA));
        for i=1:length(THETA)

            if strcmp(priorPDF{1,i},'normal')
            meanpriorPDF   = priorPDF_par{1,i}(1);
            sigmapriorPDF  = priorPDF_par{1,i}(2);
            ratioPrior(i)  = normpdf(new_THETA(i), meanpriorPDF, sigmapriorPDF)...
                /normpdf(THETA(i), meanpriorPDF, sigmapriorPDF);

            elseif strcmp(priorPDF{1,i},'lognormal')
            medianpriorPDF = priorPDF_par{1,i}(1);
            betapriorPDF   = priorPDF_par{1,i}(2);
            ratioPrior(i)  = lognpdf(new_THETA(i), log(medianpriorPDF), betapriorPDF)...
                /lognpdf(THETA(i), log(medianpriorPDF), betapriorPDF);

            elseif strcmp(priorPDF{1,i},'uniform')
            ratioPrior(i)  = 1.0;
            end

        end
        ratioPrior=prod(ratioPrior);
    else
        priorNormal_par = cell2mat(priorPDF_par');
        ratioPrior = calculateMVN(new_THETA,priorNormal_par(:,1),priorNormal_par(:,2) ...
            ,eye(length(new_THETA)),priorPDF{1,1})/calculateMVN(THETA,priorNormal_par(:,1) ...
            ,priorNormal_par(:,2),eye(length(new_THETA)),priorPDF{1,1});
    end
else
    samplePrior = cell2mat(priorPDF_par);    
    ratioPrior  = calculateKernel(new_THETA,samplePrior,weights_prior)...
        /calculateKernel(THETA,samplePrior,weights_prior);
end

%% calculate the likelihood function 
if cal_K
    new_THETA = [new_THETA(1); tempK; new_THETA(2:end)];
    THETA = [THETA(1); tempK; THETA(2:end)];
end
[Likelihoodnew,logLikelihoodnew] = likelihoodFunction(data,new_THETA);
[Likelihood,logLikelihood]       = likelihoodFunction(data,THETA);

%% calculate likelihood ratio
if (isnan(Likelihoodnew/Likelihood) || Likelihoodnew/Likelihood==0)
    pratio =  exp(logLikelihoodnew-logLikelihood);
else
    pratio =  (Likelihoodnew/Likelihood);
end
%% calculate the acceptance probability
r=pratio*proposal_ratio*ratioPrior;
alpha=min([1 r]); 
u=rand();                                    
if u <= alpha
    THETA=new_THETA; 
    accept=1;
else
    accept=0;
end
        
end

        