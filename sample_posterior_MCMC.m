% =========================================================================
%                        Component-Wise MCMC procedure
% =========================================================================
% The proposal distribution for each component is assumed to be a lognormal 
% distribution; The prior PDF is a lognormal distribution too;
% cal_K should be done after each parameter updated;
% so r=p(likelihood ratio)*p(prior ratio)*q(proposal ratio)
% =========================================================================
function [THETA,accept] = sample_posterior_MCMC...
(THETA,rank,proposalPDF,proposalPDF_par,priorPDF,priorPDF_par,likelihoodFunction,data)
% current θ (6-8 length)
theta = THETA(rank);

%% calculate proposal ratio
% Sampling new_theta from proposa PDF q(θ) or calculate Proposal Ratio
% q(θ)/q(θ*)
while 1
if strcmp(proposalPDF,'normal')
    
    sigma_theta    = proposalPDF_par(1);
    new_theta      = normrnd(theta,sigma_theta);  
    if (rank==5 || rank==7)   % p or q
        while  new_theta<=1.00
            new_theta  = normrnd(theta,sigma_theta); 
        end
    end
    proposal_ratio = 1.0;
    
elseif strcmp(proposalPDF,'lognormal')
    
    beta_theta     = proposalPDF_par(1);
    new_theta      = lognrnd(log(theta),beta_theta);  
    if (rank==5 || rank==7)   % p or q should >1
        while  new_theta<=1.00
            new_theta  = lognrnd(log(theta),beta_theta);  
        end
    end
    proposal_ratio = lognpdf(theta, log(new_theta), beta_theta)/lognpdf(new_theta, log(theta), beta_theta);
    
elseif strcmp(proposalPDF,'uniform')
    
    thetamin  = proposalPDF_par(1);
    thetamax  = proposalPDF_par(2);
    new_theta = unifrnd(thetamin,thetamax); 
    if (rank==5 || rank==7)   % p or q
        while  new_theta<=1.00
            new_theta  = unifrnd(thetamin,thetamax); 
        end
    end
    proposal_ratio = 1.0;
    
elseif strcmp(proposalPDF,'kernel')   % more information refer to supplement  
    
    seeds     = proposalPDF_par{1,1}(rank,:);
    weights   = proposalPDF_par{1,2};
    new_theta = sampleTheta_kernel(seeds,weights);
    if (rank==5 || rank==7)   % p or q
        while  new_theta<=1.00
            new_theta = sampleTheta_kernel(seeds,weights);
        end
    end    
    proposal_ratio = calculateKernel(theta,seeds,weights)/calculateKernel(new_theta,seeds,weights);

end
% make sure positive
if new_theta > 0
    break
end

end

%% calculate priors ratio

if strcmp(priorPDF,'normal')
    meanpriorPDF   = priorPDF_par(1);
    sigmapriorPDF  = priorPDF_par(2);
    ratioPrior     = normpdf(new_theta, meanpriorPDF, sigmapriorPDF)/normpdf(theta, meanpriorPDF, sigmapriorPDF);
    
elseif strcmp(priorPDF,'lognormal')
    medianpriorPDF = priorPDF_par(1);
    betapriorPDF   = priorPDF_par(2);
    ratioPrior     = lognpdf(new_theta, log(medianpriorPDF), betapriorPDF)...
    /lognpdf(theta, log(medianpriorPDF), betapriorPDF); 
    
elseif strcmp(priorPDF,'uniform')
    ratioPrior     = 1.0;
end

%% calculate the likelihood function 
% update θ one-by-one
THETA(rank) = new_theta;
[Likelihoodnew,logLikelihoodnew]=likelihoodFunction(data,THETA);
THETA(rank) = theta;
[Likelihood,logLikelihood]=likelihoodFunction(data,THETA);

%% calculate likelihood ratio
if (isnan(Likelihoodnew/Likelihood) || Likelihoodnew/Likelihood==0)
    pratio =  exp(logLikelihoodnew-logLikelihood);
else
    pratio =  (Likelihoodnew/Likelihood);
end

%% calculate the acceptance probability
r=pratio*ratioPrior*proposal_ratio;
alpha=min([1 r]); 
u=rand();                                    
if u <= alpha
    theta=new_theta; 
    accept=1;
else
    accept=0;
end
        
THETA(rank)=theta;

end        