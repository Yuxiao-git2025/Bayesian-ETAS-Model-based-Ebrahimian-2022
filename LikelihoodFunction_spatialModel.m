%                       (New) Likelihood Function
function [lik,loglik] = LikelihoodFunction_spatialModel (data, theta)
%% Read parameters
M      = data{1,1}; T      = data{1,2};
tstart = data{1,3}; Ml     = data{1,4};
r      = data{1,5}; rxy    = data{1,6};
dA     = data{1,7}; cal_K  = data{1,8};
Appro  = data{1,9}; Nbo    = data{1,10};
mu_xiyi_Ml= data{1,11};

%% Estimate the integral for distance
if Appro
    Ir=[];
else
    Ir=calculate_Ir (rxy, theta(6:end), dA, M);
end

%% Calculate K 
if cal_K
    theta(2)=calculate_Kseq (M, T, Ir, tstart, Ml, theta, Nbo);
end

%% Main Loop
loglam = 0;
intlambda = 0;
% calculate the λ in each time-point, and <log>, then sum them up
for i=1:length(T)
    loglam = loglam+log(lambdaETAS_likelihood (M, T, r, M(i), T(i), Ml, theta, mu_xiyi_Ml(i)));
    if (i>=2 && ~cal_K)
        intlambda = intlambda+int3LambdaETAS (M, T, Ir, T(i-1), T(i), Ml, theta);
    end
end

if cal_K
    lik = exp(loglam);
    loglik = loglam;
else
    intlambda = intlambda+int3LambdaETAS (M, T, Ir, T(end), tstart, Ml, theta);
    lik = exp(loglam)*exp(-intlambda);
    loglik = loglam-intlambda;
end

end
        