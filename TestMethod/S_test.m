% S-Test based on Likelihood-based tests for evaluating earthquake
% forecasts.
% Method1: Standard
% Method1: Simulation-based
function [prob1,prob2] = S_test(Nobsxy,Nforexy,Nsimxy)
% Eliminate the prediction of the total number <Nfore> compared with the
% total number of observed earthquakes <Nobs> 
% This summing and normalization procedure removes the effect of the rate 
% and magnitude components of the original forecast 
Nforexy = Nforexy*sum(Nobsxy)/sum(Nforexy); % equation(20) in Zechar.2010

% The forecast in each bin follows a Poisson distribution:Pr(w|¦Ë)=¦Ë^w/w!*exp(-¦Ë)
% Joint likelihood (for all bins) of the observation given the forecast:
% Pr1*Pr2*...Prn, then we log it(logPr(w|¦Ë)), and sum log-likelihood in each bins
Likelihood_obs = sum(-Nforexy+Nobsxy.*log(Nforexy)-log(factorial(Nobsxy)));

Nfore = 1:length(Nforexy);

%% Method 1 (Zechar et al. 2010)

Sims = 1000;
lik1 = zeros(Sims,1);
Fx = cumsum(Nforexy/sum(Nforexy)); % CDF
Nobs = sum(Nobsxy);
for i=1:Sims
    Nxy_sim = zeros(length(Nforexy),1);
    for j=1:Nobs
        sampleX=Nfore(find(Fx > rand, 1, 'first'));
        Nxy_sim(sampleX)=Nxy_sim(sampleX)+1; 
    end
    lik1(i) = sum(-Nforexy+Nxy_sim.*log(Nforexy)-log(factorial(Nxy_sim)));
end    

%% Method new 
lik2 = zeros(size(Nsimxy,2),1);
for i=1:size(Nsimxy,2)
    Nxy_sim = zeros(length(Nforexy),1);
    Nsimxy(:,i) = Nsimxy(:,i)*sum(Nobsxy)/sum(Nsimxy(:,i));
    Fx = cumsum(Nsimxy(:,i)/sum(Nsimxy(:,i)));
    for j=1:Nobs
        sampleX = Nfore(find(Fx>rand,1,'first'));
        Nxy_sim(sampleX) = Nxy_sim(sampleX)+1; 
    end
    lik2(i) = sum(-Nforexy+Nxy_sim.*log(Nforexy)-log(factorial(Nxy_sim)));

end

%% Calculate the probabilities
% equation(23) in Zechar.2010
prob1 = length(find(lik1<=Likelihood_obs))/length(lik1);
prob2 = length(find(lik2<=Likelihood_obs))/length(lik2);
end

