function [del1, del2] = N_Poitest(Nfore, Nobs, alpha)
% N-test based on poisson numbers test
nf=length(Nfore);
if nargin<3 || isempty(alpha)
    alpha=0.05;
end
del2 = poisscdf(Nobs, Nfore);
del1 = 1-poisscdf(Nobs-1, Nfore);
alphaeff=alpha/2;
% consistent = ~(del1 < alphaeff || del2 < alphaeff);
fprintf('α= %.4f \n',alpha);
fprintf('δ1= %.2f, δ2= %.2f  \n',del1,del2);
end