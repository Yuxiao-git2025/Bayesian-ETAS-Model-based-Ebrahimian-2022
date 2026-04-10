% ==========================================================================
%             Sampling θ from a Kernel function κ
function sample_theta = sampleTheta_kernel(seeds,weights)
% (cumsum is calculated the sum value in each column in current lines)
% Calculate the cumulative sum of the weights and generate the CDF of the
% discrete probability distribution
F = cumsum(weights./sum(weights));

% Generate a uniformly random number within the range [0, 1), 
% and select the component index through inverse transformation sampling
iKernel = find(F>=rand(),1,'first');
% ==========================================================================
Mu = seeds(:,iKernel);

% Notice: weights=(w*λj)^2, and Sj=(w*λj)^2*S 
% (S is the covariance matrix of the samples)
S = (weights(iKernel))^2 * cov(seeds'); 

% Cholesky decomposition
L = chol(S)';
% L=chol(S,"lower");
% L is lower triangular matrix, then: L*L'=S
% Matlab <chol> function gives the upper triangular matrix L'
% ==========================================================================
% Generate the new sample: θ=μ+L*Z (where Z is standard Guassian)
sz=size(Mu);
Z=randn(sz);
sample_theta = Mu+L*Z;

return
%nSeeds = size(seeds,2);
%F = (1:nSeeds)/nSeeds;
%%% Plot
% stairs(1:size(seeds,2),F,'-','color','k','Linewidth',3)
% xlabel('sample','fontsize',14)
% ylabel('CDF','fontsize',14)
% set(gca,'fontsize',12)
% set(gca,'Ytick',0:0.1:1)
% xlim([0 size(seeds,2)])
% ylim([0 1])
% grid on
end