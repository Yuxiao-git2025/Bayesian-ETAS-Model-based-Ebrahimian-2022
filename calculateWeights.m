% =========================================================================
% Function for calculating adaptive/non-adaptive weights in Kernel Density
%                      (Refercence: Beck.2002)
% =========================================================================
% Method:Construct a kernel density function <κ> as the weighted sum (average) 
% of n-dimensional Gaussian PDFs centered among these samples;
% =========================================================================
function weights = calculateWeights(seeds, type)

[n, nSeeds] = size(seeds);

%% Calculate width (fixed-width) and Md (the number of distinct samples)

% temp = seeds;
% nCopies = 0;
% while size(temp,2)>1
%     temp2 = [];
%     for i = 2:size(temp,2)
%         if isequal(temp(:,i),temp(:,1))
%             nCopies = nCopies+1;
%         else
%             temp2 = [temp2 temp(:,i)];
%         end
%     end
%     temp = temp2;
%     clear temp2;
% end

nCopies = size((unique(seeds','rows'))',2);
if nSeeds==nCopies
    nCopies = 0;
end
% Md is the number of distinct samples (≤Nseed) 
Md = nSeeds-nCopies;
width = (4/((n+2)*Md))^(1/(n+4));

%% Calculate the local bandwidth factors, lambda (New Method)
% The kernel function <κ> can be viewed as a PDF consisting of bumps at θi,
% Therefore, a large value of <wi> tends to over-smooth the kernel density,
% while a small value may cause noise-shaped bumps; In the adaptive kernel
% method, the idea is to use a larger width in regions of lower probability
% density, so wi=width*λi, where λi is expressed:
if strcmp(type,'adaptive')
    % <Abramson.1982> showed that varying the bandwidth proportional to
    % κ^(-1/2) (i.e., omiga=0.50) provides good performance for the kernel
    % smoothing density.
    omiga = 0.50;   % also [1/n], 0<=omiga<=1
    kp = zeros(nSeeds,1);
    S = cov(seeds');
    disp('# Calculate weights for the samples') 

    for i = 1:nSeeds
%         disp([' Calculate weight for sample ',num2str(i),'/',num2str(nSeeds)]) 
        X = diag((repmat(seeds(:,i),1,nSeeds)-seeds)'/S*(repmat(seeds(:,i),1,nSeeds)-seeds));
        kp(i) = mean(exp(-0.5*X/width^2)/width^n*1/(sqrt(det(S)*(2*pi)^n)));
    end
    % Multikp could be considered as normalized factor
    Multikp = prod(kp.^(1/nSeeds));
    lambda = (Multikp./kp).^omiga; 
else
    lambda = ones(nSeeds,1);
end
% Adaptive weights is :
weights = width*lambda;

return
end
%%  Calculate the local bandwidth factors, lambda (OLD)

% if strcmp(type,'adaptive')
%     alpha = 0.50;   % also [1/n]
%     kp = zeros(nSeeds,1);
%     for i = 1:nSeeds
%         kp(i) = calculateKernel(seeds(:,i),seeds,width*ones(nSeeds,1));
%     end
%     prodkp = prod(kp.^(1/nSeeds));
%     lambda = (prodkp./kp).^alpha;
% else
%     lambda = ones(nSeeds,1);
% end