% This scripts uses MCMC simulation with Metropolis-Hastings procedure to
% estimate the probability P(Data|Model)
% input: (indexSeq,M,time_MS,tstart,Mc,r,rxy,dA,Cal_K,Appro,Nbo,Muxiyi,PriPDF,theta)

function samples=sample_MainMCMC(indexSeq,Muxiyi)
global  M time_MS tstart Mc r rxy dA ...
        Cal_K Appro Nbo  PriPDF theta


% Initialize the Metropolis-Hastings sampler
% Pre-burning period: Eliminate the initial transient effect of the Markov chain
% Multi-chain strategy: Improve sampling efficiency and avoid local optimality
Maxchain1 = 520;   
Maxchain2 = 1000;  
burnin    = 20;               
numChain  = 6;                % number of chain
numUP     = length(theta);    % number of Uncertain Parameters (7 or 8)
kernelType= 'adaptive';       % kernelType = 'adaptive' / 'nonadaptive'
nbins = 30;                  % number of bins for the histogram

%% Define the PDF for Prior and Proposal Distributions 
% Prior: Using log-normal distribution
priorPDF = {};
priorPDFfunction = {};
for i=1:numUP 
    priorPDF = [priorPDF,'lognormal'];
    % x is value and par is [mu,sigma]
    priorPDFfunction = [priorPDFfunction,{@(x,par) lognpdf(x,log(par(1)),par(2))}];  
end
              
% Proposal Distributions
proposalPDF = {};
proposalPDF_par = {};
for i=1:numUP 
    proposalPDF = [proposalPDF,'lognormal'];
    proposalPDF_par = [proposalPDF_par, 0.30];    % COV
end
weights_prior = [];

%% likelihoodFunction definition
likelihoodFunction = @LikelihoodFunction_spatialModel;
DATA = {M(indexSeq),time_MS(indexSeq),tstart,Mc,r,...
    rxy(indexSeq,:),dA,Cal_K,Appro,Nbo,Muxiyi};
fprintf('Totally Summation %d events before <Tstart> \n\n', length(M(indexSeq)));
%%                              Start sampling
% Using Adaptive MH algorithm (introduces a sequence of intermediate
% candidate evolutionary PDFs that resemble more and more the target PDF)
nchain=1;
while nchain <= numChain
    waithandel=waitbar(0,sprintf ...
        (['Loading process in ' num2str(nchain) ' chains \n']));
    fprintf('======================= \n');
    fprintf('||                   || \n');
    fprintf(['# Chain Number = ',num2str(nchain) '\n']);
    fprintf('||                   || \n');
    fprintf('======================= \n\n');

    if nchain==1
        NumberIter = Maxchain1;
    else
        NumberIter = Maxchain2;
    end
    
    state  = zeros(numUP,NumberIter);  % Storage space for our samples
    accept = zeros(numUP,NumberIter);  % Storage space for accept decisions
    fprintf('Please waitting... \n')

    for iter = 1:NumberIter
        waitbar(iter/NumberIter,waithandel);
        if mod(iter,100)==0 
            fprintf(['>> Iter = ',num2str(nchain),'  ', ...
                num2str(iter),'/',num2str(NumberIter) '\n'] )
        end
        % =================================================================
        %                       In the first chain
        if nchain==1       %(Metropolis-Hasting Method)
            for n = 1:numUP
                if (~Cal_K || (Cal_K && n~=2)) % calculate K value solely
                    % Here instead of using <block-wise> updating scheme, 
                    % adapt to <component-wise> updating scheme; cause 
                    % <block-wise> scheme has low accept raito 
                    [theta,accept(n,iter)] = sample_posterior_MCMC...
                        (theta,n,proposalPDF{1,n},proposalPDF_par{1,n}, ...
                        priorPDF{1,n},PriPDF{1,n},likelihoodFunction,DATA);
                end
            end
            
        else
            % =============================================================
            %                       In the other chains
            [theta,accept(:,iter)] = sample_posterior_MCMC_updated...
                (theta,seeds,weights,priorPDF,PriPDF,weights_prior, ...
                likelihoodFunction,DATA);
        
        end
        % After all parameters are updated, then calculate K value(param(2))
        if Cal_K
            if Appro
                Ir = [];
            else
                Ir = calculate_Ir (rxy(indexSeq,:), theta(6:end), dA, M(indexSeq));
            end
            theta(2) = calculate_Kseq (M(indexSeq),time_MS(indexSeq), ...
                Ir, tstart, Mc, theta, Nbo);
        end
        % renew the theta into state
        state(:,iter) = theta;      
    end % iter end here! 
    delete(waithandel);
        
    %% Samples we take for further analysis
    % The seeds of each chain serve as the initial sample library for the
    % next chain, while the last(6) chain(1e3 iter) output the final samples
    if nchain==1
        % Burn some samples and record down (component-wise)
        seeds=state(:,burnin+1:NumberIter);
    else
        % Identify all different parameter combinations (block-wise)
        [~,ind] = unique(state','rows');
        seeds = state(:,sort(ind));
    end

    if Cal_K
        weights = calculateWeights(seeds([1,3:end],:), kernelType);
    else
        weights = calculateWeights(seeds, kernelType);
    end    
 
    nchain = nchain + 1;
    
end % while end here !
% In each chains the state and accept will update but the seeds are unique
% finally we storage the seeds into samples
samples = seeds;
end

