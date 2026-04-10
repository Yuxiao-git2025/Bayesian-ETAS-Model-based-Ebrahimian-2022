%                      || Main Algorithm ||
% if do_calculate_K=1, Calculate K in a algorithm  
% if do_calculate_K=0, Leave it to MCMC sampling
% if consider_approximate=1, Approximate the spatio-Int  
% if consider_approximate=0, Exact estimating spatio-Int
function [Cal_K,Appro]=Bayesian_SetMethod(Method)
switch Method
    case 'Fast'
        Cal_K = 0;       
        Appro = 1;   

    case 'Semi-Fast' 
        Cal_K = 1;       
        Appro = 1;

    case 'Slow' % default setting 
        Cal_K = 1;        
        Appro = 0; 
    otherwise 
        warning('No method for analysis is defined!')     
end

end