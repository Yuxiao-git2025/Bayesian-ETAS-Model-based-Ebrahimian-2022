%                    ETAS initial parameters
%                     [beta,K,alfa,c,p,d,q]
%       (the order of parameters is important so do not change)
function [PriPDF,theta]=Bayesian_SetParams(use_kernel_distance_m)

beta_ini  = 1.0*log(10);   
K_ini     = 5.0;          
% Notice: if do_calculate_K = 1, you can assign K an arbitrary value
% (because it is not considered)
alpha_ini = 1.0*log(10);   
c_ini     = 10^(-1.53);
p_ini     = 1.10;
d_ini     = 1.00;
q_ini     = 1.50;
gamma_ini = 0.20;


%                     COV (default setting 0.5)
cov_beta  = 0.5;cov_K     = 1.0;  
cov_alpha = 0.5;cov_c     = 0.5;     
cov_p     = 0.5;cov_d     = 0.5;
cov_q     = 0.5;cov_gamma = 0.5;

%                        Vector of ETAS parameters
vec_beta  = 0.01:0.01:4.00;vec_K     = 0.1:0.1:20;
vec_alpha = 0.01:0.01:4.00;vec_c     = 0.001:0.001:0.50; 
vec_p     = 0.01:0.01:4.00;vec_d     = 0.01:0.01:15.00;
vec_q     = 0.01:0.01:5.00;vec_gamma = 0.01:0.01:1.00;

PriPDF = {[beta_ini,cov_beta],...
                [K_ini,cov_K],...
                [alpha_ini,cov_alpha],...
                [c_ini,cov_c],...
                [p_ini,cov_p],...
                [d_ini,cov_d],...
                [q_ini,cov_q]};
    
theta = [beta_ini;K_ini;alpha_ini;c_ini;p_ini;d_ini;q_ini];    

if use_kernel_distance_m == 1
    PriPDF = [PriPDF,[gamma_ini,cov_gamma]];
    theta = [theta;gamma_ini];
end

end