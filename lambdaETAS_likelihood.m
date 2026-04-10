% =========================================================================
% Calculation Likelihood Function of events with M=m at t=time 
% λ(t,x,y,M=m|θ,seq,Ml)=β*exp(-β(m-Ml))*λ(t,x,y|θ,seq,Ml)
% =========================================================================
% Notice there is no int term in lambda expression, because when we choose
% <cal_K>, the integral is equal to No, where No is number of events in the 
% area(No=length(indexSeq))
% Thus when employing the Bayesian updating framework, exp(-int(λ)) is 
% eliminated when divided
% =========================================================================
function lambda = lambdaETAS_likelihood (Mi, Ti, ri, m, time, Ml, theta, mu_xiyi_Ml)
% Ti is time sequence and time is current time0
beta  = theta(1);
K     = theta(2);
alpha = theta(3);
c     = theta(4);
p     = theta(5);
d     = theta(6);
q     = theta(7);

index = find(Ti < time);  

if length(theta)==8
    gamma = theta(8);
    d = d*exp(gamma*Mi(index));
end

if ~isempty(index)
    Krt = (q-1)/pi*d.^(2*(q-1))*(p-1)*c^(p-1); % normalized factor
    lambda = beta*exp(-beta*(m-Ml))*(mu_xiyi_Ml+sum((K*exp(alpha*(Mi(index)-Ml)) ...
        ./((time-Ti(index)+c).^p)).*(Krt./(ri{index(end),1}.^2+d.^2).^q)));
else
    lambda = beta*exp(-beta*(m-Ml));
end
    
end

        