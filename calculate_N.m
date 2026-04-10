% Calculation of integral of rate of events with M=m in the interval [tstart,tend]
% Mi & Ti  : the vectors of events
% m : mag sequence greater than m is what we calculated
function N = calculate_N (Mi, Ti, ri, m, tstart, tend, Mc, theta, ...
    Intlam_t, dA, mu_xy_Ml)
beta = theta(1);   K= theta(2);
alpha= theta(3);   c= theta(4);
p    = theta(5);   d= theta(6);
q    = theta(7);

index = find(Ti < tend);    
if length(theta)==8
    gamma=theta(8);
    d=d*exp(gamma*Mi(index));
end
Krt=(q-1)/pi*d.^(2*(q-1))*(p-1)*c^(p-1);

if ~isempty(Intlam_t) 
    if p == 1    
        Intcp = log((tend-Ti(index)+c)./(Ti(end)-Ti(index)+c));     
    else
        % The integral of [tstart tend] but we also consider the timespan
        % in [Ti(end) tstart] so exactly the range is [ti(end) tend]
        Intcp = ((tend-Ti(index)+c).^(1-p)-(Ti(end)-Ti(index)+c).^(1-p))/(1-p);
    end
    % calculate the integral of λ(forecast seqg)
    Intlam_fore = sum((K*exp(alpha*(Mi(index)-Mc))).*(Intcp).*(Krt./(ri(index).^2+d.^2).^q));
    % Already known Int λ in [tstart tend] note as Intlam_t 
    % Intlam_t is calculated by <calculate_pxy> function, and the form is
    % similar to Intlam_fore
    % Int E[N] dxdy is calculated as:(Refer Zhuang.2011)
    N=exp(-beta*(m-Mc)) * (mu_xy_Ml*(tend-tstart) + sum(Intlam_t) + Intlam_fore);
    N=N*dA;


else
    if p == 1
        Intcp = log((tend-Ti(index)+c)./(tstart-Ti(index)+c));     
    else
        Intcp = ((tend-Ti(index)+c).^(1-p)-(tstart-Ti(index)+c).^(1-p))/(1-p);
    end
    N=exp(-beta*(m-Mc))*(mu_xy_Ml*(tend-tstart)+sum((K*exp(alpha*(Mi(index)-Mc))) ...
        .*(Intcp).*(Krt./(ri(index).^2+d.^2).^q)));
    N=N*dA;
end

end
       