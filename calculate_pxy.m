% Joint PDF
function [pxy,Intlam_t]=...
    calculate_pxy(Mi, Ti, ri, m, tstart, tend, Ml, theta, mu_xy_Ml)

beta  = theta(1);
K     = theta(2);
alpha = theta(3);
c     = theta(4);
p     = theta(5);
d     = theta(6);
q     = theta(7);
%% Calculate lambda
index = find(Ti < tend);    
if length(theta)==8
    gamma = theta(8);
    d = d*exp(gamma*Mi(index));
end
Krt = (q-1)/pi*d.^(2*(q-1))*(p-1)*c^(p-1);
lambda=(mu_xy_Ml+sum((K*exp(alpha*(Mi(index)-Ml)) ...
    ./((tend-Ti(index)+c).^p)).*(Krt./(ri(index).^2+d.^2).^q)));

%% Calculate intLambda
if p == 1    
    Io = log((tend-Ti(index)+c)./(tstart-Ti(index)+c));     
else    
    Io = ((tend-Ti(index)+c).^(1-p)-(tstart-Ti(index)+c).^(1-p))/(1-p);
end
Intlam_t=sum((K*exp(alpha*(Mi(index)-Ml)).*Io).*(Krt./(ri(index).^2+d.^2).^q));
%% Calculate pxy (equation 31)
pxy=beta*exp(-beta*(m-Ml))*lambda*exp(-(mu_xy_Ml*(tend-tstart)+Intlam_t)) ;

end
     