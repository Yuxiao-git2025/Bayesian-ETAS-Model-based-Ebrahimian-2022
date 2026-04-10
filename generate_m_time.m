% Generate mag and time
%                     (using Thinning method)
function [Ms,ts,Num] = generate_m_time (Mi, Ti, Iri, tstart, tend, ...
    Ml, theta, Mmax, Mu_b)

% generate mag
Ms=generate_M(Ml, theta(1), Mmax);

Maxlam=int2lambdaETAS (Mi,Ti,Iri,tstart,Ml,theta,Mu_b);
accept=0;
Num=0;
% tstart vary different and keep going point-by-point from <tstart> to <tend>
ts=tstart;
while (accept == 0 && ts <= tend)
    Num=Num+1;
    % generate time
    ts=generate_time(ts, Maxlam);
    lamNew=int2lambdaETAS(Mi,Ti,Iri,ts,Ml,theta,Mu_b);
    ratio=lamNew/Maxlam;
    u=rand();
    if (u<=ratio && ratio<=1)
        accept=1; % out of loop
    else
        % λmax is going to change while generating each new ti
        Maxlam=lamNew;
    end 
end
end
