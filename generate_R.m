% Function for calculating the spatial distribution of the generated events
function [rgen,Longen,Latgen,Intlam_t] = generate_R...
    (Mi, Ti, ri, Mgen, tstart, tgen, Ml, theta, Xgrid, Ygrid, Ggrid, mu_xy_Ml)
pxy=zeros(length(Ggrid),1);
Intlam_t = zeros(1,length(Ggrid));

for j=1:length(Ggrid)
    % specify the joint PDF P(x,y)
    [pxy(j),Intlam_t(j)] = calculate_pxy...
        (Mi, Ti, ri(:,j), Mgen, tstart, tgen, Ml, theta, mu_xy_Ml(j));
end
% normalized
pxy=(reshape(pxy,length(Xgrid),length(Ygrid)))' / sum(pxy); 
% sampling the lon and lat by inverse method
[Longen,Latgen]=sample_pxy(pxy,Xgrid,Ygrid);
% calculate distance in each point to Grids
rgen=calculate_rxy(Latgen,Longen,Ggrid);

end