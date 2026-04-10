% Function for generating SEQgen
function [Ms,ts,lons,lats,rs,Intlam_seq]=...
    generateSEQ (Mi, Ti, rxy, sampletheta)
global Xcgrid Ycgrid Ggrid Mmax dA Appro Muxy2col tstart tend Mc
fprintf('# Transformed timespan: [%.4f  %.4f] \n',tstart,tend);
Ms  = [];ts = [];
lons= [];lats  = [];
rs  = [];Intlam_seq = [];

Mu_b=sum(Muxy2col)*dA;

%% Generating 1st Magnitude and time
% Estimate the integral for distance
if ~Appro
    Iri=calculate_Ir (rxy, sampletheta(6:end), dA, Mi);
else
    Iri=[];
end
[Mgen,tgen,Num]=...
    generate_m_time(Mi, Ti, Iri, tstart, tend, Mc, sampletheta, Mmax, Mu_b);
if tgen > tend
    disp(['tgen = ',num2str(tgen),' > T_end = ',num2str(tend), ...
        ', NO sequence will be generated!!!'])
else
%% Generating 1st set of lat and lon
    [rgen,Longen,Latgen,Intlam_t]=generate_R ...
        (Mi, Ti, rxy, Mgen, tstart, tgen, Mc, sampletheta, Xcgrid, Ycgrid, Ggrid, Muxy2col);
    count=0;
%% Generating more samples
    % =====================================================================
    while tgen <= tend
        disp(['time: ',num2str(tgen) ' ||  mag: ', ...
            num2str(Mgen),' ||  Thinning = ',num2str(Num)])
        count=count+1;
        Ms(count,:)  = Mgen;
        ts(count,:)  = tgen;
        lons(count,:)= Longen;
        lats(count,:)= Latgen;
        rs(count,:)  = rgen;
        Intlam_seq(count,:) = Intlam_t;

        if ~Appro
            Iri=calculate_Ir([rxy;rs], sampletheta(6:end), dA, [Mi;Ms]);
        end
        % Generate more time, mag, and so on
        [Mgen,tgen,Num] = generate_m_time ([Mi;Ms], [Ti;ts], ...
            Iri, tgen, tend, Mc, sampletheta, Mmax, Mu_b);
        [rgen,Longen,Latgen,Intlam_t] = generate_R ([Mi;Ms], [Ti;ts], ...
            [rxy;rs], Mgen, ts(end), tgen, Mc, sampletheta, Xcgrid, Ycgrid, ...
            Ggrid, Muxy2col);
    end
    % =====================================================================

end

end