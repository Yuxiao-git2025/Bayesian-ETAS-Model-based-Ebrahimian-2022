function [samples,indexSeq,Muxiyi]=Bayesian_DoMCMC(DoMCMC,output_Dir)
global time_MS time M T_start T_end Mc r ...
    latitude longitude Xcgrid Ycgrid Muxy
% Define the seq (i.e history observation)
indexSeq = find(time_MS >= 0 & time < T_start & M >= Mc);
Muxiyi = zeros(1,length(indexSeq));
samples=[];
if DoMCMC == 1
    fprintf('# Use MCMC to find the parameters of ETAS model \n');
    fprintf('# Forecast timespan is [%.4f  %.4f] \n',T_start,T_end);
%     fprintf(['History observations:'  '\n' ]);
%     fprintf('%d \n',indexSeq);
    %                           Calculate r
    % The cell array r{j-1} stores the spherical projection distances from the
    % jth event to the previous j-1 events.   ​
    r = cell(length(indexSeq)-1,1);
    for j=2:length(indexSeq)
        Grdata  = [latitude(indexSeq(j)),longitude(indexSeq(j))]*pi/180;          
        % Read the aftershock latitude and longitude, gradi a radian
        Grpdata = [latitude(indexSeq(1:j-1)),longitude(indexSeq(1:j-1))]*pi/180; 
        % Read the aftershock latitude and longitude, gradi a radian
        xv = topgeo(Grpdata(:,1),Grpdata(:,2),Grdata(1),Grdata(2));
        r{j-1,1} = sqrt(xv(:,1).^2+xv(:,2).^2);
    end

    %                       Find background Grids
    % Match the longitude and latitude of each event to the nearest grid 
    % point and extract the preset background activity rate Muxy.
    for j=1:length(indexSeq)
        [~,indx_lat] = min(abs(Ycgrid-latitude(indexSeq(j))));
        [~,indx_lon] = min(abs(Xcgrid-longitude(indexSeq(j))));
        Muxiyi(j) = Muxy(indx_lat,indx_lon);
    end
    %                       Main Loop (for MCMC)
    % =====================================================================
    samples=sample_MainMCMC(indexSeq,Muxiyi);
    % =====================================================================

    %                           save
    filename = fullfile(output_Dir,'\WY_samples.mat');
    save(filename,'samples','indexSeq','-v7.3');

else
    load(fullfile(output_Dir,'\samples.mat'));
    return
end