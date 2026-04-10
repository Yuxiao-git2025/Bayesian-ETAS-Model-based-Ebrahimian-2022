function [Pmean,Nmean,Nobs,Gridfore,Gridobs,Prob1,Prob2]=Bayesian_DoTest...
    (DoRobust,DoNStest,Output,sampleN)
global time T_start T_end  M  Mc rxy Xcgrid Ycgrid  Ggrid  vecM
if DoNStest == 1
    if ~DoRobust || DoRobust==0
        filename=fullfile(Output,'\robust estimate N.mat');
        load(filename)
    end

    %% N-test
    % Calculate the Total Number of Events
    [Nmean,N50,N16,N84,N02,N98,Nobs,Pmean,~, ...
    Nmeanxy,N50xy,N16xy,N84xy,N02xy,N98xy]=Bayesian_DoTestPrepare();
    % =====================================================================
    for k=1:length(vecM)
        disp(['# Calculate N-test with M >= ',num2str(vecM(k))])
        % Sum up total samples and obtained quantiles in each cutoff-magnitude
        [Nmean(k),N50(k),N16(k),N84(k),N02(k),N98(k)]=...
            ordered_statistic(sum(sampleN(:,:,k)));
        % Transform to Grid
        [NmeanGrid,N50Grid,N16Grid,N84Grid,N02Grid,N98Grid]=...
            ordered_statistic(sampleN(:,:,k));
        if k==1
            Nforexy = NmeanGrid; % Choose the first M=Mc using the S-test
            fprintf('Grid Mapping with M > %.1f \n',vecM(k));
        end
        Nmeanxy(:,:,k)= (reshape(NmeanGrid,length(Xcgrid),length(Ycgrid)))';
        N50xy(:,:,k)  = (reshape(N50Grid ,length(Xcgrid),length(Ycgrid)))';
        N16xy(:,:,k)  = (reshape(N16Grid ,length(Xcgrid),length(Ycgrid)))';
        N84xy(:,:,k)  = (reshape(N84Grid ,length(Xcgrid),length(Ycgrid)))';
        N02xy(:,:,k)  = (reshape(N02Grid ,length(Xcgrid),length(Ycgrid)))';
        N98xy(:,:,k)  = (reshape(N98Grid ,length(Xcgrid),length(Ycgrid)))';
        % Observed data
        Nobs(k) = length(find(time>=T_start & time<T_end & M>=vecM(k)));
        % Calculate the probability that: P(M>Ml)
        Pmean(k) = 1-exp(-Nmean(k));
        % Calculate the N-test under Poisson distribution (Classic Method)
        mapsample=sum(sampleN(:,:,k));
        [d1_poi(k),d2_poi(k)]=N_test(mapsample, Nobs(k));
        % Calculate the N-test under Bayesian work-flow (Bayesian Method)
        [d1_bay(k)]=length(find(mapsample>=Nobs(k)))/length(mapsample);
        [d2_bay(k)]=length(find(mapsample<=Nobs(k)))/length(mapsample);

        % Mapping
        figure;
        h(k)=histogram(mapsample,'Normalization','pdf','LineWidth',1.5,'FaceColor', ...
            [.8 .8 .8],'FaceAlpha',0.6,'EdgeColor','k','EdgeAlpha',0.8);
%         h(k).NumBins=8;
        hold on;
        xline(N50(k),'LineWidth',2,'Color','k','LineStyle','--');
        xline(N16(k),'LineWidth',2,'Color','b','LineStyle','--');
        xline(N84(k),'LineWidth',2,'Color','b','LineStyle','--');
        xline(N02(k),'LineWidth',2,'Color','r','LineStyle','--');
        xline(N98(k),'LineWidth',2,'Color','r','LineStyle','--');
        scatter(Nobs(k),0,200,'Marker','hexagram','MarkerFaceColor',[0 1 0],'MarkerEdgeColor', ...
            'k','LineWidth',1);
        scatter(Nmean(k),0,200,'Marker','square','MarkerFaceColor','r','MarkerEdgeColor', ...
            'k','LineWidth',1);
        title(['Model Forecast for M>' num2str(vecM(k)) ]);
        Fun_defaultAxes;
        % Message
        maxs=max(mapsample); 
        mins=min(mapsample);
        fprintf('obs\t mean\t min\t max\t 50CI\t \n');
        fprintf('%.0f\t  %.0f\t  %.0f\t  %.0f\t  %.0f\t \n',Nobs(k),Nmean(k),mins,maxs,N50(k));
        fprintf('δ1,poi = %.3f  δ2,poi = %.3f\n', d1_poi(k), d2_poi(k));
        fprintf('δ1,bay = %.3f  δ2,bay = %.3f\n', d1_bay(k), d2_bay(k));

        fprintf('==========================================\n\n')
        ylim([0 1.15*max(ylim)]);
    end
    % =====================================================================

    %% S-test (Here only test for M>Mc, others in vecM are not consider)
    disp(['# Calculate S-test with M >= ',num2str(vecM(k))])
    Nobsxy = zeros(length(Ggrid),1);
    index = find(time>=T_start & time<T_end & M>=Mc);
    for j=1:length(index)
        % Find smallest distance between EQS and Grids
        % Then belonging to the nearest Grids
        [~,ind] = min(rxy(index(j),:),[],2);
        Nobsxy(ind)=Nobsxy(ind)+1;
    end
    [Prob1,Prob2] = S_test(Nobsxy,Nforexy,sampleN(:,:,1));
    fprintf('Prob is :%.4f \n',Prob1);
    fprintf('==========================================\n\n')
    fprintf('Forecast Time: [%.4f  %.4f] \n',T_start,T_end);
    %% Resave the Grid data for drawing
    % Notice dimension lon and dimension lat
    Gridfore(:,:)=(reshape(Nforexy,length(Xcgrid),length(Ycgrid)))';
    Gridobs (:,:)=(reshape(Nobsxy,length(Xcgrid),length(Ycgrid)))';
end