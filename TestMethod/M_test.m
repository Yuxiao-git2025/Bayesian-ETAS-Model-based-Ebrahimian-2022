% =========================================================================
%                            Magnitude-Test
%                             Yu.2025/10/9
%                          Refer to Savran.2020
% =========================================================================
% Based on the statistic on a square metric computed from the difference in
% logarithms between the incremental MFDs of the so-called union catalog(1),
% individual/sample catalogs(2), and the observed catalog(3)
% =========================================================================
% Using the logarithm of bin-wise magnitude counts places greater weight on
% relatively less frequent earthquakes, particularly those of larger magnitudes.
% || Reason ||: The magnitude-frequency distribution (MFD) of earthquakes
% follows a power-law relationship, known as the Gutenberg-Richter law,
% wherein the occurrence frequency of small-magnitude events significantly
% exceeds that of large-magnitude events. When absolute error metrics—such
% as differences in event counts are used directly, minor discrepancies
% in the low-magnitude bins tend to dominate the overall statistical
% evaluation due to their high event counts, thereby obscuring substantial
% deviations in the high-magnitude range
% And add unity to each bin to prevent the singularity associated with log(0)
% =========================================================================
function [kapa, dobs, D]=M_test(Obs ,Forecast,mbin)
% ===== 1. Initialize =====
Ns=numel(Forecast);
D=zeros(Ns, 1);
Valid=true(Ns, 1);
Method=2;
switch Method
    case 1
        % ===== 2. Forecasting Hist Expectation =====
        CombMag = cell2mat(Forecast);
        Nexp=histcounts(CombMag, mbin, 'Normalization', 'pdf');

        % ===== 3. Observation Hist Expectation =====
        Nobs=histcounts(Obs, mbin, 'Normalization', 'pdf');


        % ===== 4. Statistics =====
        % CDF
        ObsCDF = cumsum(Nobs) / sum(Nobs);
        ForeCDF= cumsum(Nexp) / sum(Nexp);

        % Statistics (overall or expectation)
        dobs = sum((ObsCDF - ForeCDF).^2);

        % ===== 5. Simulation cell =====
        for i = 1:Ns
            Magcell = Forecast{i};
            if isempty(Magcell)
                Valid(i) = false;
                continue;
            end
            Ncell=histcounts(Magcell, mbin, 'Normalization', 'pdf');
            CellCDF=cumsum(Ncell) / sum(Ncell);
            % Statistics (each cell)
            D(i) = sum((CellCDF - ForeCDF).^2);
        end
        D=D(Valid);
        %%


    case 2
        % <1>.Union Statistics
        % Sum of squared logarithmic residuals between the normalized
        % observed magnitude and union histograms
        TotalForecast=cell2mat(Forecast);
        Nexp=histcounts(TotalForecast, mbin, 'Normalization', 'count');
        Nobs=histcounts(Obs, mbin, 'Normalization', 'count');
        dobs=sum( (log( length(Obs)/length(TotalForecast)*Nexp+1 )-log(Nobs+1) ).^2 );
        % <2>.Sample Statistics
        for i = 1:Ns
            Magcell = Forecast{i};
            if isempty(Magcell)
                Valid(i) = false;
                continue;
            end
            Ncell=histcounts(Magcell, mbin, 'Normalization', 'count');
            D(i)=sum( ( log( length(Obs)/length(TotalForecast)*Nexp+1 )- ...
                log( length(Obs)/length(Magcell)*Ncell+1 ) ).^2 );
        end
        D=D(Valid);
end
% ===== end. κ =====
kapa=sum(D <= dobs) / numel(D);

end