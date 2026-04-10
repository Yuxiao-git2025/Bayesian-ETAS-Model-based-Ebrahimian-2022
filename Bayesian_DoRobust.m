% =========================================================================
% The loop takes very long time and you can choose some samples to try
% When generate forecasting EQS, do not use vecM, and only calculate the
% Expsamples in each grids(N-test), it should be considered
function [Expsamples]=Bayesian_DoRobust(DoRobust,samples,indexSeq,OutDir)
global time_MS  M  Mc tstart tend rxy dA  ...
      Muxy2col  Ngrid vecM

if DoRobust == 1
fprintf('# Process Robust estimation for the Expected Number Events \n') ;

Ms=cell(1,size(samples,2));
ts=cell(1,size(samples,2));
lons=cell(1,size(samples,2));
lats=cell(1,size(samples,2));
rs=cell(1,size(samples,2));
Intlam_t=cell(1,size(samples,2));
% =========================================================================
%                                  Main Loop
% Nsample=size(samples,2);
% DIY here
Nsample0=1;
Nsample=20;   

waithandel=waitbar(Nsample0-1,sprintf(['Process samples ' num2str(Nsample-Nsample0+1)]));
Expsamples=zeros(Ngrid,(Nsample-Nsample0+1),length(vecM)); % like: 6e4*371*5
for j = Nsample0: Nsample
    waitbar((j-Nsample0+1)/(Nsample-Nsample0+1),waithandel);
    fprintf('>> Generating Seq in %d sample \n ',(j)); 
    % Generate seqg
    [Ms{1,j},ts{1,j},lons{1,j},...
        lats{1,j},rs{1,j},Intlam_t{1,j}] ...
    =generateSEQ(M(indexSeq),time_MS(indexSeq),rxy(indexSeq,:),samples(:,j));
    
    % Calculate Robust N
    fprintf('# Begin to Calculate N \n');
    for k=1:length(vecM) % mag sequence like [3 4 5 6...]
        for ngrid=1:Ngrid            
            if ~isempty(Ms{1,j})
                % If an aftershock sequence is generated, then calculate 
                % the expected number of events including the background 
                % and the aftershocks
                Expsamples(ngrid,j,k) = calculate_N ...
                    ([M(indexSeq);Ms{1,j}],[time_MS(indexSeq);ts{1,j}], ...
                    [rxy(indexSeq,ngrid);rs{1,j}(:,ngrid)],vecM(k),tstart,tend,Mc, ...
                    samples(:,j),Intlam_t{1,j}(:,ngrid),dA,Muxy2col(ngrid));
            else
                % No aftershock sequence
                Expsamples(ngrid,j,k) = calculate_N...
                    (M(indexSeq),time_MS(indexSeq),rxy(indexSeq,ngrid),vecM(k), ...
                    tstart,tend,Mc,samples(:,j),[],dA,Muxy2col(ngrid));
            end
        end % Grid end here
    end % mag-seq end here
    
end % samples end here
delete(waithandel);
filename = fullfile(OutDir,'\robust estimate N.mat');
save(filename,'Expsamples','-v7.3')
        
end