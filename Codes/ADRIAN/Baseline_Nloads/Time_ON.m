%% ========================================================================
%  function Loads_N_ON
%  by Adrian Martinez Gomez
%  March 2015
%
%  Purpose:
%  * Generate content for Table 3.2 of the written Master Thesis
%  * The number displayed is for ONE EWH
%
%  ========================================================================
function Time_ON();
for ii = 1:20
    idx = ii;
    load(['archived_data/baseline/1000EWH/res_same/sim_1000EWHs_runNumber',num2str(idx),'.mat']);
    
    Switches = squeeze(Results_comparison.urec(1,:,:)).';
    
    TimeON = nan(1,n_app);
    for ii = 1:n_app
        TimeON(ii) = sum(Switches(ii,:))/length(Switches(ii,:));
    end
    % Convert to [%]
    TimeON = TimeON.*100;
    
    meanTimeON = mean(TimeON);
    stdTimeON = std(TimeON);
    
    disp([' => Mean = ',num2str(meanTimeON),'[%] and Std = ',...
        num2str(stdTimeON),'[%]']);
end
end