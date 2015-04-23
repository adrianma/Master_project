%% ========================================================================
%  function PosNeg_system_ID
%  by Adrian Martinez Gomez
%  March 2015
%
%  Purpose:
%  Split the Data into two subsets: positive and negative external input
%  These subsests are needed to later on create the model(s).
%  ========================================================================
function [DataOut_pos,DataOut_neg] = PosNeg_system_ID(Data);

% number of data input-output pairs
N_Data = size(Data.ExperimentName,1);

% this is a mask to get the experiments with positive external inputs.
% with this we not only choose randomly the experiments, but we guarantee
% that there will be the same amount of positive and negative ones.
vMask_sign = nan(1,N_Data);
for ii = 1:N_Data
    val = sum(Data.u{ii}(:,1));
    
    if(val > 0)
        vMask_sign(ii) = 1;
    elseif(val < 0)
        vMask_sign(ii) = 0;
    else 
        vMask_sign(ii) = 2;
    end
end

vIdx_pos = find(vMask_sign == 1);
vIdx_neg = find(vMask_sign == 0);

DataOut_pos = getexp(Data,vIdx_pos);
DataOut_neg = getexp(Data,vIdx_neg);

end
