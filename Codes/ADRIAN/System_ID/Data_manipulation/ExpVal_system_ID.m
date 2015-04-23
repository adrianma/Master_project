%% ========================================================================
%  function ExpVal_system_ID
%  by Adrian Martinez Gomez
%  February 2015
%
%  Purpose:
%  Split the Data into two subsets: Experimental and Validation.
%  These subsests are needed to later on create the model(s).
%  ========================================================================
function [Experimental,Validation] = ExpVal_system_ID(Data,bNoSign);
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
vRand_pos = randperm(length(vIdx_pos));

vIdx_neg = find(vMask_sign == 0);
vRand_neg = randperm(length(vIdx_neg));

vIdx_zero = find(vMask_sign == 2);
vRand_zero = randperm(length(vIdx_zero));

% experimental data
exp_pos = getexp(Data,vIdx_pos(vRand_pos(1:(length(vIdx_pos)/2))));
exp_neg = getexp(Data,vIdx_neg(vRand_neg(1:(length(vIdx_neg)/2))));

if(~isempty(vRand_zero))
    exp_zero = getexp(Data,vIdx_zero(vRand_zero(1:(length(vIdx_zero)/2))));
end

% validation data
val_pos = getexp(Data,vIdx_pos(vRand_pos((length(vIdx_pos)/2 + 1):end)));
val_neg = getexp(Data,vIdx_neg(vRand_neg((length(vIdx_neg)/2 + 1):end)));
if(~isempty(vRand_zero))
    val_zero = getexp(Data,vIdx_zero(vRand_zero((length(vIdx_zero)/2 + 1):end)));
end

% merge the positive and negative ones
if(bNoSign)
    Experimental = merge(exp_pos,exp_neg);
    Validation = merge(val_pos,val_neg);
% or don't
else
    Experimental.pos = exp_pos;
    Experimental.neg = exp_neg;
    
    Validation.pos = val_pos;
    Validation.neg = val_neg;
end

if(~isempty(vRand_zero))
    Experimental = merge(Experimental,exp_zero);
    Validation = merge(Validation,val_zero);
end

end