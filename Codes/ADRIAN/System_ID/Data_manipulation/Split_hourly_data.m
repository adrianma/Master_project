%% ========================================================================
%  function Split_hourly_data
%  by Adrian Martinez Gomez
%  February 2015
%
%  Purpose:
%  ========================================================================
function [data_4,data_7,data_12,data_16,data_18,data_20,data_22] =...
    Split_hourly_data(Data);

% get the number of experiments in 'Data'
N_experiments = size(Data.ExperimentName,1);

% get a vector vHours to map indeces to hours
vHours = nan(1,N_experiments);
for ii = 1:N_experiments
    vHours(ii) = idxToHour(ii);
end

% get the experiments separated in hours
for ii = 1:N_experiments
    if(vHours(ii)==4)
        Data_4{ii} = getexp(Data,ii);
    elseif(vHours(ii)==7)
        Data_7{ii} = getexp(Data,ii);
    elseif(vHours(ii)==12)
        Data_12{ii} = getexp(Data,ii);
    elseif(vHours(ii)==16)
        Data_16{ii} = getexp(Data,ii);
    elseif(vHours(ii)==18)
        Data_18{ii} = getexp(Data,ii);
    elseif(vHours(ii)==20)
        Data_20{ii} = getexp(Data,ii);
    elseif(vHours(ii)==22)
        Data_22{ii} = getexp(Data,ii);
    end
end
% make note to remove the empty entries of the cell array
Data_4 = Data_4(~cellfun('isempty',Data_4));
Data_7 = Data_7(~cellfun('isempty',Data_7));
Data_12 = Data_12(~cellfun('isempty',Data_12));
Data_16 = Data_16(~cellfun('isempty',Data_16));
Data_18 = Data_18(~cellfun('isempty',Data_18));
Data_20 = Data_20(~cellfun('isempty',Data_20));
Data_22 = Data_22(~cellfun('isempty',Data_22));

% convert from cell-array into an iddata object 
% (not elegant, but it works)
data_4 = Data_4{1};
for ii = 2:size(Data_4,2)
    data_4 = merge(data_4,Data_4{ii});
end
data_7 = Data_7{1};
for ii = 2:size(Data_7,2)
    data_7 = merge(data_7,Data_7{ii});
end
data_12 = Data_12{1};
for ii = 2:size(Data_12,2)
    data_12 = merge(data_12,Data_12{ii});
end
data_16 = Data_16{1};
for ii = 2:size(Data_16,2)
    data_16 = merge(data_16,Data_16{ii});
end
data_18 = Data_18{1};
for ii = 2:size(Data_18,2)
    data_18 = merge(data_18,Data_18{ii});
end
data_20 = Data_20{1};
for ii = 2:size(Data_20,2)
    data_20 = merge(data_20,Data_20{ii});
end
data_22 = Data_22{1};
for ii = 2:size(Data_22,2)
    data_22 = merge(data_22,Data_22{ii});
end

end