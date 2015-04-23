%% ========================================================================
%  function Merge_hourly_datasets
%  by Adrian Martinez Gomez
%  February 2015
%
%  Inputs: datasets Data_1 and Data_2, which contain the full experiments
%  for all the different hours.
%
%  Output: OutData (merged Data_1, Data_2, sorted hourly).
%
%  Purpose:
%  merge 2 datasets into a new one (with double as many experiments)
%
%  ========================================================================
function OutData = Merge_hourly_datasets(Data_1,Data_2);

Data_2 = Adapt_experiment_number(Data_2);

% 1) split the data into hours
[data_4_1,data_7_1,data_12_1,data_16_1,data_18_1,data_20_1,data_22_1] = ...
    Split_hourly_data(Data_1);

% 2) split the data into hours
[data_4_2,data_7_2,data_12_2,data_16_2,data_18_2,data_20_2,data_22_2] = ...
    Split_hourly_data(Data_2);

% merge the experiments
temp_12 = merge(data_12_1,data_12_2);
temp_16 = merge(data_16_1,data_16_2);
temp_18 = merge(data_18_1,data_18_2);
temp_20 = merge(data_20_1,data_20_2);
temp_22 = merge(data_22_1,data_22_2);
temp_4 = merge(data_4_1,data_4_2);
temp_7 = merge(data_7_1,data_7_2);

OutData = merge(temp_12,temp_16);
OutData = merge(OutData,temp_18);
OutData = merge(OutData,temp_20);
OutData = merge(OutData,temp_22);
OutData = merge(OutData,temp_4);
OutData = merge(OutData,temp_7);

end