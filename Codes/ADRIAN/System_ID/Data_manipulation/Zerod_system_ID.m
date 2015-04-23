%% ========================================================================
%  function Zerod_system_ID
%  by Adrian Martinez Gomez
%  March 2015
%
%  Inputs: 
%
%  Purpose:
%  Add "N_zeros" steps BEFORE the external input change happens.
%  (Of course, at the beginning of the data)
%  ========================================================================
function zData = Zerod_system_ID(Data,N_zeros);
%% 0)
% boolean which is equal to 1 if the model has an exogenous input part
bWexo = (size(Data.u{1},2) == 2);

%% 1) load the signals from the 20 experiments + if needed exo. input

% load the power signal from the 20 (baseline) experiments
load('archived_data/baseline/1000EWH/same_baseline_data/power_baseline_same.mat',...
    'y_diffs_norm');

% load the exogenous input (averaged water draw)
if(bWexo)
    load(['archived_data/baseline/1000EWH/same_baseline_data/mean_water_1000EWH_same.mat'...
        ],'mean_water');
    exogenous_input = mean_water;
    clear mean_water;
end

%% 2)
% number of data input-output pairs
N_Data = size(Data.ExperimentName,1);
% initializations
yy = cell(1,N_Data);
uu = cell(1,N_Data);
TTss = cell(1,N_Data);

for data_idx = 1:N_Data
    
    event_hour = 3600.*idxToHour(data_idx);
    low_idx = event_hour/10 + 1;
    
    yFill = y_diffs_norm(1+mod(data_idx-1,20),((low_idx-N_zeros):(low_idx-1))).';
    
    % add zeros to the output
    yy{data_idx} = [yFill;Data.y{data_idx}];
    
    % add the zeros to the inputs. 
    % (also take into account the exogenous if present)
    if(bWexo)
        uFill = exogenous_input((low_idx - N_zeros):(low_idx - 1)).';
        
        uu{data_idx} = [[zeros(N_zeros,1);Data.u{data_idx}(:,1)],...
            [uFill;Data.u{data_idx}(:,2)]];
    else
        uu{data_idx} = [zeros(N_zeros,1);Data.u{data_idx}];
    end
    
    % copy the period from the datasets
    TTss{data_idx} = Data.Ts{data_idx};
end

% get the iddata object
zData = iddata(yy,uu,TTss);

end