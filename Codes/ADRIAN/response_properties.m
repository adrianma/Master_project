%% ===================================================================
%  function response_properties
%  by Adrian Martinez Gomez
%  October 2014
%
%  Purpose:
%   * get some properties for the responses for a given event:
%       maximum/minimum value
%       rise-time
%       (do settling time and undershoot by inspection of the plots)
%   * save/display the results in 'Properties_Results'
%  ===================================================================

function Properties_Results = response_properties(Method,NFolder);
%% Read the files names for later use

% choose the folder for the method (PS:prob. switching,
% SPV:set-point variation)
if(strcmp(Method,'ProbSwitching'))
    tFolderMethod = 'PS';
elseif(strcmp(Method,'SetPointVariation'))
    tFolderMethod = 'SPV';
end
% get the subfolder for the simulation case
if(NFolder~=0)
    tFolder = strcat('/',num2str(NFolder));
else
    tFolder = '/';
end
% get the file names
folder_files = dir(['VAGGELIS/saved_data/',tFolderMethod,tFolder,...
    '/*.mat']);
load(folder_files(1).name);

% define time vector and confidence interval
vTime = Params.t_init:Params.t_sample:Params.t_sim;
confidence = 0.05;

%% Iterate through all the files and plot the response
for ii = 2:length(folder_files)
    
    % load data-set
    load(folder_files(ii).name);
    disp(folder_files(ii).name);
    
    % store the eta (for later on use)
    Properties_Results.eta(ii) = Results.etas(1);
    
    % store the time of the event for the lower subplot
    step_event_st = Results.Event(1)/Params.t_sample;
    step_event_end = Results.Event(end)/Params.t_sample;
    
    %% get the maximum/minumum value achieved
    %
    % differentiate between positive (max) and negative (min) etas
    % add 50 steps in case the value is outside the event itself (lag)
    if(Results.etas(1) > 0)
        [val,pos] = max(Results.Prec(step_event_st:(step_event_end+50)));           
    else
        [val,pos] = max(Results.Prec(step_event_st:(step_event_end+50)));
    end
    Properties_Results.extrema(ii-1) = val;
    Properties_Results.rise_time(ii-1) = pos + (step_event_st - 1);
    
    %% matlab's built in step response tools
    Properties_Results.S(ii-1) = stepinfo(Results.Prec,vTime);
end