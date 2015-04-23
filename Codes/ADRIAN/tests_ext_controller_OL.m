%% ===================================================================
%  function tests_set_point_variation
%  by Adrian Martinez Gomez
%  October 2014
%   
%  Inputs:
%   *Method: string containing the method to be tested (ProbSwitching,
%   SetPointVariation, etc)
%   *bSimulate: boolean, set to 1 if want to simulate
%   *bAnalyze: boolean, set to 1 if want to analyze the simulated data
%   *sPlot: giving (if plots desired) what plot script to be runed
%   *NAnalysis: giving the new folder name for the simulations to be
%   stored. Set to "0" if no new simulation set needed.
%
%  Purpose:
%   *
%   *
%  ===================================================================
function tests_ext_controller_OL(Method,bSimulate,NAnalysis);

% define the simulated time (=1 day)
Params.t_sim = 86400;
% load the etas for this experiment
load('ADRIAN/etas.mat');
global eta_in;

% choose the folder for the method (PS:prob. switching,
% SPV:set-point variation)
if(strcmp(Method,'ProbSwitching'))
    tFolderMethod = 'PS';
    %eta_in = etas_prob_switching;
    eta_in = 0.1:0.2:0.9;
elseif(strcmp(Method,'SetPointVariation'))
    tFolderMethod = 'SPV';
    eta_in = etas_set_point_var;    
end

% events: where the changes take place (hours 4, 12 , 18)
%
% define events withing all hours
global runs_events Events;
Events = cell(1,1);
event_hour = 18*3600;
Params.hours = event_hour;
Params.durations = 0.2*3600;

for ii = 1:length(Events)
    Events{ii} = Params.hours(ii):1:(Params.hours(ii) + Params.durations(ii));
end
% step for 12 minutes time (starting at hour 18)
%Events{1} = 3*Params.t_sim/4 : 18.2*3600;

% define number of loads to simulate the method with
n_apps = [500];

%% Simulate
if(bSimulate)
  
    global runs_eta;
    runs_eta = 1;
    runs_events = 1;
    
    for ii = 1:length(n_apps)
        tests_main_for_Adrian(n_apps(ii),Method,tFolderMethod,NAnalysis);
    end
end
    
end