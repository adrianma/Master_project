%% ===================================================================
%  function main_No_Ext_Control
%  Adrian Martinez Gomez
%  April 2015
%
%  Inputs:
%      n_app : number of EWHs
%      N_experiments : number of simulation experiments to perform
%
%  Purpose:
%  This function generates the baseline experiments (this is no external
%  control; autonomous case). It is done for n_app EWHs and for
%  N_experiments different simulations
%  ===================================================================
function [] = main_No_Ext_Control(n_app,N_experiments)

% Set up the population
Params = build_population(n_app);
Params.bWaterDraw = 1;
% precomputes the preliminary system model
PrelimModel = precompute_system_model(Params);
% select method -> autonomous
Method = 'NoControl';

% define the path to check if folder already exists; if not, create it!
my_path = strcat(strcat('archived_data/baseline/',num2str(n_app)),'EWH');
if(~(exist(my_path,'dir')==7))
    mkdir(my_path);
    addpath(my_path);
    mkdir(strcat(my_path,'/res_same'));
    addpath(strcat(my_path,'/res_same'));
end

%% Run the N_experiments simulations for the baseline
for ii = 1:N_experiments
    
    % generate the stochastic water draws
    WaterDrawScenarioReal = build_normalized_draw_scenario(Params);
    
    % simulate the system
    Results = simulate_population(Params,PrelimModel,...
        WaterDrawScenarioReal,WaterDrawScenarioReal,Method);
    
    % save results in file for future analysis
    save(['archived_data/baseline/',num2str(n_app),'EWH/res_same/sim_',...
        num2str(n_app),'EWHs_runNumber',num2str(ii)],'Params','Results',...
        'PrelimModel','WaterDrawScenarioReal','-v7.3');
    
end

end