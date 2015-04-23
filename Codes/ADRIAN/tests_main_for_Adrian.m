%% -------------------------------------------------------------------
%  function myMainForAdrian
%  by Adrian Martinez Gomez
%  October 2014
%
%  Purpose:
%   * same idea as 'MainForAdrian', but this time with given population
%   prebuilt and water draw profile.
%  -------------------------------------------------------------------
function tests_main_for_Adrian(n_app,Method,tFolderMethod,NAnalysis);

global eta_in runs_eta runs_events Events;

tFolder = strcat('/',num2str(NAnalysis));

load(['/VAGGELIS/saved_data/',tFolderMethod,tFolder,'/COMPARISON_sim_500EWHs.mat'],...
    'Params','PrelimModel','WaterDrawScenarioReal');

%% Simulate for different \eta values
for runs_events = 1:length(Events)
    for runs_eta = 1:length(eta_in)

        % choose method (e.g. stochastic blocking; probabilistic switching)
        % for EXTERNAL CONTROLLER
        tic;
        Results = simulate_population(Params,PrelimModel,...
            WaterDrawScenarioReal,WaterDrawScenarioReal,Method);
        
        sim_time_CL = toc;
        disp([Method,' method runtime: ',num2str(sim_time_CL)]);
        
        % save workspace in file for future analysis
        save(['VAGGELIS/saved_data/',tFolderMethod,tFolder,'/sim_',...
            num2str(n_app),'EWHs_',Method,'_event',...
            num2str(runs_events),'_eta_',...
            strrep(num2str(eta_in(runs_eta)), '.', '')]);
    end
end
end