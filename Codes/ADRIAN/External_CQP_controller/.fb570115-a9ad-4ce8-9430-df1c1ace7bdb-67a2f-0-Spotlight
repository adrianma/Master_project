%% ========================================================================
%  function simulate_on_real_system
%  by Adrian Martinez Gomez
%  January 2015
%
%
%  Purpose:
%   * apply an input trajectory computed by the (C)QP optimization problem
%   to the real system
%  ========================================================================
function [ResultsOptimal,Params,SS_Prec] = ...
    simulate_on_real_system(N_zeros,ref_exo_idx,u_optimal);

global event_hour;

% convert the idx into an hour, and then the hour into second of the day
event_hour = idxToHour(ref_exo_idx);
event_hour = event_hour*3600;

experiment_number = 1 + mod(ref_exo_idx - 1,20);

% (0) 
% load the scenario needed
load(['archived_data/baseline/1000EWH/res_same/',...
    'sim_1000EWHs_runNumber',num2str(experiment_number),'.mat'],...
    'Params','PrelimModel','WaterDrawScenarioReal','Results_comparison');

% (1) store the steady-state power traj. for later analysis
SS_Prec = Results_comparison.Prec;

% (2)
% overwrite the starting and end time(s) for the simulation
%
% NOTE: take the extra zeros at the beginning/end
Params.t_sample = 10;
Params.t_init = event_hour - N_zeros*Params.t_sample;
Params.t_sim = event_hour + 1*3600 + (3*N_zeros-1)*Params.t_sample;

% get the state for this hour
start_step = (Params.t_init/Params.t_sample + 1);
Params.x_init = Results_comparison.xrec(:,start_step,:);
Params.u_init = Results_comparison.urec(:,start_step,:);

% (3)
% choose method (e.g. stochastic blocking; probabilistic switching)
% for EXTERNAL CONTROLLER
Method.name = 'SetPointVariationCL';
Method.u_optimal = u_optimal;

disp(['%%%%%%%%%%%% SIMULATE OPTIMAL %%%%%%%%%%%%']);
tic;
ResultsOptimal = simulate_population(Params,PrelimModel,...
    WaterDrawScenarioReal,WaterDrawScenarioReal,Method);
sim_time_CL = toc;
disp([Method.name,' method runtime: ',num2str(sim_time_CL)]);

end