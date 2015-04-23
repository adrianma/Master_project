%% ========================================================================
%  function main_experiments
%  by Adrian Martinez Gomez
%  November 2014
%
%  Purpose:
%   * do the experiment for multiple amplitudes as well as different files
%   (= different baseline trajectories)
%  ========================================================================
function main_SPV_Callaway(Method,hour,sigmas,Ms,NAnalysis,n_app);

global eta_in event_hour;

% % choose the folder for the method (SPV)
% if(strcmp(Method,'SetPointVariation_Callaway'))
%     tFolderMethod = 'SPV';
% end

folder_files = dir(['archived_data/baseline/',...
    num2str(n_app),'EWH/res_same/*.mat']);

% number of baseline experiments
N_E = length(folder_files);
% Folder to which the experiments are stored
% NAnalysis = 4;
count = 1;

% changes in step (how do variations go)
%
for ii = 1:length(sigmas)
    for jj = 1:length(Ms)
        
        disp(['%------------------- Sigma = ',num2str(sigmas(ii)),...
            '; M = ',num2str(Ms(jj)),' -------------------%']);
        
        % define the simulated time (=1 day)
        Params.t_sim = 24*3600;
        Params.t_sample = 10.0;
        event_hour = hour*3600;
        Params.start = event_hour;
        
        % Simulate
        tests_main_for_Adrian_hourly_Callaway(n_app,Method,N_E,...
            sigmas(ii),Ms(jj),NAnalysis);
        
        % also, remember to update the folder number and the counter
        NAnalysis = NAnalysis + 1;
        count = count + 1;
        
    end
end
end

%% -------------------------------------------------------------------
%  function tests_main_for_Adrian_hourly_Callaway
%  by Adrian Martinez Gomez
%  November 2014
%
%  Purpose:
%   * same idea as 'MainForAdrian', but this time with given population
%   prebuilt and water draw profile.
%  -------------------------------------------------------------------
function tests_main_for_Adrian_hourly_Callaway(n_app,Method,...
    N_E,Sigma,MM,NAnalysis);

global runs event_hour eta_in;

tFolder.load = num2str(0);
tFolder.save = num2str(NAnalysis);

% when using different files (with different profiles)
folder_files = dir(['archived_data/baseline/',num2str(n_app),'EWH/res_same/',...
    '/*.mat']);

%% Simulate for different \eta values
Results = cell(1,N_E);
SS_Prec = cell(1,N_E);
eta_in = cell(1,N_E);

for runs = 1:N_E
    
    load(['archived_data/baseline/',num2str(n_app),...
        'EWH/res_same/',folder_files(runs).name],'Params','PrelimModel',...
        'WaterDrawScenarioReal','Results_comparison');
    % rename to stay by convention
    % Results_comparison = Results;
    
    % (0) store the steady-state power traj. for later analysis
    SS_Prec{runs} = Results_comparison.Prec;
    
    % (1)
    % overwrite the starting and end time(s) for the simulation
    %
    % NOTE: simulate either for only 1 hour or for the remaining of the
    % day (the second one seems better(?) => to avoid off-/oversets when
    % constructing the model)
    Params.t_init = event_hour;
    Params.t_sim = event_hour + 1*3600;
    
    % get the state from this hour
    start_step = (Params.t_init/Params.t_sample + 1);
    Params.x_init = Results_comparison.xrec(:,start_step,:);
    Params.u_init = Results_comparison.urec(:,start_step,:);
    
    % water draws required
    Params.bWaterDraw = 1;
    
    % (2)
    % construct the N_E input signals (N_E: number of experiments)
    % generate the external input "ala Callaway"
    myInput = Generate_input(Sigma,MM);
    eta_in{runs} = myInput(1:361);
    close;
    
    % (3)
    % run the simulations in the population
    disp(['%%%%%%%%%%%% RUN NUMBER ',num2str(runs),' %%%%%%%%%%%%']);
    tic;
    Results{runs} = simulate_population(Params,PrelimModel,...
        WaterDrawScenarioReal,WaterDrawScenarioReal,Method);
    
    sim_time_CL = toc;
    disp([Method,' method runtime: ',num2str(sim_time_CL)]);
    
end

path_saving = ['archived_data/SPV_Callaway/new_same_2/',tFolder.save,'/sim_',...
    num2str(n_app),'EWHs_',Method,'_',num2str(N_E),'_',...
    num2str(event_hour/3600),'hour'];

% if the file is already stored, try storing in a new filename
i_addition = 1;
while(1)
    if(exist(strcat(path_saving,'.mat'),'file') == 2)
        path_saving = strcat(path_saving,num2str(i_addition));
        i_addition = i_addition + 1;
    else
        break;
    end
end

% save results in file for future analysis
save(path_saving,'Results','N_E','eta_in','event_hour','Sigma','MM',...
    'SS_Prec','-v7.3');
end