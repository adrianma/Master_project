%% -------------------------------------------------------------------
%  function tests_main_for_Adrian_hourly
%  by Adrian Martinez Gomez
%  November 2014
%
%  Purpose:
%   * same idea as 'MainForAdrian', but this time with given population
%   prebuilt and water draw profile.
%  -------------------------------------------------------------------
function tests_main_for_Adrian_hourly(n_app,Method,tFolderMethod,N_E,...
    bDifferentFiles,Amplitudes,NAnalysis,bSame);

global eta_in runs runs_events runs_eta event_hour Events;

tFolder.load = num2str(0);
tFolder.save = num2str(NAnalysis);

% regard the case for different water draws
if(~bDifferentFiles)
    load(['VAGGELIS/saved_data/',tFolderMethod,'/',tFolder.load,...
        '/COMPARISON_sim_500EWHs.mat'],'Params','PrelimModel',...
        'WaterDrawScenarioReal','Results_comparison');
else
    % when using different files (with different profiles)
    if(bSame)
        folder_files = dir(['archived_data/baseline/',num2str(n_app),...
            'EWH/res_same/','/*.mat']);
    else
        folder_files = dir(['archived_data/baseline/',num2str(n_app),...
            'EWH/res/','/*.mat']);
    end
end

%% Simulate for different \eta values
Results = cell(1,N_E);
SS_Prec = cell(1,N_E);

for runs = 1:N_E
    
    if(bDifferentFiles)
        if(bSame)
          load(['archived_data/baseline/',num2str(n_app),...
            'EWH/res_same/',folder_files(runs).name],'Params','PrelimModel',...
            'WaterDrawScenarioReal','Results_comparison');  
        else
        load(['archived_data/baseline/',num2str(n_app),...
            'EWH/res/',folder_files(runs).name],'Params','PrelimModel',...
            'WaterDrawScenarioReal','Results_comparison');
        end
        % rename to stay by convention
        % Results_comparison = Results;
        
        % (0) store the steady-state power traj. for later analysis
        SS_Prec{runs} = Results_comparison.Prec;
    else
        SS_Prec{runs} = 0;
    end
    
    runs_eta = runs;
    runs_events = runs;
    % (1)
    % overwrite the starting and end time(s) for the simulation
    %
    % NOTE: simulate either for only 1 hour or for the remaining of the
    % day (the second one seems better(?) => to avoid off-/oversets when 
    % constructing the model)
    Params.t_init = event_hour;
    if(1)
        Params.t_sim = event_hour + 1*3600;
    else
        % Params.t_sim = 24*3600 - Params.t_init;
        Params.t_sim = 24*3600;
    end
    
    % get the state from this hour
    start_step = (Params.t_init/Params.t_sample + 1);
    Params.x_init = Results_comparison.xrec(:,start_step,:);
    Params.u_init = Results_comparison.urec(:,start_step,:);
    
    % water draws required
    Params.bWaterDraw = 1;
    
    % (3)
    % choose method (e.g. stochastic blocking; probabilistic switching)
    % for EXTERNAL CONTROLLER
    disp(['%%%%%%%%%%%% RUN NUMBER ',num2str(runs),' %%%%%%%%%%%%']);
    tic;
    Results{runs} = simulate_population(Params,PrelimModel,...
        WaterDrawScenarioReal,WaterDrawScenarioReal,Method);
    
    sim_time_CL = toc;
    disp([Method,' method runtime: ',num2str(sim_time_CL)]);
    
end

if(bSame)
    % save results in file for future analysis
    save(['archived_data/',tFolderMethod,'/',num2str(n_app),'EWH/new/',...
        tFolder.save,'/sim_',...
        num2str(n_app),'EWHs_',Method,'_',num2str(N_E),'_',...
        num2str(event_hour/3600),'hour'],'Results','N_E',...
        'eta_in','Events','event_hour','SS_Prec','Amplitudes','-v7.3');
else
    % save results in file for future analysis
    save(['archived_data/',tFolderMethod,'/',num2str(n_app),'EWH/new/',...
        tFolder.save,'/sim_',...
        num2str(n_app),'EWHs_',Method,'_',num2str(N_E),'_',...
        num2str(event_hour/3600),'hour'],'Results','N_E',...
        'eta_in','Events','event_hour','SS_Prec','Amplitudes','-v7.3');
end
end