%% ===================================================================
%  Prepares and runs EWH population simulation with Direct temperature
%  feedback control. It also produces the assessment figures.
%
%  Evangelos Vrettos
%  June 2012
%
%  Modified in October 2014, for Powertech paper simulations
%  ===================================================================
function MainForAdrian(bTests,n_app,bComp,count,Method);
% clearvars -except count bComp Method;

%% Set up simulation
Params = build_population(n_app);
% precomputes the preliminary system model
PrelimModel = precompute_system_model(Params);
% generate the stochastic water draws
WaterDrawScenarioReal = build_normalized_draw_scenario(Params);

% choose method (e.g. stochastic blocking; probabilistic switching, etc.)
% for EXTERNAL CONTROLLER
tic;
Results = simulate_population(Params,PrelimModel,WaterDrawScenarioReal,...
    WaterDrawScenarioReal,Method);
sim_time_CL = toc;

% for comparison with the no control (only INTERNAL CONTROLLERS) method
if(bComp)
    tic;
    Results_comparison = simulate_population(Params,PrelimModel,....
        WaterDrawScenarioReal,WaterDrawScenarioReal,'NoControl');
    sim_time_NC = toc;
    disp(['NoControl runtime: ',num2str(sim_time_NC)]);
end

disp([Method,' method runtime: ',num2str(sim_time_CL)]);

%% Save workspace
% if outside loop involved
if(~exist('count','var') || count==0)
    count = 1;
end

% if tests are bein performed
if(bTests)
    global eta_in runs_eta runs_events;
    if(strcmp(Method,'ProbSwitching'))
        save(['VAGGELIS/saved_data/PS/sim_',num2str(Params.n_app),...
            'EWHs_',Method,'_',num2str(count),'_event',...
            num2str(runs_events),'_eta_',...
            strrep(num2str(eta_in(runs_eta)), '.', '')]);
    elseif(strcmp(Method,'SetPointVariation'))
        save(['VAGGELIS/saved_data/SPV/refined/sim_',num2str(Params.n_app),...
            'EWHs_',Method,'_',num2str(count),'_event',...
            num2str(runs_events),'_eta_',...
            strrep(num2str(eta_in(runs_eta)), '.', '')]);
    end
else
    save(['VAGGELIS/saved_data/open_loop_sim_',num2str(Params.n_app),...
        'EWHs_',Method,'_',num2str(count)]);
end

end