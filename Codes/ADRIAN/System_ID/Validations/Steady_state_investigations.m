%% ========================================================================
%  function Steady_state_investigations.m
%  by Adrian Martinez Gomez
%  February 2015
%
%  Output:
%  Input:
%
%  Purpose:
%  ========================================================================
function Steady_state_investigations(ss,Nmodel);

%% 1)
% ...(change)
if(0)
    % load the state-space linear models
    load('archived_data/SPV/1000EWH/Models/MODELS_SUM_AVERAGED_1000_02022015.mat',...
        'ss4_n4sid_sim');
    load('archived_data/SPV/1000EWH/Models/MODELS_SUM_AVERAGED_1000_04022015.mat',...
        'ss4_n4sid_sim_ext');
end

% boolean which is equal to 1 if the model has an exogenous input part
bWExo = (size(ss.b,2) == 2);

% load the exogenous input signal
load('archived_data/baseline/1000EWH/mean_water_P_agg_1000EWH.mat',...
    'mean_water');
exogenous_input = mean_water.signal.';
clear mean_water;

% Construct input (with/without exogenous) u
if(bWExo)
    u_inputs = [zeros(8641,2)];
    u_inputs_w = [zeros(8641,1),exogenous_input];
else
    u_inputs = [zeros(8641,1)];
end

% Simulate the linear system for the contidions and see output
t = 0:10:3600*24;
x0 = zeros(size(ss.a,1),1);
[y_,~,~] = lsim(ss,u_inputs,t,x0);
if(bWExo)
    [y_w,~,~] = lsim(ss,u_inputs_w,t,x0);
end

% plot the responses
figure;
subplot(2,1,1);
plot(y_);
grid on;
xlabel('Step index k');
ylabel(['Model ',num2str(Nmodel)]);
legend('y_{ss}','Location','Best');
title('External input u_{k}=0  and exogenous input w_{k}=0');

if(bWExo)
    subplot(2,1,2);
    plot(y_w);
    grid on;
    xlabel('Step index k');
    ylabel(['Model ',num2str(Nmodel)]);
    legend('y_{b}','Location','Best');
    title('External input u_{k}=0  and exogenous input w_{k}\neq0');
    
    % plot the baseline of the model + the averaged baseline consumption 
    % from the real system
    load('archived_data/baseline/1000EWH/power_baseline.mat', 'y_baseline_norm');
    figure;
    plot(y_w.' + y_baseline_norm);
    grid on;
    xlabel('Step index k');
    ylabel(['Model ',num2str(Nmodel)]);
    title(['Baseline of Model ',num2str(Nmodel),...
        ' (u_{k}=0 w_{k}\neq 0) + averaged baseline consumption']);
    ylim([0,0.3]);
end


disp('DONE part 1)!');

%% 2)
% simulate a new scenarion (takes time to simulate real system ~20 min
if(0)
    if(0)
        if(0)
            n_app = 1000;
            Params = build_population(n_app);
            PrelimModel = precompute_system_model(Params);
            WaterDrawScenarioReal = build_normalized_draw_scenario(Params);
        else
            idx = randi(20);
            load(['archived_data/baseline/1000EWH/res/sim_1000EWHs_runNumber',num2str(idx),'.mat'],...
                'Params','PrelimModel','WaterDrawScenarioReal','Results_comparison');
        end
        
        Method = 'NoControl';
        Results = simulate_population(Params,PrelimModel,...
            WaterDrawScenarioReal,WaterDrawScenarioReal,Method);
    end
    
    % load the scenario already calculated
    if(0)
        load('archived_data/SPV/matlab.mat','Results','Params','idx');
        load(['archived_data/baseline/1000EWH/res/sim_1000EWHs_runNumber',num2str(idx),'.mat'],...
            'Params','PrelimModel','WaterDrawScenarioReal','Results_comparison');
    end
    
    figure;
    subplot(2,1,1);
    plot(Results.Prec./sum(Params.P1_el));
    grid on;
    xlabel('Step index k');
    ylabel('Real system');
    title('External input u_{k}=0  and exogenous input w_{k}=0');
    
    subplot(2,1,2)
    plot(Results_comparison.Prec./sum(Params.P1_el));
    grid on;
    xlabel('Step index k');
    ylabel('Real system');
    title('External input u_{k}=0  and exogenous input w_{k}\neq0');
    
    
    disp('DONE part 2)!');
end

end