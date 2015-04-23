function Results = simulate_population(Params, PrelimModel, WaterDrawScenarioReal, WaterDrawScenarioAgg, Method)
%% =====================================================================
%  Simulate the water heater population and save results
%  Stephan Koch
%  Oct 17, 2010
%  Modified by Martin Pfeiffer, June 2011
%
%  Modified by Evangelos Vrettos
%  April 2012
%
%  Modified by Evangelos Vrettos
%  October 2014
% 
%  Modified by Adrian Martinez Gomez
%  October 2014
%
%  Modified by Adrian Martinez Gomez
%  January 2015
%  =====================================================================

%% Initialize the records
timeSteps   = length(Params.t_init:Params.t_sample:Params.t_sim);
numStates   = size(Params.x_init,1);
numInputs   = size(Params.u_init,1);
xrec        = nan(numStates,timeSteps,Params.n_app);                        % states (temperatures)
SOC         = nan(Params.n_app,timeSteps);
Umirec      = zeros(numStates-2,timeSteps,Params.n_app);
urec        = nan(numInputs,timeSteps,Params.n_app);                        % external control signals
Psetrec     = nan(1,timeSteps);
Prec        = nan(1,timeSteps);
mdotrec     = nan(1,timeSteps,Params.n_app);
step        = 1;

PrelimModelAgg = PrelimModel;

% load the optimal control (if provided) from the SQP
if(isstruct(Method))
    temp_name = Method.name;
    temp_u_optimal = Method.u_optimal;
    clear Method;
    
    Method = temp_name;
    u_external_optimal = temp_u_optimal;
    clear temp_name temp_u_optimal;
end

%%
if(0)
if(strcmp(Method,'StochBlocking') || strcmp(Method,'ProbSwitching') || ...
        strcmp(Method,'SetPointVariation'))
    blocking_old = nan(1,Params.n_app);
    n_rand = rand(1,Params.n_app);
   
    % needed for script for baseline power(?)
    % WaterDrawScenarioPred = build_normalized_draw_scenario(Params);
    %
    % get the desired setpoint trajectory
    % [baseline_power,~] = create_baseline_power('N',Params,...
    %    WaterDrawScenarioReal,WaterDrawScenarioPred,PrelimModel);
    
    vTime = Params.t_init:Params.t_sample:Params.t_sim;
    [avP_set] = runs_P_plot(0,0);
    
    [fitresult,gof] = createFit_P_base(vTime,avP_set);
    baseline_power = fitresult(vTime);
    
    % create the LFC signal, which is base_part + variable part
    %
    % TODO: variable part changes (suggested by Maryam)!
    LFC_signal = create_LFC_signal(baseline_power,Params);
    
    P_set = LFC_signal;
    
    fprintf('hello\n');
else
    blocking_old = nan(1,Params.n_app);
    n_rand = rand(1,Params.n_app);
    
    [P_set,~] = construct_P_set(Params.n_app,0);
    
    % random number for 'ProbSwitching' or 'SetPointVariation'
    etas_rnd = nan(1,(Params.t_sim-Params.t_init)*Params.t_sample + 1);
end
else
    etas_rnd = nan(1,(Params.t_sim-Params.t_init)*Params.t_sample + 1);
end

%% Running the simulation over simulation period
for t = Params.t_init:Params.t_sample:Params.t_sim
    
    % Output sim time
    if mod(t/3600, 0.1) < 0.0001,
        disp(['Time: ', num2str(t/3600), ' h'])
    end
    
    % Load the water draw (WaterDrawScenario.flow_rates is in kg/h, so here it is transformed to kg/sec)
    mdot = WaterDrawScenarioReal.flow_rates(:,(t/Params.t_sample + 1))/3600;
    % with this bool variable one can select the case with no water draws
    if(~Params.bWaterDraw)
        mdot = zeros(size(mdot));
    end
    
    % Build the final model now that the water draw is known
    % ...(change)
    if(0)
        rho = Params.rho;
        n = Params.n;
        k = Params.k;
        c = Params.c;
        d = Params.d;
        x_current = Params.xcurrent;
        diam = Params.diam;
        turbMixing = Params.turbMixing;
        h = Params.h;
        g = Params.g;
    end
    Model = build_final_system_model(Params, PrelimModel, mdot);
    
    % Calculate the on/off mode for each EWH (differentiate depending on control strategy)
    switch Method
        case 'NoControl'
            unew = NoControl(Params, step, t, xrec, urec);
        case 'StochBlocking'
            [unew,blocking_new] = StochBlocking(Params,step,t,xrec,...
                urec,P_set,blocking_old,n_rand);
            blocking_old = blocking_new;
        case 'ProbSwitching'
            [unew,eta] = ProbSwitching(Params,step,t,xrec,urec);
            etas_rnd(step) = eta;
        case 'SetPointVariation'
            [unew,eta] = SetPointVariation(Params,step,t,...
                xrec,urec);
            etas_rnd(step) = eta;
        case 'SetPointVariationCL'
            [unew,eta] = SetPointVariationCL(Params,step,t,...
                xrec,urec,u_external_optimal);
            etas_rnd(step) = eta;
        case 'SetPointVariation_Callaway'
            [unew,eta] = SetPointVariation_Callaway(Params,step,t,...
                xrec,urec);
            etas_rnd(step) = eta;
    end
    
    % Update the states of the population
    [xnew, Uminew] = updateStates(step, Params, Model, xrec, unew, Umirec);
    
    % Update records
    xrec(:,step,:)      = xnew(:,1,:);
    Umirec(:,step,:)    = Uminew(:,1,:);
    urec(:,step,:)      = unew(:,1,:);
    mdotrec(:,step,:)   = mdot(:,:);
    tmp1                = sum(Params.P1_el(urec(1,step,:) == 1));
    tmp2                = sum(Params.P1_el(urec(2,step,:) == 1));
    Prec(step)          = tmp1+tmp2;
    
    % for now, we don't need the state of charge (SOC)
    if(0)
        for i = 1:Params.n_app
            SOC(i,step) = SOC_calc(Params, i, xrec(:,step,:));
        end
    end
    
    % Update time index
    step = step + 1;
end

%% Assemble result structure
Results.urec = urec;
Results.xrec = xrec;
Results.SOC = SOC;
Results.mdotrec = mdotrec;
Results.Umirec = Umirec;
Results.Prec = Prec;
Results.Psetrec = Psetrec;
% for stochastic blocking
if(strcmp(Method,'StochBlocking'))
    Results.n_rand = n_rand;
end

global runs_events Events;

% for probabilistic switching
if(strcmp(Method,'ProbSwitching'))
    temp = etas_rnd(~isnan(etas_rnd));
    Results.etas = temp;
    %Results.P_set = P_set;
    Results.Event = Events{runs_events};
end
% for set-point variation
if(strcmp(Method,'SetPointVariation'))
    temp = etas_rnd(~isnan(etas_rnd));
    Results.etas = temp;
    %Results.P_set = P_set;
    Results.Event = Events{runs_events};
end
% for closed-loop set-point variation
if(strcmp(Method,'SetPointVariationCL'))
    temp = etas_rnd(~isnan(etas_rnd));
    Results.etas = temp;
end
% for set-point variation "ala Callaway"
if(strcmp(Method,'SetPointVariation_Callaway'))
    temp = etas_rnd(~isnan(etas_rnd));
    Results.etas = temp;
end
end