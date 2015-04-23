function [baseline_power, Results_Benchmark] = create_baseline_power(ImportResults, Params, WaterDrawScenarioReal, WaterDrawScenarioPred, PrelimModel)
%% ========================================================================
%  Creates the baseline power consumption of the EWH population without
%  external control. It is based on the daily water draw probability
%  profile and the total energy consumed by the population during the day.
%
%  Evangelos Vrettos
%  June 2012
%
%  Modified by Evangelos Vrettos
%  October 2014
%  ========================================================================

% Decide how to calculate the baseline power: (1) based on the known water draw profile;
% or (2) based on a prediction for the water draw realizations for each EWH
choice = 1;

if ImportResults == 'Y'
    load M:\Results\DataForImport\Params
    load M:\Results\DataForImport\Results_Benchmark
else
    if choice == 1
        Results_Benchmark = simulate_population(Params, PrelimModel, WaterDrawScenarioReal, WaterDrawScenarioReal, 'NoControl');
    elseif choice == 2
        Results_Benchmark = simulate_population(Params, PrelimModel, WaterDrawScenarioPred, WaterDrawScenarioPred, 'NoControl');
    end
end

% Calculate the baseline consumption using the total water draw profile
if choice == 1
    P_total = sum(Results_Benchmark.Prec(1:end-1));              % in kW
    E_total = sum(P_total * (Params.t_sample/3600));              % in kWh
    
    sumDraws = sum(WaterDrawScenarioReal.flow_rates*(Params.t_sample/3600),1);
    hourlyDrawProf = NaN(24,1);
    i=1;
    for k = 1:24
        st = 3600/Params.t_sample;
        hourlyDrawProf(k) = sum(sumDraws(i:i+st-1));
        i = i + st;
    end
    
    normHourlyDrawProf = hourlyDrawProf/sum(hourlyDrawProf);
    Hourly_power_profile = normHourlyDrawProf * E_total;               % kWh/hour
    
    NoOfIntSteps = 3600/Params.t_sample;                               % the number of simulation steps within an hour
    Power_profile_per_t_sample = NaN(1,24*NoOfIntSteps);
    
    idx = 1;
    for i = 1:24
        tmp = Hourly_power_profile(i);
        for j=0:NoOfIntSteps-1
            Power_profile_per_t_sample(idx+j) = tmp;
        end
        idx = idx+NoOfIntSteps;
    end
    
    baseline_power = [Power_profile_per_t_sample'; Power_profile_per_t_sample(end)];
    
elseif choice == 2
    Hourly_power_profile = NaN(24,1);
    i=1;
    for k = 1:24
        st = 3600/Params.t_sample;
        Hourly_power_profile(k) = mean(Results_Benchmark.Prec(i:i+st-1));
        i = i + st;
    end
    
    NoOfIntSteps = 3600/Params.t_sample;                               % the number of simulation steps within an hour
    Power_profile_per_t_sample = NaN(1,24*NoOfIntSteps);
    
    idx = 1;
    for i = 1:24
        tmp = Hourly_power_profile(i);
        for j=0:NoOfIntSteps-1
            Power_profile_per_t_sample(idx+j) = tmp;
        end
        idx = idx+NoOfIntSteps;
    end
    
    baseline_power = [Power_profile_per_t_sample'; Power_profile_per_t_sample(end)];
end

