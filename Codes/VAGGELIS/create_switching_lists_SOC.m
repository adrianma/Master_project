function Switching_lists = create_switching_lists_SOC(x, Params, u, q, P_diff, Method)
%% =====================================================================
%  This function divides the total population into subpopulations that will
%  be used by the controller to make switching decision.
%  It also includes  the ranking algorithm for deciding which water
%  heaters should be switched/blocked preferentially.
%
%  Martin Pfeiffer, June 2011
%
%  Modified by Evangelos Vrettos
%  April 2012
%
%  Note: many variable names include the word "boiler", which is the German
%  word for "water heater".
%
%  The selection of the subpopulation for control is based on Te1-Tmin1.
%  The ranking algorithm is based on the EWH SOC.
%
%  Modified by Evangelos Vrettos
%  October 2014
%  =====================================================================

%% List of all water heaters
all_boilers = (1:Params.n_app)';

%% ON/OFF lists
% List of water heaters that are on
boilers_ON = find(u(1,:,:)==1);

% List of water heaters that are off
boilers_OFF = find(u(1,:,:)==0);

%% Temperature lists
if strcmp(Method,'IAFSF')
    % Get temperatures at relevant location
    T_meas1 = x(Params.sensor1_location,:,:);
    
    % List of water heaters above deadband
    boilers_above_t_max = find( T_meas1 > Params.T_max1(:,:,:));
    
    % List of water heaters below deadband
    boilers_below_t_min = find( T_meas1 < Params.T_min1(:,:,:));
    
    % List of water heaters inside deadband
    if length(boilers_above_t_max)==1 && length(boilers_below_t_min)==1
        boilers_inside_deadband = setdiff(all_boilers,(union(boilers_above_t_max,boilers_below_t_min))');
    else
        boilers_inside_deadband = setdiff(all_boilers,union(boilers_above_t_max,boilers_below_t_min));
    end
    
    % List of water heaters that are a certain percentage below upper upper deadband. Used to
    % avoid immediate active re-switch ON after internal switch OFF. (Set to zero by default)
    boilers_x_percent_below_t_max = find( T_meas1 < (Params.T_max1(:,:,:) - 0.0 * (Params.T_max1(:,:,:)-Params.T_min1(:,:,:))));
    
    % List of water heaters that are a certain percentage above lower deadband. Used to
    % avoid immediate active re-switch OFF after internal switch ON. (Set to zero by default)
    boilers_x_percent_above_t_min = find( T_meas1 > (Params.T_min1(:,:,:) + 0.0 * (Params.T_max1(:,:,:)-Params.T_min1(:,:,:))));
elseif strcmp(Method,'IAPSF')
    % List of water heaters above deadband
    boilers_above_t_max = find(q(1,:,:)==-1);
    
    % List of water heaters below deadband
    boilers_below_t_min = find(q(1,:,:)==1);
    
    % List of water heaters inside deadband
    boilers_inside_deadband(:,1) = find(q(1,:,:)==0);
end

%% Lists of water heaters requiring internal control
% Requiring internal switch ON
Switching_lists.internal_ON = intersect(boilers_below_t_min, boilers_OFF);
if isempty(Switching_lists.internal_ON) == 0,
    % preallocating for speed
    power_of_internal_ON = nan(1,length(Switching_lists.internal_ON));
    for i = 1:length(Switching_lists.internal_ON),
        power_of_internal_ON(i) = Params.P1_el(:,:,Switching_lists.internal_ON(i));
    end
    power_of_internal_ON_cum = sum(sum(power_of_internal_ON));
else
    power_of_internal_ON_cum = 0;
end

% Requiring internal switch OFF
Switching_lists.internal_OFF = intersect(boilers_above_t_max, boilers_ON);
if isempty(Switching_lists.internal_OFF) == 0,
    % preallocating for speed
    power_of_internal_OFF = nan(1,length(Switching_lists.internal_OFF));
    for i = 1:length(Switching_lists.internal_OFF),
        power_of_internal_OFF(i) = Params.P1_el(:,:,Switching_lists.internal_OFF(i));
    end
    power_of_internal_OFF_cum = sum(sum(power_of_internal_OFF));
else
    power_of_internal_OFF_cum = 0;
end

%% Lists of water heaters for external control
if strcmp(Method,'IAFSF')
    % water heaters that can be actively switched without violating temperature limits
    active_switch_ON_candidates =  intersect(intersect(boilers_inside_deadband,boilers_x_percent_below_t_max), boilers_OFF);
    active_switch_OFF_candidates = intersect(intersect(boilers_inside_deadband,boilers_x_percent_above_t_min), boilers_ON);
elseif strcmp(Method,'IAPSF')
    active_switch_ON_candidates =  intersect(boilers_inside_deadband, boilers_OFF);
    active_switch_OFF_candidates = intersect(boilers_inside_deadband, boilers_ON);
end

% Initializing variables that can indicate whether switching limits have been reached
Switching_lists.ON_limit_reached = 0;
Switching_lists.OFF_limit_reached = 0;

%% Power difference that needs to be achieved through controller action
P_diff_effective = P_diff + power_of_internal_OFF_cum - power_of_internal_ON_cum;

%% Create active switch ON list
if P_diff_effective > 0 && isempty(active_switch_ON_candidates) == 0,
    % preallocating variables for speed
    SOC_1 = nan(1,length(active_switch_ON_candidates));
    power_of_active_switch_ON_candidates = nan(1,length(active_switch_ON_candidates));
    
    for i = 1:length(active_switch_ON_candidates),
        % lowest SOC = coldest heater = high switch on priority
        idx = active_switch_ON_candidates(i);
        SOC_1(i) = SOC_calc(Params, idx, x);
                
        % find power of each water heater in the unsorted priority list
        power_of_active_switch_ON_candidates(i) = Params.P1_el(:,:,active_switch_ON_candidates(i));
    end
    
    %sort active switch ON candidates according to SOC
    sorted = sortrows([SOC_1', active_switch_ON_candidates, power_of_active_switch_ON_candidates'],1);
    
    %the second column of "sorted" (in REGULAR order) is the water heater priority list
    active_switch_ON_candidates_ranked = sorted(:,2);
    
    %the third column of "sorted" (in regular order) are the corresponding water heater powers
    power_of_active_switch_ON_candidates_ranked = sorted(:,3);
    
    % create vector with commulative power of water heaters in priority list
    power_of_active_switch_ON_candidates_cum = cumsum(power_of_active_switch_ON_candidates_ranked);
    
    % adding the option of zero water heaters and zero power
    power_of_active_switch_ON_candidates_cum = [0;power_of_active_switch_ON_candidates_cum(:)];
    
    % find number of water heaters that need to be switched ON
    [~, number_of_boilers] = min(abs(power_of_active_switch_ON_candidates_cum - P_diff_effective));
    
    % subtracting 1 one from number of water heaters since the index that the "min" function produces is offset by 1
    number_of_boilers = number_of_boilers - 1;
    
    if number_of_boilers == 0,
        Switching_lists.active_ON = [];
    else
        Switching_lists.active_ON = active_switch_ON_candidates_ranked(1:(number_of_boilers));
    end
elseif P_diff_effective > 0 && isempty(active_switch_ON_candidates) == 1,
    Switching_lists.active_ON = [];
    Switching_lists.ON_limit_reached = 1;
else
    Switching_lists.active_ON = [];
end

%% Create active switch OFF list
if P_diff_effective < 0 && isempty(active_switch_OFF_candidates) == 0,
    % preallocating variables for speed
    SOC_2 = nan(1,length(active_switch_OFF_candidates));
    power_of_active_switch_OFF_candidates = nan(1,length(active_switch_OFF_candidates));
    
    for i = 1:length(active_switch_OFF_candidates),
        % highest SOC = hottest heater = high switch off priority
        idx = active_switch_OFF_candidates(i);
        SOC_2(i) = SOC_calc(Params, idx, x);
                
        %find power of each water heater in the unsorted priority list
        power_of_active_switch_OFF_candidates(i) = Params.P1_el(:,:,active_switch_OFF_candidates(i));
    end
    
    %sort active switch OFF candidates according to SOC
    sorted2 = sortrows([SOC_2', active_switch_OFF_candidates, power_of_active_switch_OFF_candidates'],1);
    
    %the second column of "sorted2" (in REVERSE order) is the water heater priority list
    active_switch_OFF_candidates_ranked = flipud(sorted2(:,2));
    
    %the third column of "sorted2" (in REVERSE order) are the corresponding water heater powers
    power_of_active_switch_OFF_candidates_ranked = flipud(sorted2(:,3));
    
    % create vector with commulative power of water heaters in priority list
    power_of_active_switch_OFF_candidates_cum = cumsum(power_of_active_switch_OFF_candidates_ranked);
    power_of_active_switch_OFF_candidates_cum = [0;power_of_active_switch_OFF_candidates_cum(:)];
    
    % find number of water heaters that need to be switched ON
    [~, number_of_boilers] = min(abs(P_diff_effective + (power_of_active_switch_OFF_candidates_cum)));
    
    % subtracting 1 one from number of water heaters since the index that the "min" function produces is offset by 1
    number_of_boilers = number_of_boilers - 1;
    
    if number_of_boilers == 0,
        Switching_lists.active_OFF = [];
    else
        Switching_lists.active_OFF = active_switch_OFF_candidates_ranked(1:(number_of_boilers));
    end
elseif P_diff_effective < 0 && isempty(active_switch_OFF_candidates) == 1,
    Switching_lists.active_OFF = [];
    Switching_lists.OFF_limit_reached = 1;
else
    Switching_lists.active_OFF = [];
end

%% Creating additional lists for evaluation purposes
% The following lists are not used to give switching signals but are used for
% determing the E_th_rel of groups of certain groups of water heaters
Switching_lists.boilers_ON = boilers_ON;
Switching_lists.boilers_OFF = boilers_OFF;
Switching_lists.active_switch_ON_candidates = active_switch_ON_candidates;
Switching_lists.active_switch_OFF_candidates = active_switch_OFF_candidates;

end
