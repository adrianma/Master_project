function WaterDrawScenario = construct_water_draw_scenario(Params, n_hh, daily_profile, n_draw_dist, n_draw_dist_param, shower_dist, shower_dist_param, mediumLoad_dist, mediumLoad_dist_param, shortLoad_dist, shortLoad_dist_param, flow_rate_shower_dist, flow_rate_shower_dist_param, flow_rate_mediumLoad_dist, flow_rate_mediumLoad_dist_param, flow_rate_shortLoad_dist, flow_rate_shortLoad_dist_param, t_sample, t_sim)
%% ========================================================================
%  Import the water draw profiles
%  Stephan Koch
%  EEH - Power Systems Laborator, ETH Zurich
%  September 14, 2010
%
%  Modified by Evangelos Vrettos
%  October 2014
%  ========================================================================

%% Create the number of draws
n_draw = nan(1,n_hh);

% Read the upper and lower limits for number of draws
ll = nan(1,n_hh);
ul = nan(1,n_hh);
for i = 1:1:n_hh
    switch Params.m(:,:,i)
        case Params.m_categories(1)
            ll(i) = n_draw_dist_param(1,1);
            ul(i) = n_draw_dist_param(1,2);
        case Params.m_categories(2)
            ll(i) = n_draw_dist_param(2,1);
            ul(i) = n_draw_dist_param(2,2);
        case Params.m_categories(3)
            ll(i) = n_draw_dist_param(3,1);
            ul(i) = n_draw_dist_param(3,2);
        case Params.m_categories(4)
            ll(i) = n_draw_dist_param(4,1);
            ul(i) = n_draw_dist_param(4,2);
    end
    n_draw(i) = random(n_draw_dist, ll(i), ul(i), 1, 1);
end
n_draw(n_draw < 1) = 1;

%% Create the average number of events each hour
average_draw_events_per_hour = daily_profile * n_draw;

%% Create the actual event times (draw starting times)
draw_events = cell(1, n_hh);
max_draw_events = 1;
for i = 1:1:n_hh, % calculate this per appliance
    
    % Repeat until the number of draws is between the lower and upper limits
    while (length(draw_events{i}) < ll(i)) || (length(draw_events{i}) > ul(i))
        draw_events{i} = [];
     
        current_av_draw_events_per_hour = average_draw_events_per_hour(:,i);
        t = 0.001;
        current_av_draw_events = current_av_draw_events_per_hour(ceil(t));
        % progress in time with exponential distribution (draw events arriving exponentially)
        t = t + random('exp', 1/current_av_draw_events);
        while t < t_sim
            % log draw event
            draw_events{i} = [draw_events{i}, t];
            % assemble the times of the draw events with rate of that hour
            current_av_draw_events = current_av_draw_events_per_hour(ceil(t));
            % progress in time with exponential distribution (draw events arriving exponentially)
            t = t + random('exp', 1/current_av_draw_events);
            
        end
        max_draw_events = max(max_draw_events, length(draw_events{i}));
        if mod(i,100) == 0,
            disp(['Household ', num2str(i), ': draw events at ', num2str(draw_events{i}), ' ...']);
        end
    end
end

%% Decide whether a water draw is long or short and assign duration
% assumption: when more water is being drawn, long draws are more likely

draw_event_types = cell(n_hh, length(max_draw_events));
draw_event_duration = cell(n_hh, length(max_draw_events));
flow_rate =  cell(n_hh, length(max_draw_events));

% Define the peak period, during any water draw is either 'shower' or
% 'mediumLoad' with probability 70%, and 'shortLoad' with probability 30%
peakHrs = (7:23);
probPeak = 0.7;

% During the peak period, and if a water draw is either 'shower' or 'mediumLoad', 
% then the water draw is 'shower' with probability 14%, and 'mediumLoad' with probability 86%
probShower = 0.14;

for i = 1:1:n_hh, % calculate this per appliance
    for j = 1:1:length(draw_events{i}), % and per event

        current_hour = ceil(draw_events{i}(j));        
   
        if ismember(current_hour,peakHrs)
            if random('unif', 0, 1) < probPeak
                if random('unif', 0, 1) < probShower
                    draw_event_types{i, j} = 'shower';
                    draw_event_duration{i}(j) = random(shower_dist, shower_dist_param(1), shower_dist_param(2));
                    flow_rate{i}(j) = random(flow_rate_shower_dist, flow_rate_shower_dist_param(1), flow_rate_shower_dist_param(2));
                else
                    draw_event_types{i, j} = 'mediumLoad';
                    draw_event_duration{i}(j) = random(mediumLoad_dist, mediumLoad_dist_param(1), mediumLoad_dist_param(2));
                    flow_rate{i}(j) = random(flow_rate_mediumLoad_dist, flow_rate_mediumLoad_dist_param(1), flow_rate_mediumLoad_dist_param(2));
                end
            else
                draw_event_types{i, j} = 'shortLoad';
                draw_event_duration{i}(j) = random(shortLoad_dist, shortLoad_dist_param(1), shortLoad_dist_param(2));
                flow_rate{i}(j) = random(flow_rate_shortLoad_dist, flow_rate_shortLoad_dist_param(1), flow_rate_shortLoad_dist_param(2));
            end
        else
            if random('unif', 0, 1) < (1-probPeak)
                if random('unif', 0, 1) < probShower
                    draw_event_types{i, j} = 'shower';
                    draw_event_duration{i}(j) = random(shower_dist, shower_dist_param(1), shower_dist_param(2));
                    flow_rate{i}(j) = random(flow_rate_shower_dist, flow_rate_shower_dist_param(1), flow_rate_shower_dist_param(2));
                else
                    draw_event_types{i, j} = 'mediumLoad';
                    draw_event_duration{i}(j) = random(mediumLoad_dist, mediumLoad_dist_param(1), mediumLoad_dist_param(2));
                    flow_rate{i}(j) = random(flow_rate_mediumLoad_dist, flow_rate_mediumLoad_dist_param(1), flow_rate_mediumLoad_dist_param(2));
                end
            else
                draw_event_types{i, j} = 'shortLoad';
                draw_event_duration{i}(j) = random(shortLoad_dist, shortLoad_dist_param(1), shortLoad_dist_param(2));
                flow_rate{i}(j) = random(flow_rate_shortLoad_dist, flow_rate_shortLoad_dist_param(1), flow_rate_shortLoad_dist_param(2));
            end
        end
    end
    if mod(i,100) == 0,
        disp(['Household ', num2str(i), ': assigning draw event durations ...']);
    end
end
        
%% Sample the draw profiles according to the sampling time
% create a matrix (households x time_steps) of water draw in per unit of flow rate (60 kg/h)
temp_flow_rate_over_time = zeros(n_hh, t_sim/t_sample);

for i = 1:1:n_hh, % calculate this per appliance
    time_step = 0;
    if mod(i,10) == 0,
        disp(['Sampling the draw profile of household ', num2str(i), ' ...'])
    end
    for t = 0:t_sample:t_sim,
        time_step = time_step + 1;
        % cycle through the event list: find the most recent (potentially
        % still active) draw event
        passed_events = find(t > draw_events{i});
        if ~isempty(passed_events),
            most_recent_event = passed_events(end);
        else
            most_recent_event = [];
        end
        % if there the time went beyond an event in the list, but not
        % beyond its duration
        if t <= draw_events{i}(most_recent_event) + draw_event_duration{i}(most_recent_event),
            active_event = most_recent_event;
        else
            active_event = [];
        end
        
        if ~isempty(active_event),
            temp_flow_rate_over_time(i,time_step) = flow_rate{i}(active_event);
        else
            temp_flow_rate_over_time(i,time_step) = 0;
        end
    end
end

%% Calculate the total household consumption
% Minimum value of flow rate is zero
flow_rate_over_time = temp_flow_rate_over_time;
flow_rate_over_time(isnan(flow_rate_over_time)) = 0;
flow_rate_over_time = max(0,flow_rate_over_time);

% Normalize based on t_sample
flowRate = flow_rate_over_time*t_sample;
totCons = sum(flowRate,2);

%% Export the water draw profiles to the workspace
WaterDrawScenario.flow_rates = flow_rate_over_time;
WaterDrawScenario.time = 0:t_sample:t_sim;
WaterDrawScenario.total_draw = totCons;             % in liters per day

%% Recycle Bin
% for i = 1:1:n_hh
%     numDraws = length(draw_event_duration{i});
%     totDrawTime = sum(draw_event_duration{i});
%     avgRate = daily_consumption(i)/totDrawTime;
%     random(flow_rate_dist,avgRate,60,numDraws,1)
% 
% end
