function u = determine_heating_element_state_population(u, t, x, Params,i)
%% =======================================================================
%  Determine the state of the heating element based on the GUI settings
%  Stephan Koch
%  October 26, 2010
%
%  Modified by Evangelos Vrettos
%  October 2014
%  =======================================================================

% get the measured temperature
T_meas1 = x(Params.sensor1_location); 
T_meas2 = x(Params.sensor2_location);

%% Determine the heating element states based on the deadbands,
%  distinguish whether element 1 or 2 should be dominant
if Params.dominant_element == 1,
    % toggle the state of the first heating element if necessary
    if T_meas1 > Params.T_max1(:,:,i),
        u(1) = 0;
    elseif T_meas1 < Params.T_min1(:,:,i),
        u(1) = 1;
    end
    
    % override if it does not exist
    if Params.element1_present == 0,
        u(1) = 0;
    end
    
    % check if the first element is on, in this case block second
    % if off, operate second element according to deadband
    if u(1) == 1,
        % block the non-dominant element
        u(2) = 0;
    elseif u(1) == 0,
        % operate the non-dominant element according to deadband
        if T_meas2 > Params.T_max2(:,:,i),
            u(2) = 0;
        elseif T_meas2 < Params.T_min2(:,:,i),
            u(2) = 1;
        end
    end
    
elseif Params.dominant_element == 2,
    % toggle the state of the second heating element if necessary
    if T_meas2 > Params.T_max2(:,:,i),
        u(2) = 0;
    elseif T_meas2 < Params.T_min2(:,:,i),
        u(2) = 1;
    end
    
    % override if it does not exist
    if Params.element2_present == 0,
        u(2) = 0;
    end
    
    % check if the second element is on, in this case block first
    % if off, operate lower element according to deadband
    if u(2) == 1,
        % block the non-dominant element
        u(1) = 0;
    elseif u(2) == 0,
        % operate lower according to deadband
        if T_meas1 > Params.T_max1(:,:,i),
            u(1) = 0;
        elseif T_meas1 < Params.T_min1(:,:,i),
            u(1) = 1;
        end
    end
end

%% Override u to be equal to zero if within blocking hour
time_of_day = mod(t/3600, 24);

begin_time1 = Params.blocking_hours1(1);
end_time1 = Params.blocking_hours1(2);
if time_of_day >= begin_time1 && time_of_day <= end_time1,
    u = [0; 0];
end

begin_time2 = Params.blocking_hours2(1);
end_time2 = Params.blocking_hours2(2);
if time_of_day >= begin_time2 && time_of_day <= end_time2,
    u = [0; 0];
end


%% Override u to be equal to zero if the element does not exist
if Params.element1_present == 0,
    u(1) = 0;
end

if Params.element2_present == 0,
    u(2) = 0;
end
