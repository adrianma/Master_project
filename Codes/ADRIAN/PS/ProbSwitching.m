%% ===================================================================
%  function ProbSwitching
%  by Adrian Martinez Gomez
%  October 2014
%   
%  Purpose:
%   *call the external controller for probabilistic switching
%   *external controller is used if it is in deadband
%  ===================================================================

function [unew,eta] = ProbSwitching(Params,step,t,xrec,urec);

%% Internal controllers

unew = NoControl(Params, step, t, xrec, urec);
eta = nan;

%% External controller (random switching)
%
global eta_in runs_eta runs_events Events bInput;

% to get rid off cases where bInput is not defined elsewhere
if(isempty(bInput))
    bInput = 0;
end

% events occur at at some hours 
if(bInput || ...
        ((t>=min(Events{runs_events})) && (t<=max(Events{runs_events}))))
    % draw random number between 0 and 1
    %
    if(~exist('eta_in','var'))
        eta = random('unif',-1,1);
    else
        eta = eta_in(runs_eta);
    end
    
    % perform random switching
    %
    for ii = 1:Params.n_app
        
        % get the measured temperature
        %
        % also, handle the case where we take hourly and we need to
        % initialize (to avoid errors)
        if(step-1 == 0)
            x = Params.x_init;
        else
            x = xrec(:,step-1,:);
        end
        
        T_meas1 = x(Params.sensor1_location);
        T_meas2 = x(Params.sensor2_location);
        
        % there is only a change by the probabilistic switching if 
        % temperature is between the deadband!
        bInDeadband1 = (T_meas1 <= Params.T_max1(:,:,ii))&&...
            (T_meas1 >= Params.T_min1(:,:,ii))&&Params.element1_present;
        bInDeadband2 = (T_meas2 <= Params.T_max2(:,:,ii))&&...
            (T_meas2 >= Params.T_min2(:,:,ii))&&Params.element2_present;
        bInDeadband = bInDeadband1 || bInDeadband2;
        
        if((eta >= rand) && (eta > 0) && bInDeadband)
            % turn device on
            unew(:,1,ii) = 1;
        elseif((eta <= -rand) && (eta < 0) && bInDeadband)
            % turn device off
            unew(:,1,ii) = 0;
        end
    end
end

end