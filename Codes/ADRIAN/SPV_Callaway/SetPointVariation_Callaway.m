%% ===================================================================
%  function SetPointVariation
%  by Adrian Martinez Gomez
%  October 2014
%
%  Purpose:
%   *call the external controller for probabilistic switching
%   *external controller is used if it is in deadband
%  ===================================================================

function [unew,eta] = SetPointVariation_Callaway(Params,step,t,xrec,...
    urec);

eta = nan;

%% External controller (random switching)
%
global eta_in runs;
% global runs_events Events;

% time is within the Event interval (eather a single step unit or an entire
% time interval)
% if((t>=min(Events{runs_events})) && (t<=max(Events{runs_events})))
% draw random number between 0 and 1

eta = eta_in{runs}(step);

% perform set-point variation
%
for ii = 1:Params.n_app
    % all devices DECREASE temp. set-point by eta(%) of their deadband
    if(eta < 0)
        temp_T_set = Params.T_set-abs(eta)*Params.T_dead;
        % all devices INCREASE temp. set-point by eta(%) of their deadband
    elseif(eta >= 0)
        temp_T_set = Params.T_set+abs(eta)*Params.T_dead;
    end
    
    % [Update] Option of setting fixed deadband for all heaters
    Params.T_min1 = temp_T_set - 1/2*Params.T_dead;
    Params.T_max1 = temp_T_set + 1/2*Params.T_dead;
    
    Params.T_min2 = temp_T_set - 1/2*Params.T_dead;
    Params.T_max2 = temp_T_set + 1/2*Params.T_dead;
end
% end

%% Internal controllers
unew = NoControl(Params, step, t, xrec, urec);

%% Update the max min setpoint temperatures after they are used in the
% event
Params.T_min1 = Params.T_set - 1/2*Params.T_dead;
Params.T_max1 = Params.T_set + 1/2*Params.T_dead;

Params.T_min2 = Params.T_set - 1/2*Params.T_dead;
Params.T_max2 = Params.T_set + 1/2*Params.T_dead;

end