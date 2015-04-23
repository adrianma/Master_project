%% ===================================================================
%  function SetPointVariationCL
%  by Adrian Martinez Gomez
%  January 2015
%   
%  Purpose:
%   * apply the optimal control trajectory (for SPV) computed through the
%   QP optimization problem (CLOSED LOOP).
%  ===================================================================

function [unew,eta] = SetPointVariationCL(Params,step,t,xrec,...
    urec,u_external_optimal);
%% External controller (temperature set-point variation)        

if(step < length(u_external_optimal))
    eta = u_external_optimal(step);
else
    eta = 0;
end

% change the set-point
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

%% Internal controllers
unew = NoControl(Params, step, t, xrec, urec);

%% Update the max min setpoint temperatures after they are used in the
% event
Params.T_min1 = Params.T_set - 1/2*Params.T_dead;
Params.T_max1 = Params.T_set + 1/2*Params.T_dead;

Params.T_min2 = Params.T_set - 1/2*Params.T_dead;
Params.T_max2 = Params.T_set + 1/2*Params.T_dead;

end