%% ===================================================================
%  function Present_Data_system_ID
%  by Adrian Martinez Gomez
%  January 2015
%
%  Purpose:
%   * rearrange the zero'd my_ref_signal and exogenous_input to use it in
%   the main_CL_simulation file.
%  ===================================================================
function [my_ref_signal,exo_input] = Present_Data_system_ID(data);

NN = length(data.y);

% pre-allocate
my_ref_signal = nan(NN,8641);
exo_input = nan(NN,8641);
idx_hour = nan;

for ii = 1:NN
    % get the right index out from the iddata experimental object
    % code is unclean, but is paired with how the file Data_system_ID.m was
    % constructed
    if(1<=mod(ii,140) && mod(ii,140)<=20)
        idx_hour = 12;
    elseif(21<=mod(ii,140) && mod(ii,140)<=40)
        idx_hour = 16;
    elseif(41<=mod(ii,140) && mod(ii,140)<=60)
        idx_hour = 18;
    elseif(61<=mod(ii,140) && mod(ii,140)<=80)
        idx_hour = 20;
    elseif(81<=mod(ii,140) && mod(ii,140)<=100)
        idx_hour = 22;
    elseif(101<=mod(ii,140) && mod(ii,140)<=120)
        idx_hour = 4;
    elseif(121<=mod(ii,140) && mod(ii,140)<=140)
        idx_hour = 7;
    end
    
    % store the zero'd reference signal
    my_ref_signal(ii,:) = [zeros(1,idx_hour*360),...
        data.y{ii}.',...
        zeros(1,8640 - (idx_hour+1)*360)];
    % store the zero'd exogenous input
    exo_input(ii,:) = [zeros(1,idx_hour*360),...
        data.u{ii}(:,2).',...
        zeros(1,8640 - (idx_hour+1)*360)];
    
    % to generate the "short" (361 points = 1 hour) signals
    if(0)
        my_ref_signal(ii,:) = data.my_ref_signal(ii,(idx_hour*360):(idx_hour*360+360));
        exo_input(ii,:) = data.exo_input(ii,(idx_hour*360):(idx_hour*360+360));
    end
end

end
