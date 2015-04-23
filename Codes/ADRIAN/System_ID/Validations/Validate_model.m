%% ========================================================================
%  function Check_model_properties
%  by Adrian Martinez Gomez
%  February 2015
%
%  Purpose:
%   * checks stability, controllability and observability for a state-space
%   system ss, and shows on the console the properties.
%  ========================================================================
function [Fit_vector] = Validate_model(ValData,ss)

% set the inital condition x0 to zero (as it is what the model was
% estimated for)
if(1)
    opt_comp = compareOptions;
    opt_comp.InitialCondition = 'z';
end

% calculate comparison fit vector 'Fit_vector' to validation data 'ValData'
[Y_response,Fit_vector,~] = compare(ValData,ss,opt_comp);
% get the fit vector ready
if(1)
    Fit_vector = cell2mat(Fit_vector);
else
    mY_response = zeros(381,length(ValData.ExperimentName));
    for ii = 1:70
        mY_response(:,ii) =  Y_response{ii}.y;
    end
    tilde_Y = cell2mat(ValData.y);
    
    % ...(change)
    % Here would be an alternative definition for Maryam for the fit signal
end

% plot the fit for every experiment
figure;
stem(Fit_vector);
grid on;
xlabel('Experiment number');
ylabel('Fit to validation data [%]');

% display mean and standard deviation for the fit for all the validation
% experiments
disp(['Average fit for validation data: ',num2str(mean(Fit_vector)),' %']);
disp(['     Standard deviation : ',num2str(std(Fit_vector)),' %']);

if(0)
    f_v_4_pos = Validate_model(val_4_pos,My_Models.my_mod_4_pos); close;
    f_v_7_pos = Validate_model(val_7_pos,My_Models.my_mod_7_pos); close;
    f_v_12_pos = Validate_model(val_12_pos,My_Models.my_mod_12_pos); close;
    f_v_16_pos = Validate_model(val_16_pos,My_Models.my_mod_16_pos); close;
    f_v_18_pos = Validate_model(val_18_pos,My_Models.my_mod_18_pos); close;
    f_v_20_pos = Validate_model(val_20_pos,My_Models.my_mod_20_pos); close;
    f_v_22_pos = Validate_model(val_22_pos,My_Models.my_mod_22_pos); close;
end

end

