function [zData,zExperimental,zValidation] = ...
    Data_Analysis_Callaway(NFolder,bSplit);
%% 1) Loading and options
%
%
folder_files = dir(['archived_data/SPV_Callaway/new_same_2/',...
    num2str(NFolder),'/*.mat']);
% load necessary stuff
load('archived_data/baseline/1000EWH/same_baseline_data/mean_water_1000EWH_same.mat');
load('archived_data/baseline/1000EWH/same_baseline_data/power_baseline_same.mat',...
    'y_baseline');
load('archived_data/baseline/1000EWH/same_baseline_data/Params_experiments_1000EWH_same.mat');

% select the options
Output_option = 2;
Exogenous_option = 2;
n_app = 1000;
Params.t_sample = 10;

%% 2)
%
zData_temp = cell(1,length(folder_files));
for jj = 1:length(folder_files)
    load(['archived_data/SPV_Callaway/new_same_2/',...
        num2str(NFolder),'/',folder_files(jj).name]);
    
    disp(['===== Loading folder : ',num2str(NFolder),' =====']);
    
    for kk = 1:N_E
        % resolve the indices
        low_idx = event_hour/Params.t_sample + 1;
        high_idx = (event_hour + 1*3600)/Params.t_sample + 1;
        
        % this is how the output is constructed, we have 3 choices:
        %
        % 1) take the aggregate power as is (NOT RECOMMENDED!)
        if(Output_option == 1)
            y{kk} = Results{kk}.Prec;
            % 2) substract the baseline corresponding to the experiment kk
        elseif(Output_option == 2)
            y{kk} = Results{kk}.Prec - SS_Prec{kk}(low_idx:high_idx);
            % 3) substract the averaged baseline
        elseif(Output_option == 3)
            y{kk} = Results{kk}.Prec - y_baseline(low_idx:high_idx);
            % 4) from Callaway's paper
        elseif(Output_option == 4)
            
            % ...(change)
            mMask = cell(1,length(Results{kk}.Prec));
            dead_active = nan(1,length(Results{kk}.Prec));
            for ii = 1:361
                % not really used
                % ...(change)
                mMask{ii} = Results{kk}.urec(1,ii,:);
                % calculate the sum of active (switch ON) temperature
                % deadbands, to scale the output "ala" Callaway
                dead_active(ii) = sum(Params_cell{kk}.T_dead(Results{kk}.urec(1,ii,:) ~= 0));
                
            end
            
            % scale by the temperature deadbands and the efficiency coeff
            y{kk} = Results{kk}.Prec./(dead_active*Params_cell{kk}.eta);
        end
        
        % normalize the output signal
        if(n_app == 500)
            y{kk} = y{kk}./sum(Params_cell{mod(counter_file+...
                NFolders*(counter_folder-1),length(Params_cell))}.P1_el);
        elseif(n_app == 1000)
            y{kk} = y{kk}./sum(squeeze(Params_cell{kk}.P1_el));
        end
        
        % get the external input AND the period ready for the iddata object
        u{kk} = eta_in{kk};
        
        period{kk} = Params.t_sample;
        
        % transpose to get the desired sizes for iddata objects
        y{kk} = y{kk}.';
        u{kk} = u{kk}.';
        
        % calculate the exogenous input signal
        %
        % (1) it's the sum of the water draws
        if(Exogenous_option == 1)
            u_exogenous{kk} = sum(squeeze(Results{kk}.mdotrec).').';
        end
        % (2) to avoid the high variance, use the averaged signal for the
        % sum (load it from file)
        if(Exogenous_option == 2)
            u_exogenous{kk} = mean_water(low_idx:high_idx).';
        end
        % (3) take the averaged steady-state aggregate power signal (aka
        % the blue line) (load it from file)
        if(Exogenous_option == 3)
            u_exogenous{kk} = mean_P_agg(low_idx:high_idx).';
        end
        
    end
    u = u(~cellfun('isempty',u));
    y = y(~cellfun('isempty',y));
    period = period(~cellfun('isempty',period));
    
    % select as 1 to consider the exogenous input
    if(Exogenous_option == 4)
        fprintf(' => No exogenous input selected\n');
    else
        u_exogenous = u_exogenous(~cellfun('isempty',u_exogenous));
        
        for kk = 1:length(u)
            u{kk} = [u{kk},u_exogenous{kk}];
        end
    end
    
    % store the data
    data = iddata(y,u,period);
    
    % postprocessing
    zData_temp{jj} = Zerod_system_ID(data,20);
    
end

%% 3)
%
% merge the experiments
zData = zData_temp{1};
for ii = 2:length(folder_files)
    zData = merge(zData,zData_temp{ii});
end

% split the data into experimental and validation
if(bSplit)
    % number of data input-output pairs
    N_Data = size(zData.ExperimentName,1);
    
    % randomize the experiments
    vRand = randperm(N_Data);
    
    % split the data accordingly
    zExperimental = getexp(zData,vRand(1:N_Data/2));
    zValidation = getexp(zData,vRand((N_Data/2 + 1):end));
    
    % ...(change)
    % [zExperimental,zValidation] = ExpVal_system_ID(zData,1);
else
    zExperimental = nan;
    zValidation = nan;
end

end