%% ===================================================================
%  function Data_system_ID
%  by Adrian Martinez Gomez
%  November 2014
%
%  Last modified:
%  December 2014
%
%  Inputs:
%   * Method: string with the method selected (prob. switching or set-point
%   variation.
%   * NFolders: the number of folder from which we want to read files
%	  (starting with 1,2,...,NFolders)
%	* n_app: number of loads for the population
%   * Exogenous_option: number giving options for the different
%     exogenous input options (1:water draw sum, 2:averaged water draw sum,
%     3:averaged steady-state aggregate power, 4:no exogenous input).
%   * Output_option: number giving options for the different output options
%       1: take the aggregate power as is (NOT RECOMMENDED!)
%       2: substract the baseline corresponding to the experiment kk
%       3: substract the averaged baseline
%   * bSame: boolean; select 1 if population P is the same (RECOMMENDED!)
%
%  Purpose:
%   * compute the iddata object experiments out from the Amplitudes
%     experiments out from 'main_experiments'.
%  ===================================================================
function Data = Data_system_ID(Method,NFolders,n_app,Exogenous_option,...
    Output_option,bSame);
% load the Params cell-array for normalization of the output signal
if(n_app == 500)
    load('VAGGELIS/saved_data/baseline/Params_experiments_500EWH.mat');
elseif(n_app == 1000)
    if(bSame)
        load('archived_data/baseline/1000EWH/same_baseline_data/Params_experiments_1000EWH_same.mat');
    else
        load('archived_data/baseline/1000EWH/Params_experiments_1000EWH.mat');
    end
else
    error('That number of loads was not considered! Try again');
end

% check if bExogenous has the allowed form for the script
% ...(change) remove this piece
if(0)
if(sum(bExogenous)~=1)
	error('Wrong value for bExogenous. Only one entry can be 1!');
end
end

% initialize the output counter
output_counter = 1;

% choose the folder for the method (PS:prob. switching,
% SPV:set-point variation)
if(strcmp(Method,'ProbSwitching'))
    tFolderMethod = 'PS';
elseif(strcmp(Method,'SetPointVariation'))
    tFolderMethod = 'SPV';
end

for jj = NFolders(1):NFolders(2)
    disp(['%%%%%%% Loading folder number ',num2str(jj),' %%%%%%%']);
    if(bSame)
        temp_path = ['archived_data/',tFolderMethod,'/',num2str(n_app),'EWH/new_same_4/',num2str(jj)];
    else
        temp_path = ['archived_data/',tFolderMethod,'/',num2str(n_app),'EWH/',num2str(jj)];
    end
    
    folder_files = dir([temp_path,'/*.mat']);
    for ii = 1:length(folder_files)
        % select the according path for the file
        if(n_app == 500)
            path_file = strcat('VAGGELIS/saved_data/',tFolderMethod,'/',...
                num2str(jj));
            path_file = strcat(path_file,'/');
            path_file = strcat(path_file,folder_files(ii).name);
        elseif(n_app == 1000)
            path_file = strcat(temp_path,'/');
            path_file = strcat(path_file,folder_files(ii).name);
        end
        
        [data] = params_get(path_file,Params_cell,ii,jj,...
            NFolders(2),n_app,Exogenous_option,Output_option,bSame);
        
        % store the values calculated into the Output struct
        DataOut{output_counter} = data;
        
        % update the counter number
        output_counter = output_counter + 1;
        
        clearvars -except folder_files ii jj Params_cell tFolderMethod ...
            NFolders orders orders_count fits fits_count DataOut...
            output_counter n_app temp_path Exogenous_option Output_option bSame;
    end
end

% store the data in a single iddata element
Data = DataOut{1};
for ii = 2:size(DataOut,2)
    Data = merge(Data,DataOut{ii});
end

end

function [data] = params_get(path_file,Params_cell,counter_file,...
    counter_folder,NFolders,n_app,Exogenous_option,Output_option,bSame);
% load the files (path given as argument)
load(path_file);
% this one is only needed if variant (2) for exogenous input is selected
if(n_app == 500)
    load('archived_data/baseline/500EWH/mean_water_P_agg_500EWH.mat');
end
if(n_app == 1000)
    if(bSame)
        load('archived_data/baseline/1000EWH/same_baseline_data/mean_water_1000EWH_same',...
            'mean_water');
        load('archived_data/baseline/1000EWH/same_baseline_data/power_baseline_same.mat',...
            'y_baseline');
    else
        load('archived_data/baseline/1000EWH/mean_water_P_agg_1000EWH.mat');
        temp = mean_water.signal; clear meant_water;
        mean_water = temp; clear temp;
        load('archived_data/baseline/1000EWH/power_baseline.mat',...
            'y_baseline');
    end
    
end
% display the amplitude range in the console
if(counter_file==1)
    Amplitudes
end

% initialize/preallocate
data = nan;

%%
vIdx = find(~isnan(eta_in));
Params.t_sample = 10;

% preallocate the cell-arrays
u = cell(1,N_E);
u_exogenous = cell(1,N_E);
y = cell(1,N_E);
period = cell(1,N_E);

for kk = 1:N_E
    if(sum(ismember(vIdx,kk))==1)
        
        % resolve the indices
        low_idx = event_hour/Params.t_sample + 1;
        high_idx = (event_hour + 1*3600)/Params.t_sample + 1;
        low_event_idx = Events{kk}(1)/Params.t_sample + 1 - low_idx;
        high_event_idx = Events{kk}(end)/Params.t_sample + 1 - low_idx;
        
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
        u{kk} = [Results{kk}.etas(1).*ones(1,high_event_idx - low_event_idx + 1),zeros(1,(high_idx-low_idx)- high_event_idx)];
        
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
        
    else
        error('Error with the input amplitudes!');
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

data = iddata(y,u,period);
if(0)
    data = detrend(data);
end

%% BIN 1
%
% convert the experiment names, to avoid problems when merging
if(0)
    
    for ii = 1:size(Data_2.ExperimentName,1)
        temp = str2num(Data_2.ExperimentName{ii}(4:end));
        Data_2.ExperimentName{ii} = strcat('Exp',num2str(temp+size(Data_2.ExperimentName,1)));
    end
    clear ii temp;
    
end

end