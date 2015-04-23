%% ========================================================================
%  function Offline_histogram
%  by Adrian Martinez Gomez
%  March 2015
%
%  Purpose:
%   Do a Markov bin analysis like described in Maryam's E-mail.
%   1) Divide the deadband for each EWH into N_bin/2 (= 50 for example)
%   2) Then have N_bin (= 100) to represent both normalized temperature
%   from deadband and the ON/OFF state.
%
%  ========================================================================
function Offline_Markov(idx_baseline);

path_file = ['archived_data/Markov_Chain/sim_',...
        num2str(idx_baseline),'_bins.mat'];
% if this idx_baseline was not already computed
if(~(exist(path_file,'file') == 2))
    
    % load the file required
    disp(['===== Loading file number : ',num2str(idx_baseline),' =====']);
    
    load(['archived_data/baseline/1000EWH/res_same/sim_1000EWHs_runNumber',...
        num2str(idx_baseline),'.mat']);
    
    % number of bins
    N_bins = 20;
    % length of time vector simulation
    N_sim = length(Results_comparison.Prec);
    
    % predefinitions
    Temperatures = cell(1,N_sim);
    ON_OFF_state = cell(1,N_sim);
    vvv = cell(1,N_sim);
    www = cell(1,N_sim);
    
    %% Core loop(s)
    for ii = 1:N_sim
        
        % 1.1) get the temperature evolution at the disk of the sensor
        %      for each EWH
        Temperatures{ii} = squeeze(Results_comparison.xrec(...
            Params.sensor1_location,ii,:));
        
        % 1.2) ON/OFF state of the internal switch for an entire day
        %      for each EWH
        ON_OFF_state{ii} = Results_comparison.urec(1,ii,:);
        
        % 1.3) iterate through all EWHS of the population
        %
        vvv{ii} = [];
        www{ii} = [];
        for jj = 1:n_app
            % for each EWH get the limits for the bins of the deadband
            Edges_deadband = Params.T_set(:,:,jj) + ...
                linspace(-Params.T_dead(:,:,jj)./2,Params.T_dead(:,:,jj)./2,N_bins + 1);
            
            % bins lower edge
            lower_limits = Edges_deadband(1:end-1);
            % bins upper edge
            upper_limits = Edges_deadband(2:end);
            % bins center
            bin_centers = (lower_limits + upper_limits)./2;
            
            % assign values to bins
            [~,binIdx] = histc(Temperatures{ii}(jj),[lower_limits upper_limits(end)]);
            % store the temperature value
            vvv{ii} = [vvv{ii},binIdx];
            % store also if the given indeces are ON/OFF
            www{ii} = [www{ii},ON_OFF_state{ii}(jj)]; 
        end
        
        
        disp([' %%%% RUN NUMBER : ',num2str(ii),' %%%%']);
    
    end
    
    % 1.4) convert to matrices
    ON_OFF_mat = cell2mat(ON_OFF_state.');
    vvv_mat = cell2mat(vvv.');
    www_mat = cell2mat(www.');
    
    % get the sub-populations for ON and OFF matrices
    njs_ON = cell(1,N_sim);
    njs_OFF = cell(1,N_sim);
    for ii = 1:N_sim
        njs_ON{ii} = vvv_mat(ii,www_mat(ii,:) == 1);
        njs_OFF{ii} = vvv_mat(ii,www_mat(ii,:) == 0);
    end
    
    % 1.5) save the data
    %
    save(['archived_data/Markov_Chain/sim_',num2str(idx_baseline),...
        '_bins.mat']);
    
% if it was already computed, simply tell the user and terminate
else
    % load the data
    disp(' Change the idx_baseline value; this one was already computed!');
end

fprintf(' => Script terminated!\n');

end