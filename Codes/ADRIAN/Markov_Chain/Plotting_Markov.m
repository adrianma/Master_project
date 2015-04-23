%% ========================================================================
%  function Plotting_histogram
%  by Adrian Martinez Gomez
%  March 2015
%
%  Inputs:
%   *idx_baseline: the file to be loaded. 
%       (Already pre-computed with the "Offline_Markov" script)
%   *bSimulate: boolean; choose "1" for a simulation and "0" for sampled
%   examples (every "intervals" minutes)
%   *intervals: sampling [minutes]
%
%  Purpose:
%   Plot the number of EWHs in each bin and ON/OFF status over one day
%   measured at X minute intervals.
%
%  ========================================================================
function Plotting_Markov(idx_baseline,bSimulate,intervals);
%% Load the corresponding file
%
path_file = ['archived_data/Markov_Chain/sim_',...
        num2str(idx_baseline),'_bins.mat'];
load(path_file);
    
% if that idx_baseline 
if(~exist(path_file,'file') == 2)
    error('idx_baseline was not computed previously (run Offline_Markov)');
end

%% Do the simulation/plotting of the figures
%
% if bSimulate is 1, then do a simulation;
% else, choose some examples (sampled every "intervals" minutes)
if(bSimulate)
    
    figure;
    for ii = 1:N_sim
        subplot(1,2,1);
        hist(njs_ON{ii},N_bins);
        xlabel('Temperature deadband bins','FontSize',12);
        ylabel('Number of ON EWHs in the bin','FontSize',12);
        title(['Time step = ',num2str(ii)],'FontSize',12);
        xlim([0,N_bins + 1]);
        grid on;
        
        subplot(1,2,2);
        hist(njs_OFF{ii},N_bins);
        xlabel('Temperature deadband bins','FontSize',12);
        ylabel('Number of OFF EWHs in the bin','FontSize',12);
        xlim([0,N_bins + 1]);
        grid on;
        
        pause(0.005);
    end
    
else
    % choose every how many minutes we would like to plot the Markov state
    % (if not chosen take 30 minutes)
    if(intervals ~= 0)
        intervals = 6*intervals;
    else
        intervals = 6*30;
    end
    
    count = 1;
    for ii = 1:N_sim
        % generate a figure every "intervals" minutes
        if(mod(ii,intervals) == 1)
            
            temp_N_ON = sum(Results_comparison.urec(1,ii,:));
            temp_N_OFF = n_app - temp_N_ON;
            
            figure;
            subplot(1,2,1);
            hist(njs_ON{ii},N_bins);
            xlabel('Temperature deadband bins','FontSize',12);
            ylabel('Number of ON EWHs in the bin','FontSize',12);
            title(['Time step = ',num2str(ii),'; Hour = ',...
                num2str(ii/360)],'FontSize',12);
            h_legend = legend(['N_{ON} = ',num2str(temp_N_ON)],...
                'Location','Best');
            set(h_legend,'FontSize',12);
            xlim([0,N_bins + 1]);
            grid on;
            
            subplot(1,2,2);
            hist(njs_OFF{ii},N_bins);
            xlabel('Temperature deadband bins','FontSize',12);
            ylabel('Number of OFF EWHs in the bin','FontSize',12);
            h_legend = legend(['N_{OFF} = ',num2str(temp_N_OFF)],...
                'Location','Best');
            set(h_legend,'FontSize',12);
            xlim([0,N_bins + 1]);
            grid on;
            
            % save the figure
            saveas(gcf,['ADRIAN/Markov_Chain/my_fig',num2str(count)],...
                'png');
            count = count + 1;
        end
    end
end

fprintf(' => Script terminated!\n');

end