%% ========================================================================
%  function main_SOC
%  by Adrian Martinez Gomez
%  April 2015
%
%  Input: 
%   bDisplay: is a vector of booleans to select the plotting that we want
%   to see.
%
%  Purpose:
%   Do calculations and plotting of Carsten Heinrich's definition for the
%   State Of Charge (SOC).
%
%  ========================================================================
function main_SOC(bDisplay);
%% 1) Load and calculate Carsten's SOC
%
path_file = 'archived_data/Markov_Chain/Carsten/SOC_calculated.mat';

% if it was already calculated, then just load the file
if(exist(path_file,'file') == 2)
    load(path_file);
else
    my_SOC = cell(1,20);
    for ii = 1:20
        disp(['@@@@@@@@@@@ Loading number : ',num2str(ii),' @@@@@@@@@@@']);
        my_SOC{ii} = SOC_Carsten_3(ii);
    end
    
    % save the file to its place
    save(path_file,'my_SOC');
end

% define the time vector
xx = linspace(0,24,size(my_SOC{1},1));
% matrix for the SOC
idx_SOC = randi(20);
mSOC = my_SOC{idx_SOC};
disp(['  Index for the SOC : ',num2str(idx_SOC)]);

%% 2) Display the normalized input
if(bDisplay(1))
    figure;
    plot(xx,100.*sum(my_SOC{idx_SOC}.')./1000);
    grid on;
    xlabel('Time [h]','FontSize',12);
    ylabel('Aggregate SOC [%]','FontSize',12);
    xlim([0,24]);
end

%% 3) Run a simulation
if(bDisplay(2))
    figure;
    for ii = 1:8641
        hist(mSOC(ii,:),20);
        pause(0.001);
        grid on;
        xlabel('Number of EWHs in that bin state','FontSize',12);
        ylabel('Bin states for SOC','FontSize',12);
        xlim([0,1]);
        ylim([0,200]);
    end
end

%% 4) Plot at certain intervals
if(bDisplay(3))
    % choose every how many minutes we would like to plot the SOC bins
    intervals = input('Choose every how many minutes we want to plot the SOC bins:\n>');
    % convert from minutes to steps
    intervals = 6*intervals;
    
    count = 1;
    for ii = 1:8641
        % generate a figure every "intervals" minutes
        if(mod(ii,intervals) == 1)
            figure;
            hist(mSOC(ii,:),20);
            grid on;
            xlabel('Bin states for SOC','FontSize',12);
            ylabel('Number of EWHs in that bin state','FontSize',12);
            xlim([0,1]);
            ylim([0,250]);
            
            title(['Hour = ',num2str(floor(ii*10/3600))],'FontSize',12);
            
            if(0)
                saveas(gcf,['/Users/adrianma/Documents/ETH/Master_project/Weekly_meetings/Update_07042015/Figures/SOC_Carsten_bins/SOC_bins',num2str(count),'.fig']);
                saveas(gcf,['/Users/adrianma/Documents/ETH/Master_project/Weekly_meetings/Update_07042015/Figures/SOC_Carsten_bins/SOC_bins',num2str(count),'.png']);
                count = count + 1;
            end
        end
    end
end

fprintf(' => End of Carstens definitions! \n');

end