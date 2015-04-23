%% ========================================================================
%  function Loads_N_ON
%  by Adrian Martinez Gomez
%  March 2015
%
%  Purpose:
%  * Extract the subpopulations of N_ON from the N_E = 20 experiments/files
%  * Plot the populations (first std then max/min)
%
%  ========================================================================
function Loads_ON(bCalculate,bNormalize);
%% Either calculate or load the N_ON subpopulation for all the N_E=20 exp
%
if(bCalculate)
    folder_dir = dir('archived_data/baseline/1000EWH/res_same/*.mat');
    
    % gather the data
    %
    N_ON = cell(1,length(folder_dir));
    for ii = 1:length(folder_dir)
        load(strcat('archived_data/baseline/1000EWH/res_same/',folder_dir(ii).name),'Results_comparison');
        N_ON{ii} = sum(squeeze(Results_comparison.urec(1,:,:)).').';
    end
else
    load(['ADRIAN/Baseline_Nloads/N_ON_populations.mat']);
    n_app = 1000;
end

%% Convert into matrix; take mean and standard deviation
mN_ON = cell2mat(N_ON);
meanN_ON = mean(mN_ON.');
stdN_ON = std(mN_ON.');

lowN_ON = meanN_ON - stdN_ON;
highN_ON = meanN_ON + stdN_ON;
maxN_ON = max(mN_ON.');
minN_ON = min(mN_ON.');

% if selected, normalize the max/min and mean with the number of EWHs n_app
% (and convert to [%] !)
if(bNormalize)
    meanN_ON = meanN_ON./n_app .* 100;
    lowN_ON = lowN_ON./n_app .* 100;
    highN_ON = highN_ON./n_app .* 100;
    maxN_ON = maxN_ON./n_app .* 100;
    minN_ON = minN_ON./n_app .* 100;
end

%% Plotting
% time vector
xx = linspace(0,24,8641);

% plot with mean and standard deviations
figure;
title(['Plot of the N_{ON} set for N_{E}=',num2str(size(mN_ON,2)),...
    ' experiments'],'FontSize',12);
hold on;
fill([xx,fliplr(xx)],[highN_ON,fliplr(lowN_ON)],...
    [0.7,0.7,0.7],'linestyle','none');
plot(xx,meanN_ON,'LineWidth',2);
hold off;
grid on;
xlim([0,24]);
xlabel('Time [h]','FontSize',12);
ylabel('N_{ON} [%]','FontSize',12);
h_legend = legend('Std','Average','Location','Best');
set(h_legend,'FontSize',12);

% plot with mean and max/mins
figure;
title(['Plot of the N_{ON} set for N_{E}=',num2str(size(mN_ON,2)),...
    ' experiments'],'FontSize',12);
hold on;
fill([xx,fliplr(xx)],[maxN_ON,fliplr(minN_ON)],...
    [0.7,0.7,0.7],'linestyle','none');
plot(xx,meanN_ON,'LineWidth',2);
hold off;
grid on;
xlim([0,24]);
xlabel('Time [h]','FontSize',12);
ylabel('N_{ON} [%]','FontSize',12);
h_legend = legend('Min/Max','Average','Location','Best');
set(h_legend,'FontSize',12);

end