%% ========================================================================
%  script Chapter_2_6
%  by Adrian Martinez Gomez
%  February 2015
%
%  Modified: March 2015
%
%  Purpose:
%  Generate Figure 2.3 of the Thesis (dynamics of 1 EWH for 1 day)
%  ========================================================================
%% 1) Load the baseline experiment of choice
%
% load the 1000 EWH experiments and select an EWH from the population
if(1)
    % index of the file we want to load
    idx_file = 3;
    disp(['===== Loading file number : ',num2str(idx_file),' =====']);
    
    load(['archived_data/baseline/1000EWH/res_same/sim_1000EWHs_runNumber',...
        num2str(idx_file),'.mat']);
    
    % Choose one of the EWHs
    if(0)
        idx_EWH = randi(1000);
    else
        idx_EWH = 10;
    end
    disp(['   With EWH number : ',num2str(idx_EWH)]);
end
% select an EWH from the individual experiments (which have different
% initial temperatures T_0)
if(0)
    idx_file = 2;
    disp(['===== Loading file number : ',num2str(idx_file),' =====']);
    
    load(['archived_data/baseline/different_T0_1EWH/sim_1EWHs_runNumber',...
        num2str(idx_file),'.mat']);
    
    % Choose one of the EWHs
    idx_EWH = 1;
    disp(['   With EWH number : ',num2str(idx_EWH)]);
end
% load the experiments WITHOUT water draws
if(0)
    % idx_file = 2;
    disp(['===== Loading file number : ',num2str(idx_file),' =====']);
    
    load(['archived_data/baseline/different_T0_1EWH/no_draws/sim_1EWHs_runNumber',...
        num2str(idx_file),'.mat']);
    
    % Choose one of the EWHs
    idx_EWH = 1;
    disp(['   With EWH number : ',num2str(idx_EWH)]);
end

%% 2) Get all the needed signals
%
% display the initial temperature vector
disp('    Initial temperature vector T_{0}:')
Params.x_init(:,:,idx_EWH)

% vector of water draws for that given tank (idx)
zz11 = Results_comparison.mdotrec(1,:,idx_EWH).*60; 
% matrix of the temperature (layer-wise) within the tank 
zz22 = flip(Results_comparison.xrec(2:11,:,idx_EWH));
% vector with the switch state 
zz33 = Results_comparison.urec(1,:,idx_EWH);
% vector with the temperature at the sensor
zz44a = Results_comparison.xrec(Params.sensor1_location,:,idx_EWH);
zz44b = Results_comparison.xrec(Params.element1_location,:,idx_EWH);
% time vector
xx = linspace(0,24,size(zz11,2));

%% 3) Plot the water draws for the day
% 
% number of subplots
nSubPlots = 4;

figure;

ah1 = subplot(nSubPlots,1,1);
plot(xx,zz11);
grid on;
ylabel('Flow rate [liters/min]','FontSize',12);
title(['T_{set} = ',num2str(Params.T_set(:,:,idx_EWH)),...
    '[Celsius] and \delta = ',num2str(Params.T_dead(:,:,idx_EWH)),...
    '[Celsius]'],'FontSize',12);
xlim([xx(1),xx(end)]);

ah2 = subplot(nSubPlots,1,2);
% set colormap
colormap(jet(100));   
% draw image and scale colormap to values range
imagesc(zz22); 
colorbar;
ax = gca;
ax.YLabel.String = 'Y Axis';
ax.YLabel.FontSize = 12;
ylabel('Layer number of the EWH','FontSize',12);

ylabels = 10:-1:1;
xlabels = [0,4,8,12,16,20,24];
% Change x-axis ticks
set(gca,'YTick',1:length(ylabels),'XTick',1:length(xlabels)); 
% Change x-axis ticks labels
set(gca,'YTickLabel',ylabels,'XTickLabel',xlabels); 

ah3 = subplot(nSubPlots,1,3);
area(xx,zz33);
grid on;
ylim([0,2]);
ylabel('Binary switch state','FontSize',12);
% only put the time label if we are plotting 3 subplots, else it comes with
% subplot number 4
if(nSubPlots == 3)
    xlabel('Time [h]','FontSize',12);
end
legend('ON','Location','Best','FontSize',12);

% find current position [x,y,width,height]
pos1 = get(ah1,'Position');
pos2 = get(ah2,'Position');
pos3 = get(ah3,'Position');

% set width of first/third axes equal to second
pos1(3) = pos2(3);
pos3(3) = pos2(3);
set(ah1,'Position',pos1);
set(ah3,'Position',pos3);

if(nSubPlots == 4)
    high_temp = ones(1,length(zz44a)).*...
        (Params.T_set(:,:,idx_EWH) + 0.5.*Params.T_dead(:,:,idx_EWH));
    low_temp = ones(1,length(zz44b)).*...
        (Params.T_set(:,:,idx_EWH) - 0.5.*Params.T_dead(:,:,idx_EWH));
    
    ah4 = subplot(nSubPlots,1,4);
    hold on;
    plot(xx,zz44a);
    plot(xx,zz44b,'r');
    % also plot the deadband
    plot(xx,low_temp,'--k');
    plot(xx,high_temp,'--k');    
    hold off;
    grid on;
    legend(['T_{sensor}=',num2str(Params.sensor1_location - 1)],...
        ['T_{heat}=',num2str(Params.element1_location - 1)],'T_{min}',...
        'T_{max}','Location','bestoutside');
    xlim([xx(1),xx(end)]);
    ylabel('Temperature [Celsius]','FontSize',12);
    xlabel('Time [h]','FontSize',12);
    
    % find current position [x,y,width,height]
    pos4 = get(ah4,'Position');
    % set width of first/third axes equal to second
    pos4(3) = pos2(3);
    set(ah4,'Position',pos4);
end