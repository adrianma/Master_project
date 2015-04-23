%% ========================================================================
%  script Chapter_3_2
%  by Adrian Martinez Gomez
%  April 2015
%
%  Purpose:
%  Generate Figure 3.1 of the Thesis (population baseline for 1 day)
%  ========================================================================
idx = 8;
load(['archived_data/baseline/1000EWH/res_same/sim_1000EWHs_runNumber',...
    num2str(idx),'.mat']);

% define the time vector
xx = linspace(0,24,8641);
% define the aggregate power vector
yy = Results_comparison.Prec;
% get the water draws
zz = sum(WaterDrawScenarioReal.flow_rates);
% adjust from kg/hour to kg/step, where period is 10 seconds
zz = zz./(60*10);

% do the plot
figure;

subplot(2,1,1);
bar(xx,zz);
grid on;
xlim([0,24]);
xlabel('Time [h]');
ylabel('Water draws volume [liters/step]');
legend('Sampling time = 10 seconds','Location','Best');

subplot(2,1,2);
plot(xx,yy);
grid on;
xlim([0,24]);
xlabel('Time [h]');
ylabel('Aggregate Power [kW]');
legend('N_{app} = 1000','Location','Best');