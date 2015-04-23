%% ========================================================================
%  script Chapter_2_5
%  by Adrian Martinez Gomez
%  February 2015
%
%  Modified: March 2015
%
%  Purpose:
%  Generate Figure 2.1 and 2.2 of the Thesis (probability of draw event and 
%  instance with cumulated sum)
%  ========================================================================
%% 1)
% The probability distribution from Lutz's paper.
prob_dist = [0.0129 0.0091 0.0091 0.0098 0.0144 0.0263 0.0491 0.054...
    0.054 0.0538 0.0515 0.0492 0.047 0.043 0.042 0.043 0.0492 0.061 ...
    0.0695 0.0685 0.0623 0.0528 0.0435 0.025]';

figure; 
bar(0:1:23,prob_dist.*100);
xlim([0,24]);
grid on;
xlabel('Time [h]');
ylabel('Probability of draw event [%]');

%% 2)
% Load the baseline experiment of choice
load(['archived_data/baseline/1000EWH/res/','sim_1000EWHs_runNumber',...
    num2str(3),'.mat'],'Results_comparison');
% matrix of water draws
idx = randi(1000);
yy = Results_comparison.mdotrec(1,:,idx);
% time vector
xx = linspace(0,24,size(yy,2));

figure;
subplot(2,1,1);
plot(xx,yy);
grid on;
xlabel('Time [h]');
ylabel('Flow rate [liters/min]');

subplot(2,1,2);
plot(xx,cumsum(yy));
grid on;
xlabel('Time [h]');
ylabel('Accumulated water flow [liters]');

disp([' => Size of the tank = ',num2str(Params.m(:,:,idx)),' [liters]']);