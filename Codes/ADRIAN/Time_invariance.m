%% ========================================================================
%  script Time_invariance.m
%  by Adrian Martinez Gomez
%  March 2015
%
%  Purpose:
%  Generate the figure for Chapter 4.2 of the written Master Thesis where
%  the time variance of the system is shown with a counterexample!
%
%  ========================================================================
%% 1)
idx = 8;
load(['archived_data/baseline/1000EWH/res_same/sim_1000EWHs_runNumber',...
    num2str(idx),'.mat']);
load('archived_data/baseline/1000EWH/same_baseline_data/power_baseline_same.mat');
% store the steady-state power traj. for later analysis
SS_Prec = Results_comparison.Prec;

global eta_in Events runs_events runs_eta;

Input.amplitude = 0.3;
Input.duration = 10;

event_hour = [4,18,20].*3600;
eta_in = [Input.amplitude,Input.amplitude,Input.amplitude];

% define the simulated time (=1 day)

Params.t_sample = 10.0;
Params.duration = Input.duration*60/Params.t_sample*Params.t_sample;
% water draws required
Params.bWaterDraw = 1;

Results = cell(1,length(event_hour));
for ii = 1:length(event_hour)
    runs_events = ii;
    runs_eta = ii;
    
    Params.t_init = event_hour(ii);
    Params.t_sim = event_hour(ii) +1*3600;
    
    Events{ii} = Params.t_init:1:(Params.t_init + Params.duration);
    
    Results{ii} = simulate_population(Params,PrelimModel,...
        WaterDrawScenarioReal,WaterDrawScenarioReal,'SetPointVariation');
end

%% 2)
low_idx_1 = event_hour(1)/Params.t_sample + 1;
high_idx_1 = (event_hour(1) + 1*3600)/Params.t_sample + 1;
low_idx_2 = event_hour(2)/Params.t_sample + 1;
high_idx_2 = (event_hour(2) + 1*3600)/Params.t_sample + 1;

% outputs (normalized!)
y_1 = Results{1}.Prec./sum(squeeze(Params.P1_el)) - y_baseline_norm(low_idx_1:high_idx_1);
y_2 = Results{2}.Prec./sum(squeeze(Params.P1_el)) - y_baseline_norm(low_idx_2:high_idx_2);

low_event_idx_1 = Events{1}(1)/Params.t_sample + 1 - low_idx_1;
high_event_idx_1 = Events{1}(end)/Params.t_sample + 1 - low_idx_1;
low_event_idx_2 = Events{2}(1)/Params.t_sample + 1 - low_idx_2;
high_event_idx_2 = Events{2}(end)/Params.t_sample + 1 - low_idx_2;

% external input
u_1 = [Input.amplitude.*ones(1,high_event_idx_1 - low_event_idx_1 + 1),zeros(1,(high_idx_1-low_idx_1)- high_event_idx_1)];
u_2 = [Input.amplitude.*ones(1,high_event_idx_2 - low_event_idx_2 + 1),zeros(1,(high_idx_2-low_idx_2)- high_event_idx_2)];

% add 20 "zeros" at the beginning of the signals
N_zeros = 20;
y_1 = [y_baseline_norm((low_idx_1-N_zeros):(low_idx_1-1)),y_1];
y_2 = [y_baseline_norm((low_idx_2-N_zeros):(low_idx_2-1)),y_2];
u_1 = [zeros(1,N_zeros),u_1];
u_2 = [zeros(1,N_zeros),u_2];

%% 3)
time_1 = linspace((event_hour(1) - N_zeros*Params.t_sample)/3600,...
    (event_hour(1))/3600 + 1,length(y_1));
time_2 = linspace((event_hour(2) - N_zeros*Params.t_sample)/3600,...
    (event_hour(2))/3600 + 1,length(y_2));

figure;
subplot(2,2,1);
plot(time_1,u_1);
grid on;
xlim([time_1(1),time_1(end)]);
ylabel('Temperature set-point variation');

subplot(2,2,2);
plot(time_2,u_2);
grid on;
xlim([time_2(1),time_2(end)]);

subplot(2,2,3);
plot(time_1,y_1);
grid on;
xlim([time_1(1),time_1(end)]);
ylabel('Normalized aggregate power output')

subplot(2,2,4);
plot(time_2,y_2);
grid on;
xlim([time_2(1),time_2(end)]);

if(0)
    N_ON_1 = sum(squeeze(Results{1}.urec(1,:,:)).').';
    N_ON_2 = sum(squeeze(Results{2}.urec(1,:,:)).').';
    N_ON_3 = sum(squeeze(Results{3}.urec(1,:,:)).').';
    
    figure;
    subplot(1,3,1);
    plot(N_ON_1,'o');
    
    subplot(1,3,2);
    plot(N_ON_2,'o');
    
    subplot(1,3,3);
    plot(N_ON_3,'o');
end