%% ========================================================================
%  script Chapter_3_3
%  by Adrian Martinez Gomez
%  March 2015
%
%  Purpose:
%  Calculate some examples (and generate plot) for Probability switching
%  This corresponds to Section 3.4.1 of the written Master Thesis report.
%
%  ========================================================================
%% 1) Loading and initializations
%
% load the file with baseline data
if(0)
    idx = randi(20);
else
    idx_file = 3;
end
load(['archived_data/baseline/1000EWH/res_same/sim_1000EWHs_runNumber',...
    num2str(idx),'.mat']);
disp([' =====> The idx of the loaded baseline is ',num2str(idx),' <=====']);
fprintf('\n');

% define gobal variables needed later on
global eta_in runs runs_events runs_eta Events event_hour;
% external control signal is % of probability switching
Method = 'ProbSwitching';
% take hour 18 as start of the event
event_hour = 18*3600;
% population of 1000 EWHs
n_app = 1000;

% take an example of 
eta_in = [0.3,-0.5];
input_duration = 10;

% define the simulated times (initial, simulation, sampling, duration of 
% pulse)
Params.t_sample = 10.0;
Params.t_init = event_hour;
if(1)
    Params.t_sim = event_hour + 1*3600;
else
    Params.t_sim = 24*3600;
end
Params.duration = input_duration*60;

start_step = (Params.t_init/Params.t_sample + 1);
Params.x_init = Results_comparison.xrec(:,start_step,:);
Params.u_init = Results_comparison.urec(:,start_step,:);

Events{1} = Params.t_init:1:(Params.t_init+Params.duration);
Events{2} = Events{1};

% water draws required
Params.bWaterDraw = 1;

%% 2) Simulate
%
% struct where the simulations with probability switching are stored
Results = cell(1,2);
for runs = 1:2
    runs_events = runs;
    runs_eta = runs;
    
    disp(['%%%%%%%%%%%% RUN NUMBER ',num2str(runs),' %%%%%%%%%%%%']);
    tic;
    Results{runs} = simulate_population(Params,PrelimModel,...
        WaterDrawScenarioReal,WaterDrawScenarioReal,Method);
    toc;
end

%% 3) Preparing the signals
%
% define low and high indeces for the event
low_idx = event_hour/Params.t_sample + 1;
high_idx = (event_hour + 1*3600)/Params.t_sample + 1;
low_event_idx = Events{1}(1)/Params.t_sample + 1 - low_idx;
high_event_idx = Events{1}(end)/Params.t_sample + 1 - low_idx;
% define N_zeros time steps before the pulse is applied
N_zeros = 10;

% prepare the aggregate power
yy = cell(1,2);
yy{1} = [Results_comparison.Prec((low_idx-N_zeros):(low_idx-1)),Results{1}.Prec];
yy{2} = [Results_comparison.Prec((low_idx-N_zeros):(low_idx-1)),Results{2}.Prec];
% get the reference signal
yy_ref = Results_comparison.Prec((low_idx-N_zeros):(high_idx));
% normalize
yy{1} = yy{1}./sum(Params.P1_el);
yy{2} = yy{2}./sum(Params.P1_el);
yy_ref = yy_ref./sum(Params.P1_el);

% prepare the external input
u{1} = [eta_in(1).*ones(1,high_event_idx - low_event_idx + 1),zeros(1,(high_idx-low_idx)- high_event_idx)];
u{1} = [zeros(1,N_zeros),u{1}];
u{2} = [eta_in(2).*ones(1,high_event_idx - low_event_idx + 1),zeros(1,(high_idx-low_idx)- high_event_idx)];
u{2} = [zeros(1,N_zeros),u{2}];

%% 4) Plotting the signals
%
% get the number of loads
% xx = linspace(event_hour/3600 - N_zeros*Params.t_sample/3600,...
%     event_hour/3600 + 1,length(yy{1}));
%
% define the time vector
xx = linspace(event_hour/3600,event_hour/3600 + 1,length(yy{1}));

figure;
subplot(2,2,1);
plot(xx,u{1});
grid on;
ylim([0,0.6]);
ylabel('Probabilistic switching input','FontSize',12);

subplot(2,2,2);
plot(xx,u{2});
grid on;
ylim([-0.6,0]);

subplot(2,2,3);
hold on;
plot(xx,yy{1},'r');
plot(xx,yy_ref);
hold off;
grid on;
ylabel('Normalized aggregate power','FontSize',12);

subplot(2,2,4);
hold on;
plot(xx,yy{2},'r');
plot(xx,yy_ref);
hold off;
grid on;
legend('Probabilistic switching','Reference baseline','Location','Best');