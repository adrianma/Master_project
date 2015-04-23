%% ========================================================================
%  function main_Appendix_A
%  by Adrian Martinez Gomez
%  March 2015
%
%  Inputs:
%   *input_hours: hours to experiments with (good choice is [4,7,12,20])
%
%   *test_type: 
%       1: additivity 
%       2: homogeneity 
%       3: additivity with u_1 + u_2 = 0 (proposed by Maryam) 
%
%   *opt_sign: the sign of the external inputs (u_1,u_2)
%       1: both positive (PosPos)
%       2: both negative (NegNeg)
%       3: one positive and one negative (PosNeg)
%
%  Purpose:
%  Simulate examples to later on in script "errorsAppendix_A" assess if the
%  system is linear, or if this assumption is made, what level of error is
%  to be expected.
%  ========================================================================
function main_Appendix_A(input_hours,test_type,opt_sign);

%% 1) Core definitions
% number of simulations
N_sim = 100;
N_sim_big = 30;

% select 10 minutes as the duration of the input pulses
myDuration = 10;
% hours for the events
myHours = input_hours;

%% 2) Sampling the input space
%
% the first snippet does sampling for additivity
% the second snippet does sampling for homogeneity
%
if(test_type == 1)
    % create the inputs
    u_1_inputs = random('unif',0.2,0.45,[1,N_sim]);
    u_2_inputs = random('unif',0.2,0.45,[1,N_sim]);
    % add also some cases where 1 input is very big and other very small
    u_1_inputs = [u_1_inputs,random('unif',0.5,0.75,[1,N_sim_big])];
    u_2_inputs = [u_2_inputs,random('unif',0.1,0.2,[1,N_sim_big])];
    
    if(opt_sign == 1)
        disp('The inputs are positive (PosPos)');
    elseif(opt_sign == 2)
        disp('The inputs are negative (NegNeg)');
        % revert the sign
        u_1_inputs = -u_1_inputs;
        u_2_inputs = -u_2_inputs;
    elseif(opt_sign == 3)
        disp('The inputs are positive and negative (PosNeg)');
        % revert the sign of some of them
        u_1_inputs(1:N_sim) = -u_1_inputs(1:N_sim);
        
        % from the N_sim_big other experiments revert some of them with a
        % mask (to  ensure that not only the "big" inputs are reverted, but
        % rather a mix of "small" and "big")
        vMask = round(rand([1,N_sim_big]));
        u_1_inputs(N_sim + find(vMask == 1)) = -u_1_inputs;
        u_2_inputs(N_sim + find(vMask == 0)) = -u_2_inputs;
    end
    
    % add both inputs
    u_3_inputs = u_1_inputs + u_2_inputs;
    
    % display the sampling
    figure;
    subplot(2,1,1);
    scatter(u_1_inputs,u_2_inputs,'o');
    xlabel('u_1','FontSize',12);
    ylabel('u_2','FontSize',12);
    grid on;
    subplot(2,1,2);
    stem(u_3_inputs);
    ylabel('u_3 = u_1 + u_2','FontSize',12);
    grid on;
    
elseif(test_type == 2)
    u_1_inputs = random('unif',0,1,[1,5.*N_sim]);
    % these are the values for homogeneity
    alphas = random('unif',0.1,3,[1,5.*N_sim]);
    
    % select experiments where resulting input is <= 1 (requirement!)
    vIdx = find(alphas.*u_1_inputs <= 1);
    % shuffle them
    vIdxShuffled = randperm(numel(vIdx));
    vIdx_final = vIdx(vIdxShuffled(1:N_sim));
    % randomly select N_sim of the "allowed" data
    u_1_inputs = u_1_inputs(vIdx_final);
    alphas = alphas(vIdx_final);
    
    u_2_inputs = nan;
    
    u_3_inputs = alphas.*u_1_inputs;
    
    % display the sampling
    figure;
    subplot(2,1,1);
    scatter(u_1_inputs,alphas,'o');
    xlabel('u_1','FontSize',12);
    ylabel('\alpha','FontSize',12);
    grid on;
    subplot(2,1,2);
    stem(u_3_inputs);
    ylabel('u_3 = u_1 + u_2','FontSize',12);
    grid on;
    
elseif(test_type)
    % the idea is that both external inputs should add to zero
    % (u_1 + u_2 = 0)
    u_1_inputs = random('unif',0.2,0.6,[1,N_sim/2]);
    u_2_inputs = -u_1_inputs;
    u_3_inputs = nan;
    
    % display the sampling
    figure;
    scatter(u_1_inputs,u_2_inputs,'o');
    xlabel('u_1','FontSize',12);
    ylabel('u_2','FontSize',12);
    grid on;
    
    % redefine the opt_sign to have the folder be "PosNeg"
    opt_sign = 3;
else
    error('Select correct test_type value!');
end

%% 3) Load the necessary data for the simulations
%
if(0)
    idx = randi(20);
else
    idx = 3;
end
load(['archived_data/baseline/1000EWH/res_same/sim_1000EWHs_runNumber',...
    num2str(idx),'.mat']);
load('archived_data/baseline/1000EWH/same_baseline_data/power_baseline_same.mat',...
    'y_diffs_norm');
% load the water draws
load('archived_data/baseline/1000EWH/same_baseline_data/mean_water_1000EWH_same.mat',...
    'mean_water');

% define input struct
rInputStruct.idx = idx;
rInputStruct.mean_water = mean_water;
rInputStruct.y_diffs_norm = y_diffs_norm;
rInputStruct.Params = Params;
rInputStruct.PrelimModel = PrelimModel;
rInputStruct.WaterDrawScenarioReal = WaterDrawScenarioReal;
rInputStruct.Results_comparison = Results_comparison;

%% 4) Simulate the scenarios
%
yy = cell(length(myHours),length(u_3_inputs));
uu = cell(length(myHours),length(u_3_inputs));

if(test_type == 1)
    for ii = 1:(N_sim + N_sim_big)
        
        disp(['@@@@@@@@@@@@@@@@@@ ',num2str(ii),' @@@@@@@@@@@@@@@@@@']);
        myAmplitudes = [u_1_inputs(ii),u_2_inputs(ii),u_3_inputs(ii)];
        
        for jj = 1:length(myHours)
            [yy{jj,ii},uu{jj,ii}] = Appendix_A(rInputStruct,...
            myAmplitudes,myDuration,myHours(jj),0,0);
        end
    end
elseif(test_type == 2)
    for ii = 1:N_sim
        
        disp(['@@@@@@@@@@@@@@@@@@ ',num2str(ii),' @@@@@@@@@@@@@@@@@@']);
        myAmplitudes = [u_1_inputs(ii),u_3_inputs(ii)];
        
        for jj = 1:length(myHours)
            [yy{jj,ii},uu{jj,ii}] = Appendix_A(rInputStruct,...
                myAmplitudes,myDuration,myHours(jj),0,0);
        end
    end
elseif(test_type == 3)
    for ii = 1:N_sim
        
        disp(['@@@@@@@@@@@@@@@@@@ ',num2str(ii),' @@@@@@@@@@@@@@@@@@']);
        myAmplitudes = [u_1_inputs(ii),u_2_inputs(ii)];
        
        for jj = 1:length(myHours)
            [yy{jj,ii},uu{jj,ii}] = Appendix_A(rInputStruct,...
                myAmplitudes,myDuration,myHours(jj),0,0);
        end
    end
else
    error('Select correct test_type value!');
end

%% 5) Save the data accordingly
%
if(opt_sign == 1)
    t_opt_sign = 'PosPos';
elseif(opt_sign == 2)
    t_opt_sign = 'NegNeg';
elseif(opt_sign == 3)
    t_opt_sign = 'PosNeg';
end

if(test_type == 1)
    path_saving = ['archived_data/Linearity_SPV/',...
        t_opt_sign,'/my_sim_additivity_more_data.mat'];
elseif(test_type == 2)
    path_saving = ['archived_data/Linearity_SPV/',...
        t_opt_sign,'/my_sim_hom.mat'];
elseif(test_type == 3)
    path_saving = ['archived_data/Linearity_SPV/',...
        t_opt_sign,'/my_sim_zero.mat'];
end
% save results in file for future analysis
save(path_saving,'yy','uu','u_1_inputs','u_2_inputs','u_3_inputs',...
    'myDuration','myHours','opt_sign','test_type','idx','-v7.3');

fprintf(' ------> End of the script! <------\n');

end

%  Helping function Appendix_A_additivity
%
%  Inputs:
%   *myAmplitudes: 2x1 vector containing the amplitudes for the linearity
%   *myDuration: duration of the experiment
%   *myHour: containing the hour of the Event
%
%  Purpose:
%  Simulate Appendix cases
%  ========================================================================
function [yy,uu] = Appendix_A(rInputStruct,myAmplitudes,...
    myDuration,myHour,bSave,bPlot);

%% 1) Realize the SPV experiments
%
% store the steady-state power traj. for later analysis
SS_Prec = rInputStruct.Results_comparison.Prec;

global eta_in Events runs_events runs_eta;

Input.amplitude = myAmplitudes;
Input.duration = myDuration;

event_hour = repmat(myHour.*3600,1,length(Input.amplitude));
eta_in = Input.amplitude;

% define the sampling time and the duration of the simulation
rInputStruct.Params.t_sample = 10.0;
rInputStruct.Params.duration = (Input.duration*60/rInputStruct.Params.t_sample)*rInputStruct.Params.t_sample;
% water draws required
rInputStruct.Params.bWaterDraw = 1;

Results = cell(1,length(event_hour));
for ii = 1:length(Input.amplitude)
    disp(['------------ ',num2str(ii),' ------------']);
    
    runs_events = ii;
    runs_eta = ii;
    
    rInputStruct.Params.t_init = event_hour(ii);
    rInputStruct.Params.t_sim = event_hour(ii) + 1*3600;
    
    % get the state from this hour
    start_step = (rInputStruct.Params.t_init/rInputStruct.Params.t_sample + 1);
    rInputStruct.Params.x_init = rInputStruct.Results_comparison.xrec(:,start_step,:);
    rInputStruct.Params.u_init = rInputStruct.Results_comparison.urec(:,start_step,:);
    
    Events{ii} = rInputStruct.Params.t_init:1:(rInputStruct.Params.t_init + rInputStruct.Params.duration);
    tic;
    Results{ii} = simulate_population(rInputStruct.Params,rInputStruct.PrelimModel,...
        rInputStruct.WaterDrawScenarioReal,rInputStruct.WaterDrawScenarioReal,'SetPointVariation');
    toc;
end

%% 2) Gather the data and define the signals of interest
%
% define indeces
low_idx = event_hour/rInputStruct.Params.t_sample + 1;
high_idx = (event_hour + 1*3600)/rInputStruct.Params.t_sample + 1;
low_event_idx = Events{1}(1)/rInputStruct.Params.t_sample + 1 - low_idx;
high_event_idx = Events{1}(end)/rInputStruct.Params.t_sample + 1 - low_idx;

% number of simulation points
N_sim = length(Results{1}.Prec);
% number of zeros to "fill"
N_zeros = 20;

yy = zeros(N_sim + N_zeros,3);
uu_ext = zeros(N_sim,3);
u_exogenous = zeros(N_sim,3);
uu = zeros(N_sim + N_zeros,3*2);
for ii = 1:length(Input.amplitude)
    % output signal
    yy((N_zeros+1):end,ii) = Results{ii}.Prec - SS_Prec(low_idx(ii):high_idx(ii));
    % normalize
    yy((N_zeros+1):end,ii) = yy((N_zeros+1):end,ii)./sum(squeeze(rInputStruct.Params.P1_el));
    % get the external input signal
    uu_ext(:,ii) = [Results{ii}.etas(1).*ones(1,high_event_idx(ii) - low_event_idx(ii) + 1),zeros(1,(high_idx(ii)-low_idx(ii))- high_event_idx(ii))].';
    % get the exogenous input signal
    u_exogenous(:,ii) = rInputStruct.mean_water(low_idx(ii):high_idx(ii));
end

% "filling" for the output
yFill = rInputStruct.y_diffs_norm(rInputStruct.idx,...
    ((low_idx-N_zeros):(low_idx-1)));
% "filling" for the (exogenous) input
uFill = rInputStruct.mean_water((low_idx - N_zeros):(low_idx - 1)).';

for ii = 1:length(Input.amplitude)
    yy(:,ii) = [yFill.';yy(21:end,ii)];
end
for ii = 1:3
    % gather both inputs (with zeros)
    uu(:,ii:ii+1) = [[zeros(N_zeros,1);uu_ext(:,ii)],...
        [uFill;u_exogenous(:,ii)]];
end

%% 3) Save the calculated data
%
if(bSave)
    if(myAmplitudes(1)>0 && myAmplitudes(2)>0)
        tFolderSave = 'PosPos';
    elseif((myAmplitudes(1)>0 && myAmplitudes(2)<0) || (myAmplitudes(1)<0 && myAmplitudes(2)>0))
        tFolderSave = 'PosNeg';
    elseif(myAmplitudes(1)<0 && myAmplitudes(2)<0)
        tFolderSave = 'NegNeg';
    end
    
    path_saving = ['archived_data/Linearity_SPV/',tFolderSave,'/sim_',...
        num2str(myDuration),'duration_',strrep(num2str(myAmplitudes(1)),'.',''),'_',...
        strrep(num2str(myAmplitudes(2)),'.',''),'_',strrep(num2str(myAmplitudes(3)),'.',''),'_',...
        num2str(myHour),'hour'];
    
    % if the file is already stored, try storing in a new filename
    i_addition = 1;
    while(1)
        if(exist(strcat(path_saving,'.mat'),'file') == 2)
            path_saving = strcat(path_saving,num2str(i_addition));
            i_addition = i_addition + 1;
        else
            break;
        end
    end
    
    % save results in file for future analysis
    save(path_saving,'yy','uu','myAmplitudes','myDuration','myHour','-v7.3');
end

%% 4) Plot the experiment
%
% define the time vector
if(bPlot)
    xx_time = linspace((myHour*3600-N_zeros*Params.t_sample)/3600,...
        myHour + 1,length(yy{1}));
    
    figure;
    subplot(2,1,1);
    hold on;
    plot(xx_time,uu{1}(:,1) + uu{2}(:,1));
    plot(xx_time,uu{3}(:,1),'r');
    hold off;
    grid on;
    legend(['u_1 = ',num2str(uu{1}(N_zeros + 2,1)),' and u_2 = ',...
        num2str(uu{2}(N_zeros + 2,1))],['u_3 = ',...
        num2str(uu{3}(N_zeros + 2,1))],'Location','Best');
    % get the y-axis limits correctly
    if(uu{3}(N_zeros + 2,1) > 0)
        ylim([-0.1,uu{3}(N_zeros + 2,1)+0.1]);
    else
        ylim([uu{3}(N_zeros + 2,1)-0.1,0.1]);
    end
    xlim([(myHour*3600-N_zeros*Params.t_sample)/3600,myHour + 1]);
    ylabel('External input','FontSize',12);
    
    subplot(2,1,2);
    hold on;
    plot(xx_time,yy{1} + yy{2});
    plot(xx_time,yy{3},'r');
    hold off;
    grid on;
    legend('y_1 + y_2','y_3','Location','Best');
    xlim([(myHour*3600-N_zeros*Params.t_sample)/3600,myHour + 1]);
    xlabel('Time [h]','FontSize',12);
    ylabel('Normalized output signal','FontSize',12);
end

fprintf(' -> End of the run!\n');

end