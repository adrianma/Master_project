%% ========================================================================
%  function main_CL_simulation
%  by Adrian Martinez Gomez
%  December 2014
%
%  Modified by Adrian Martinez Gomez
%  March/April 2015
%
%  Output: none
%  Input:
%   *ss: the state-space system (linear model)
%   *Data: iddata object from which the reference signal is taken
%   with the input trajectory calculated.
%   *rRefining: struct with different simulation entries(basically choices)
%   (N_horizon,p_terminal,bHalf,idx,CostFunctionCoeff,bSave,bArtificial
%   bRealSystem).
%
%  Purpose:
%   * calculate the solutions for the GUROBI solver for 1 time step with
%   MPC and (C)QP.
%   * matrices are defined here (this is the most difficult part by far)
%  ========================================================================
function main_CL_simulation(ss,Data,rRefining);
%% 1) Initializations and definitions

% 1.1) get the idx for the reference signal
ref_exo_idx = rRefining.idx;

% 1.2) load the exogenous input signal (averaged water draw)
load('archived_data/baseline/1000EWH/same_baseline_data/mean_water_1000EWH_same.mat','mean_water');

% 1.3) arrange the signals (exogenous input and reference) + display them
%
% convert the idx into an hour, and then the hour into step of the day
event_hour = idxToHour(ref_exo_idx);
event_hour = event_hour*3600;
disp([' => Hour of the event: ',num2str(event_hour./3600)]);
disp('');

% get the indeces for the exogenous input, where the event starts
low_idx = event_hour/10 + 1;
high_idx = (event_hour + 1*3600)/10 + 1;
% number of "filling zeros" is where there is no external input u_k
N_zeros = 20;
% shift the indeces to account for the extra zeros at the beginning
low_idx = low_idx - N_zeros;
high_idx = high_idx - N_zeros;
    
% try artificial signals as suggested by Maryam (ones that don't actuallty
% come from the system itself)
if(rRefining.bArtificial)
    % ============================================================== %
    ref_signal = Different_Reference(rRefining.ArtificialType);
    % ============================================================== %
else
    % load the reference signal from the validation dataset
    load(['archived_data/baseline/1000EWH/',...
        'same_baseline_data/power_baseline_same.mat'],'y_diffs_norm');
    ref_signal =  [Data.y{ref_exo_idx}.',...
        y_diffs_norm(mod(ref_exo_idx-1,20)+1,((high_idx + 1):(high_idx + ...
        4*N_zeros - 1)))];
end

% define the exogenous input
exogenous_input = mean_water(low_idx:(high_idx + 4*N_zeros - 1));
 
% 1.4) display the signals
%
% compare options
opt_comp = compareOptions('InitialCondition','z');
[y_LM,~,~] = compare(getexp(Data,rRefining.idx),ss,opt_comp);
time_vec_plotting = linspace((event_hour-N_zeros*10)/3600,...
    event_hour/3600 + 1,length(y_LM.y));

figure;
if(~rRefining.bArtificial)
    subplot(2,1,1);
    plot(time_vec_plotting,Data.u{rRefining.idx}(:,1));
    title(['Hourly interval starts at = ',num2str(event_hour./3600)])
    legend('Exogenous input w_{k}','Location','Best');
    grid on;
    % get the y-axis limits correctly
    if(Data.u{rRefining.idx}(23,1) > 0)
        ylim([-0.1,Data.u{rRefining.idx}(23,1)+0.1]);
    else
        ylim([Data.u{rRefining.idx}(23,1)-0.1,0.1]);
    end
    xlim([(event_hour-N_zeros*10)/3600,event_hour/3600 + 1]);
    
    subplot(2,1,2);
    hold on;
    plot(time_vec_plotting,ref_signal(1:length(y_LM.y)));
    plot(time_vec_plotting,y_LM.y,'r');
    hold off;
    grid on;
    legend('Reference signal ref_{k}','ARMAX model response','Location',...
        'Best');
    xlim([(event_hour-N_zeros*10)/3600,event_hour/3600 + 1]);
    xlabel('Time [h]');
else
   plot(time_vec_plotting,ref_signal(1:length(y_LM.y))); 
   grid on;
   xlim([(event_hour-N_zeros*10)/3600,event_hour/3600 + 1]);
   xlabel('Time [h]');
   legend('Reference signal ref_{k}','Location','Best');
end
close;
pause(3);

% 1.5) check the system properties (stability, observability,
% controllability)
Check_model_properties(ss);

% 1.6) other definitions
% define the time horizon (N_horizon * Ts [seconds])
N_horizon = rRefining.N_horizon;
% usually set to 0 to handle the exogenous-input part correctly
if(size(ss.b,2) == 2)
    bExogenous = 0;
else
    bExogenous = 1;
end
% number of simulation steps
%N_simulation = length(ref_signal) - N_horizon;
N_simulation = 380;
disp([' ----> The horizon is going to be ',num2str(N_horizon*10),...
    ' [seconds] <----']);
pause(2);

%% 2) Simulate
%
% pre-allocate/define needed vals
Xpred = cell(1,N_simulation);
Upred = cell(1,N_simulation);
Ypred = cell(1,N_simulation);
run_time = nan(N_simulation,1);

% store the first values
Xpred{1} = zeros(size(ss.a,1),1);

for ii = 1:N_simulation-1
    
    % call the solver to get the optimization for 1 time-step
    Results = solve_CQP_problem(ss,Xpred{ii},N_horizon,bExogenous,...
        rRefining,ref_signal,exogenous_input,ii);
    
    % store the predicted control value
    Upred{ii} = Results.upred{1};
    % store the output
    if(~bExogenous)
        Ypred{ii} = ss.c*Xpred{ii} + ss.d*[Upred{ii};...
            mean_water(ii)];
    else
        Ypred{ii} = ss.c*Xpred{ii} + ss.d*Upred{ii};
    end
    
    % calculate new state
    %   (remember that we disregarded the noise!)
    if(~bExogenous)
        Xpred{ii+1} = ss.a*Xpred{ii} + ss.b(:,1)*Upred{ii} + ...
            ss.b(:,2)*exogenous_input(ii);
    else
        Xpred{ii+1} = ss.a*Xpred{ii} + ss.b(:,1)*Upred{ii};
    end
    
    if(norm(Xpred{ii} - Results.xpred{1}) ~= 0)
        error('Careful, values not the same!!!')
    end
    
    run_time(ii) = Results.gurobi_res.runtime;
    
    % display run-time
    if(mod(ii,360) == 0)
        disp(['Optimization time of the day: ',num2str(ii/360),' h']);
    end
end
% get the last output
if(~bExogenous)
    Ypred{ii+1} = ss.c*Xpred{ii+1} + ss.d*[Upred{ii};mean_water(ii)];
else
    Ypred{ii+1} = ss.c*Xpred{ii+1} + ss.d*Upred{ii};
end

%% 3) Save the results in a struct
disp('');
disp(' => Finished the whole thing!');
% store the runtime, and predicted values into matrices
tot_run_time = sum(run_time(1:N_simulation-1));
disp('');
disp([' => Total gurobi runtime: ',num2str(tot_run_time)]);

matUpred = cell2mat(Upred);
matYpred = cell2mat(Ypred);
matXpred = cell2mat(Xpred);

% get the results into an output Struct
ResultStruct.tot_run_time = tot_run_time;
ResultStruct.matUpred = matUpred;
ResultStruct.matYpred = matYpred;
ResultStruct.matXpred = matXpred;
ResultStruct.ref_signal = ref_signal;
ResultStruct.N_horizon = N_horizon;
ResultStruct.N_simulation = N_simulation;
ResultStruct.ref_exo_idx = ref_exo_idx;
ResultStruct.rRefining = rRefining;

%% 4) Apply the optimal input trajectory to the real system (and simulate)
%
% call the external function (which then simulates the real system with the
% given external control trajectory)
if(rRefining.bRealSystem)
    [ResultsOptimal,Params,SS_Prec] = simulate_on_real_system(N_zeros,...
        ref_exo_idx,ResultStruct.matUpred);
    
    idx_start = event_hour/10;
    idx_end = idx_start + length(ref_signal) - 1;
    % normalized signal
    y_RS = (ResultsOptimal.Prec - SS_Prec(idx_start:idx_end - N_zeros));
    y_RS = y_RS./sum(squeeze(Params.P1_el));
    
    % plot a result
    if(0)
        figure;
        plot(y_RS);
        xlabel('Time index k','FontSize',12);
        ylabel('y^{*}','FontSize',12);
        grid on;
    end
    
    fprintf('\n Finished simulating for the real system! \n');
end

%% 5) Save output structs (for linear model and for real system's applied
if(rRefining.bSave &&  rRefining.bRealSystem)
    save(['archived_data/SPV/ClosedLoop/',t_Day_Saving,'Matrices/',...
        t_Model_Saving,'/Predicted_matrices2_p',num2str(...
        rRefining.p_terminal),'_N',num2str(rRefining.N_horizon),'_bHalf',...
        num2str(rRefining.bHalf),'_idx',num2str(ref_exo_idx),'_CFCoeff',...
        strrep(num2str(sum(rRefining.CostFunctionCoeff)),'.','')],'ResultStruct','ResultsOptimal');
elseif(rRefining.bSave)
    save(['archived_data/SPV/ClosedLoop/',t_Day_Saving,...
        'Matrices/',t_Model_Saving,'/Predicted_matrices1_p',num2str(...
        rRefining.p_terminal),'_N',num2str(rRefining.N_horizon),'_bHalf',...
        num2str(rRefining.bHalf),'_idx',num2str(ref_exo_idx),'_CFCoeff',...
        strrep(num2str(sum(rRefining.CostFunctionCoeff)),'.','')],'ResultStruct');
end

%% 6) Plots
% 6.1) plot figure with: U^* and ref=real_system(U^*) + l.model(U_ref)
figure;
subplot(2,1,1);
plot(ResultStruct.matUpred);
legend('U_{k}^{*}','Location','Best');
grid on;
xlabel('Index k');
xlim([0,length(ref_signal)]);

title(['N_{horizon} = ',num2str(rRefining.N_horizon),' P_{terminal} = ',...
    num2str(rRefining.p_terminal),' bHalf = ',num2str(rRefining.bHalf),...
    '_idx',num2str(ref_exo_idx),' Cost Function coeffs = ',...
    strrep(num2str(rRefining.CostFunctionCoeff(1)),'.',''),',',...
    strrep(num2str(rRefining.CostFunctionCoeff(2)),'.',''),',',...
    strrep(num2str(rRefining.CostFunctionCoeff(3)),'.',''),',',...
    strrep(num2str(rRefining.CostFunctionCoeff(4)),'.','')]);

subplot(2,1,2);
hold on;
plot(ResultStruct.matYpred,'r');
plot(ResultStruct.ref_signal);
hold off;
legend('y_k','Ref_k','Location','Best');
grid on;
xlabel('Index k');
xlim([0,length(ref_signal)]);

if(rRefining.bSave)
    % save figure
    saveas(gcf,['archived_data/SPV/ClosedLoop/',...
        t_Day_Saving,'Figures/',t_Model_Saving,'/my_fig1_p',...
        num2str(rRefining.p_terminal),'_N',num2str(rRefining.N_horizon),...
        '_bHalf',num2str(rRefining.bHalf),'_idx',num2str(ref_exo_idx),...
        '_InpCost',strrep(num2str(sum(rRefining.CostFunctionCoeff)),'.','')],'fig');
end

if(rRefining.bRealSystem)
    % 6.2) plot figure with: U^* and y^*=real_system(U^*)
    figure;
    subplot(2,1,1);
    plot(ResultStruct.matUpred);
    legend('U_{k}^{*}','Location','Best');
    grid on;
    xlabel('Index k');
    xlim([0,length(ref_signal)]);
    
    title(['N_{horizon} = ',num2str(rRefining.N_horizon),' P_{terminal} = ',...
        num2str(rRefining.p_terminal),' bHalf = ',num2str(rRefining.bHalf),...
        '_idx',num2str(ref_exo_idx),' Cost Function coeffs = ',...
        strrep(num2str(rRefining.CostFunctionCoeff(1)),'.',''),',',...
        strrep(num2str(rRefining.CostFunctionCoeff(2)),'.',''),',',...
        strrep(num2str(rRefining.CostFunctionCoeff(3)),'.',''),',',...
        strrep(num2str(rRefining.CostFunctionCoeff(4)),'.','')]);
    
    subplot(2,1,2);
    hold on;
    plot(ResultStruct.ref_signal,'k');
    plot(ResultStruct.matYpred,'r');
    plot(y_RS(1:length(ResultStruct.matYpred)));
    hold off;
    legend('Ref_{k}','y_{k}^{*}','y_{k}','Location','Best');
    grid on;
    xlabel('Index k');
    xlim([0,length(y_RS)]);
    
    if(rRefining.bSave)
        % save figure
        saveas(gcf,['archived_data/SPV/ClosedLoop/',t_Day_Saving,'Figures/',...
            t_Model_Saving,'/my_fig2_p',num2str(rRefining.p_terminal),'_N',...
            num2str(rRefining.N_horizon),'_bHalf',num2str(rRefining.bHalf),...
            '_idx',num2str(ref_exo_idx),'_InpCost',...
            strrep(num2str(sum(rRefining.CostFunctionCoeff)),'.','')],'fig');
    end
    
    % 6.3) plot figure with: U^* and U_ref
    if(0)
    figure;
    subplot(2,1,1);
    plot(ResultStruct.matUpred);
    legend('U_{k}^{*}','Location','Best');
    grid on;
    xlabel('Index k');
    xlim([0,length(ref_signal)]);
    
    title(['N_{horizon} = ',num2str(rRefining.N_horizon),' P_{terminal} = ',...
        num2str(rRefining.p_terminal),' bHalf = ',num2str(rRefining.bHalf),...
        '_idx',num2str(ref_exo_idx),' Cost Function coeffs = ',...
        strrep(num2str(rRefining.CostFunctionCoeff(1)),'.',''),',',...
        strrep(num2str(rRefining.CostFunctionCoeff(2)),'.',''),',',...
        strrep(num2str(rRefining.CostFunctionCoeff(3)),'.',''),',',...
        strrep(num2str(rRefining.CostFunctionCoeff(4)),'.','')]);
    
    % this is the input from which the reference signal was generated in the
    % real system
    reference_input = [zeros(1,N_zeros),Data.u{...
        ref_exo_idx}(:,1).',zeros(1,3*N_zeros - 1)];
    
    subplot(2,1,2);
    plot(reference_input);
    hold off;
    legend('U_{ref}','Location','Best');
    grid on;
    xlabel('Index k');
    xlim([0,length(ref_signal)]);
    end
    
    if(rRefining.bSave)
        % save figure
        saveas(gcf,['archived_data/SPV/ClosedLoop/',t_Day_Saving,'Figures/',...
            t_Model_Saving,'/my_fig3_p',num2str(rRefining.p_terminal),'_N',...
            num2str(rRefining.N_horizon),'_bHalf',num2str(rRefining.bHalf),...
            '_idx',num2str(ref_exo_idx),'_InpCost',...
            strrep(num2str(sum(rRefining.CostFunctionCoeff)),'.','')],'fig');
    end
end

%% BIN
if(0)
    %
    % possible struct for running the simulations
    
    rRefining.N_horizon = 30;
    rRefining.p_terminal = 20;
    rRefining.bHalf = 0;
    rRefining.idx = 170;
    rRefining.u_max = 0.95;
    rRefining.u_min = 0.0;
    rRefining.bArtificial = 0;
    rRefining.bRealSystem = 0;
    rRefining.CostFunctionCoeff = [1,0,0,0];
    rRefining.bSave = 0;
    
    
end
end