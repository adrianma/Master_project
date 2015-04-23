%% ===================================================================
%  function linearity_property_2 (after Vaggelis's meeting)
%  by Adrian Martinez Gomez
%  October 2014
%
%  Purpose:
%   * do a case-study for linearity of aggregate power vs. parameter \eta
%   * generate some plots and display some observations/analysis from it
%   all.
%  ===================================================================
function Output = linearity_property(Method,NFolder);
% choose the folder for the method (PS:prob. switching,
% SPV:set-point variation)
if(strcmp(Method,'ProbSwitching'))
    tFolderMethod = 'PS';
elseif(strcmp(Method,'SetPointVariation'))
    tFolderMethod = 'SPV';
end
% get the subfolder for the simulation case
tFolder = strcat('/',num2str(NFolder));

folder_files = dir(['VAGGELIS/saved_data/',tFolderMethod,tFolder,...
    '/*.mat']);
load(folder_files(1).name);

%%
Output.eta = nan(1,length(folder_files)-1);
Output.P_var = ...
    nan(length(folder_files)-1,length(Results_comparison.Prec));
for ii = 1:(length(folder_files)-1)
    
    % load data-set
    load(folder_files(ii+1).name);
    disp(folder_files(ii+1).name);
    
    Output.P_var(ii,:) = Results.Prec - Results_comparison.Prec;
    
    % collect the \eta into the results struct
    Output.eta(ii) = Results.etas(1);   
end

% for convenience, sort the results
[Output.eta,vIdx] = sort(Output.eta);
Output.P_var = Output.P_var(vIdx,:);

% get the empirical mean and standard deviation of the variable part of the
% aggregate power
Output.mean_P_var = mean(Output.P_var.');
Output.std_P_var = std(Output.P_var.');

% time vector (just in case)
Output.vTime = Params.t_init:Params.t_sample:Params.t_sim;


%% Plot variable part of P_{agg}
fprintf('Generate some plots\n')

if(1)
    figure;
    hold on;
    plot(Output.eta,Output.mean_P_var,'*');
    if(strcmp(Method,'ProbSwitching') && 0)
        x = 0:0.1:1;
        plot(x,12.*x,'r');
    end
    hold off;
    grid on;
    xlabel('\eta');
    ylabel('Mean of the Variable part of P_{agg} [kW]');
    legend('Data points','Fit y=12x','Location','Best');
    
    figure;
    hold on;
    plot(Output.eta,Output.std_P_var,'*');
    % x = 0:0.1:1;
    % plot(x,12.*x,'r');
    hold off;
    grid on;
    xlabel('\eta');
    ylabel('Standard dev. of the Variable part of P_{agg} [kW]');
    legend('Data points','Fit y=12x','Location','Best');
end

% Mask for times after the time of the event t_event
Output.vIdx_time = find(Output.vTime > Results.Event(1));

% do a 3D plot of \eta, time, variable part of P_{agg}
if(0)
    figure;
    surf(Output.vTime(Output.vIdx_time)./3600,Output.eta,Output.P_var(:,Output.vIdx_time));
    xlabel('Time [h]');
    ylabel('\eta');
    zlabel('Variable part of the aggregate power [kW]');
end

% create a few 2d plots for different time steps after the disturbance
% (for linearity property)
if(0)
    s_time = min(Output.vIdx_time);
    
    while(s_time <= Output.vTime(end))
        figure;
        plot(Output.eta,Output.P_var(:,s_time),'*');
        grid on;
        xlabel('\eta');
        ylabel('Variable part of the aggregate power [kW]');
        
        % Plot every 20 min
        s_time = s_time + 1200; 
    end
    
end

% plot the ss for some cases (to see the response of the variable part)
if(1)
    % \eta_1 = 0.1 + \eta_2 = 0.2 equal \eta_3 = 0.3 ? 
    figure;
    hold on; 
    plot(Output.vTime(Output.vIdx_time)./3600,Output.P_var(1,Output.vIdx_time) + ...
        Output.P_var(2,Output.vIdx_time));
    plot(Output.vTime(Output.vIdx_time)./3600,Output.P_var(3,Output.vIdx_time),'r');
    grid on;
    hold off;
    xlabel('Time [h]');
    ylabel('Variable part of the aggregate power P_{agg} [kW]');
    legend(['\eta_1 =',num2str(Output.eta(1)),' + \eta_2 =',...
        num2str(Output.eta(2))],['\eta_3 =',num2str(Output.eta(3))],...
        'Location','Best');
    
    % add info on how good the linearity is
    Output.RMSE_res(1) = RMSE(Output.P_var(1,Output.vIdx_time)+Output.P_var(2,Output.vIdx_time),Output.P_var(3,Output.vIdx_time));
    Output.rel_error(1) = norm(Output.P_var(1,Output.vIdx_time)+Output.P_var(2,Output.vIdx_time),2)/norm(Output.P_var(3,Output.vIdx_time),2);
    title(['RMSE in [kW]: ',num2str(Output.RMSE_res(1)),'. Rel. error:',num2str(Output.rel_error(1)),'%']);
    
    % \eta_1 = 0.4 + \eta_2 = 0.4 equal \eta_3 = 0.8 ? 
    figure;
    hold on; 
    plot(Output.vTime(Output.vIdx_time)./3600,Output.P_var(4,Output.vIdx_time) + ...
        Output.P_var(4,Output.vIdx_time));
    plot(Output.vTime(Output.vIdx_time)./3600,Output.P_var(8,Output.vIdx_time),'r');
    grid on;
    hold off;
    xlabel('Time [h]');
    ylabel('Variable part of the aggregate power P_{agg} [kW]');
    legend(['\eta_1 =',num2str(Output.eta(4)),' + \eta_2 =',...
        num2str(Output.eta(4))],['\eta_3 =',num2str(Output.eta(8))],...
        'Location','Best');
    
    % add info on how good the linearity is
    Output.RMSE_res(2) = RMSE(Output.P_var(4,Output.vIdx_time)+Output.P_var(4,Output.vIdx_time),Output.P_var(8,Output.vIdx_time));
    Output.rel_error(2) = norm(Output.P_var(4,Output.vIdx_time)+Output.P_var(4,Output.vIdx_time),2)/norm(Output.P_var(8,Output.vIdx_time),2);
    title(['RMSE in [kW]: ',num2str(Output.RMSE_res(2)),'. Rel. error:',num2str(Output.rel_error(2)),'%']);
end
if(0)
    for ii = 1:length(folder_files)-1
        figure;
        plot(Output.vTime(Output.vIdx_time)./3600,Output.P_var(ii,Output.vIdx_time));
        grid on;
        xlabel('Time [h]');
        ylabel('Variable part of the aggregate power P_{agg} [kW]');
        legend(['\eta_1 =',num2str(Output.eta(ii))],'Location','Best');
    end
end

end

% Function that calculates the Root Mean Squared Error
%
%
function RMSE_res = RMSE(vec_est,vec_real)
RMSE_res = 0;
for ii = 1:length(vec_est)
    RMSE_res = RMSE_res + (vec_est(ii) - vec_real(ii))^2;
end
RMSE_res = sqrt(RMSE_res/length(vec_est));

end