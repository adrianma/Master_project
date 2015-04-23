%% -------------------------------------------------------------------
%  function tests_Analyze()
%  by Adrian Martinez Gomez
%  October 2014
%
%  Purpose:
%   * perform reads of analysis in data
%  -------------------------------------------------------------------
function DataAnalysis = tests_Analyze(tFolderMethod,NAnalysis);
%% do some data analysis
%
global eta_in Events;

% generate a new folder for the simulation case
if(NAnalysis~=0)
    tFolder = strcat('/',num2str(NAnalysis));
else
    tFolder = '';
end

% get the file names
prob_PS_files = dir(['VAGGELIS/saved_data/',tFolderMethod,tFolder,...
    '/*.mat']);
load(prob_PS_files(1).name);

% time vector
DataAnalysis.vTime = Params.t_init:Params.t_sample:Params.t_sim;

% number of simulation data sets
DataAnalysis.sNDataSets = length(prob_PS_files);

% for settling time, converge towards 5% error from internal controllers 
% only 
DataAnalysis.epsilon = 0.05;

%
% preallocate
% (names should be self-explanatory)
DataAnalysis.urec = nan(length(DataAnalysis.vTime),DataAnalysis.sNDataSets);
DataAnalysis.urec_comp = nan(length(DataAnalysis.vTime),DataAnalysis.sNDataSets);
DataAnalysis.mean_Prec = nan(1,DataAnalysis.sNDataSets);
DataAnalysis.mean_Prec_comp = nan(1,DataAnalysis.sNDataSets);
DataAnalysis.var_Prec = nan(1,DataAnalysis.sNDataSets);
DataAnalysis.var_Prec_comp = nan(1,DataAnalysis.sNDataSets);
DataAnalysis.res_RMSE = nan(1,DataAnalysis.sNDataSets);
DataAnalysis.val_events = nan(length(Events),DataAnalysis.sNDataSets);
DataAnalysis.settling_time = nan(length(Events),DataAnalysis.sNDataSets);
DataAnalysis.N_ON = nan(length(DataAnalysis.vTime),DataAnalysis.sNDataSets);
DataAnalysis.N_ON_comp = nan(length(DataAnalysis.vTime),DataAnalysis.sNDataSets);
DataAnalysis.etas = nan(1,DataAnalysis.sNDataSets);

% load the internal-controller only file first
load(prob_PS_files(1).name);

% for all files in directory perform the analysis
%
% start loop at 3 because the first 2 files are the internal controller
% ones
for ii = 3:DataAnalysis.sNDataSets
    
    % load data-set
    load(prob_PS_files(ii).name);
    disp(prob_PS_files(ii).name);
    DataAnalysis.etas(ii) = Results.etas(1);
    
    %%
    % save the input vector(s) take. Reduce them to the sensor working
    % (always sensor 2, per Vrettos definition)
    sIdx = 1 + (sum(sum(Results.urec(2,:,:))) ~= 0);
    
    DataAnalysis.urec(:,ii) = sum(squeeze(Results.urec(sIdx,:,:)).');
    DataAnalysis.urec_comp(:,ii) = ...
        sum(squeeze(Results_comparison.urec(sIdx,:,:)).');
    
    % also, store the set of EWHs that are ON
    for jj = 1:length(DataAnalysis.vTime)
        DataAnalysis.N_ON(jj,ii) = sum(Results.urec(sIdx,jj,:));
        DataAnalysis.N_ON_comp(jj,ii) = ...
            sum(Results_comparison.urec(sIdx,jj,:));
    end
    
    %%
    % compute means
    DataAnalysis.mean_Prec(ii) = mean(Results.Prec);
    DataAnalysis.mean_Prec_comp(ii) = mean(Results_comparison.Prec);
    
    % compute variances
    DataAnalysis.var_Prec(ii) = var(Results.Prec);
    DataAnalysis.var_Prec_comp(ii) = var(Results_comparison.Prec);
    
    % calculate the RMSE
    DataAnalysis.res_RMSE(ii) = RMSE(Results.Prec,Results_comparison.Prec);
    
    % store the value(s) caused by the event(s)
    for jj = 1:length(Events)
        DataAnalysis.val_events(jj,ii) = ...
            Results.Prec(DataAnalysis.vTime == ...
            Events{jj}(1).*ones(1,length(DataAnalysis.vTime)));
    end
    
    %%
    % -----------------------------%
    % find how much time does the disturbance need to settle within epsilon
    % from the "internal-controller" only variant
    %   -> settling-time
    
    % get normalized power vectors
    norm_Prec = Results.Prec./norm(Results.Prec);
    norm_Prec_comparison = ...
        Results_comparison.Prec./norm(Results_comparison.Prec);
    
    v_diff = abs(norm_Prec - norm_Prec_comparison)./...
        abs(norm_Prec);
    
    for jj = 1:length(Events)
        
        % find the coefficients for which the differences are maintained
        coeffs_temp = find(v_diff <= DataAnalysis.epsilon);
        coeffs = coeffs_temp(coeffs_temp.*Params.t_sample > Events{jj}(1));
        
        if(isempty(coeffs))
            DataAnalysis.settling_time(jj,ii) = nan;
        else
            % settling time in seconds
            DataAnalysis.settling_time(jj,ii) = ...
                DataAnalysis.vTime(min(coeffs)) - Events{jj}(1);
        end
        
    end
    %-----------------------------%
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