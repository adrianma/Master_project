%% ===================================================================
%  function Models_system_ID
%  by Adrian Martinez Gomez
%  November 2014
%
%  Last modified:
%  December 2014
%
%  Inputs:
%   * Method: string with the method selected (prob. switching or set-point
%   variation.
%   * NFolders: the number of folder from which we want to read files
%	  (starting with 1,2,...,NFolders)
%	* n_app: number of loads for the population
%   * Exogenous_option: number giving options for the different
%     exogenous input options (1:water draw sum, 2:averaged water draw sum,
%     3:averaged steady-state aggregate power, 4:no exogenous input).
%
%  Purpose:
%   * compute the iddata object experiments out from the Amplitudes
%     experiments out from 'main_experiments'.
%  ===================================================================
function Models =  Models_system_ID(Method,n_app,Exogenous_option);

%% 1) Gather iddata object
%
% choose the folder for the method (PS:prob. switching,
% SPV:set-point variation)
if(strcmp(Method,'ProbSwitching'))
    tFolderMethod = 'PS';
elseif(strcmp(Method,'SetPointVariation'))
    tFolderMethod = 'SPV';
end

% load the file with the iddata objects from the files
switch Exogenous_option
    case 0
        error('This file (without ex. input) was not generated yet!');
    case 1
        % load corresponding file
        load(['archived_data/',tFolderMethod,'/',num2str(n_app),...
            'EWH/Data/SUM_DATA_',num2str(n_app),'_EWH.mat']);
        
        % gather the data to build the model later on
        if(n_app==500)
            data = data_sum_500;
        elseif(n_app==1000)
            data = data_sum_1000;
        end
    case 2
        % load corresponding file
        load(['archived_data/',tFolderMethod,'/',num2str(n_app),...
            'EWH/Data/SUM_AVERAGED_DATA_',num2str(n_app),'_EWH.mat']);
        
        % gather the data to build the model later on
        if(n_app==500)
            data = data_sum_averaged_500;
        elseif(n_app==1000)
            data = data_sum_averaged_1000;
        end
    case 3
        % load corresponding file
        load(['archived_data/',tFolderMethod,'/',num2str(n_app),...
            'EWH/Data/POW_AVERAGED_DATA_',num2str(n_app),'_EWH.mat']);
        
        % gather the data to build the model later on
        if(n_app==500)
            data = data_pow_averaged_500;
        elseif(n_app==1000)
            data = data_pow_averaged_1000;
        end
end

%% 2) Construct model
% detrend the data
data_d = detrend(data);
% select the order of the system
order = 3;
if(0)
    NN = struc(1:5,[1:5,1:5],[0 0]);
    V = arxstruc(data(1:length(data.Period)/2),...
        data(length(data.Period)/2+1,end),NN);
    nn = selstruc(V,0);
    m = arx(z,nn);
end
% create the state space system with N4SID
%
% setting DisturbanceModel to none sets the matrix K to zero
opt = n4sidOptions('Focus','simulation','N4Weight','auto','Display','on');
sys = n4sid(data_d,order,opt,'DisturbanceModel','none');

% compare model to real response
[YY,fit,XX0] = compare(validation_data,sys);
% find the indeces for which the fit is very good (greater than 90%)
vIdx_good = find(cell2mat(fit) > 90);
N_exp_good = validation_data.ExperimentName(vIdx_good);
disp(['There are ',num2str(length(vIdx_good)),...
    ' experiments with fit above 90%']);
% plot them all!
for ii = 1:length(vIdx_good)
    figure;
    subplot(2,1,1);
    hold on;
    plot(validation_sum_averaged_1000.u{vIdx_good(ii)}(:,1));
    plot(validation_sum_averaged_1000.u{vIdx_good(ii)}(:,2),'r');
    hold off;
    grid on;
    xlabel('Index k');
    legend('Set-point variation u_k','Exogenous input w_k','Location','Best');
    subplot(2,1,2);
    hold on;
    plot(validation_sum_averaged_1000.y{vIdx_good(ii)}.');
    plot(YY_stab_noisei{vIdx_good(ii)}.y,'r');
    hold off;
    grid on;
    legend('Real system ref_k','Linear model y_k','Location','Best');
    xlabel('Index k');
end


Models = sys;
        
end