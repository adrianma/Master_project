%% ===================================================================
%  function estimate_N_ON_populations
%  by Adrian Martinez Gomez
%  November 2014
%
%  Inputs:
%
%  Purpose:
%   * Estimate sets N_On and N_off from the aggregate power measured
%  ===================================================================
function OutputStruct = estimate_N_ON_populations(bParts)

%% set to 1 to gather the sets and aggregate powers
if(bParts(1))
    folder_files = dir('VAGGELIS/saved_data/baseline/*.mat');
    % N_E: number of experiments
    OutputStruct.N_E = length(folder_files);
    
    % loop through all the simulation files to get the aggregate power and 
    % N_on set vectors for all times t
    %
    OutputStruct.aggregate_power = cell(OutputStruct.N_E,1);
    OutputStruct.sets_N_on = cell(OutputStruct.N_E,1);
    OutputStruct.sets_N_off = cell(OutputStruct.N_E,1);
    count = 1;
    for ii = 1:OutputStruct.N_E
        % store the number of loads, the time vector and the period of the
        % system
        if(count == 1)
            OutputStruct.n_app = Params.n_app;
            OutputStruct.vTime = ...
                Params.t_init:Params.t_sample:Params.t_sim;
            OutputStruct.period = Params.t_sample;
        end
        
        % load the ii-th simulation file
        load(folder_files(ii).name);
        
        % get the aggregate power
        OutputStruct.aggregate_power{ii} = Results.Prec;
        % get the total rated power (for normalization purposes later on)
        OutputStruct.rated_power{ii} = sum(Params.P1_el);
        
        % get the measurement of the set N_on
        N_on = nan(length(Results.Prec),1);
        for jj = 1:length(Results.Prec);
            N_ON_1 = sum(Results.urec(1,jj,:) == 1);
            N_ON_2 = sum(Results.urec(2,jj,:) == 1);
            
            N_on(jj) = max(N_ON_1,N_ON_2);
            disp(num2str(jj));
        end
        % store the sets N_on and N_off
        OutputStruct.sets_N_on{ii} = N_on;
        OutputStruct.sets_N_off{ii} = Params.n_app - N_on;
        
        disp(['%%%% File ',num2str(count),' data saved! %%%%']);
        count = count + 1;
        
    end
    
    % save the output struct
    save('VAGGELIS/saved_data/baseline/result/data_sets_N_ON.mat',...
        'OutputStruct');
end

%% in this part of the code the mean and standard dev analysis is performed
if(bParts(2))
    % load the already computed data-set
    load('VAGGELIS/saved_data/baseline/result/data_sets_N_ON.mat',...
        '');
    
    % define the matrices out from the loaded struct
    matrix_P_agg = cell2mat(OutputStruct.aggregate_power);
    matrix_N_on = cell2mat(OutputStruct.sets_N_on.').';
	matrix_N_off = cell2mat(OutputStruct.sets_N_off.').';
    
    % compute means for P_agg and N_on,N_off sets
    %   
    [mmean,ddevs] = mean_and_deviations(matrix_P_agg);
    OutputStruct.P_agg.mean = mmean;
    OutputStruct.P_agg.dev = ddevs.dev;
    OutputStruct.P_agg.dev_upper = ddevs.upper_bound;
    OutputStruct.P_agg.dev_lower = ddevs.lower_bound;
    
    [mmean,ddevs] = mean_and_deviations(matrix_N_on);
    OutputStruct.N_on.mean = mmean;
    OutputStruct.N_on.dev = ddevs.dev;
    OutputStruct.N_on.dev_upper = ddevs.upper_bound;
    OutputStruct.N_on.dev_lower = ddevs.lower_bound;
    
    [mmean,ddevs] = mean_and_deviations(matrix_N_off);
    OutputStruct.N_off.mean = mmean;
    OutputStruct.N_off.dev = ddevs.dev;
    OutputStruct.N_off.dev_upper = ddevs.upper_bound;
    OutputStruct.N_off.dev_lower = ddevs.lower_bound;
    
    % save the (modified) output struct
    save('VAGGELIS/saved_data/baseline/result/data_sets_N_ON.mat',...
        'OutputStruct');
    
    % plot the mean and standard devs for aggregate power and set N_on
    plot_stdevs(OutputStruct,0);
    plot_stdevs(OutputStruct,1);
end

%% in this part define a basic ARX model between P_agg and N_on set
if(bParts(3))
    % load the already computed data-set
    load('VAGGELIS/saved_data/baseline/result/data_sets_N_ON.mat');
    
    % prepare the data for the 'iddata' call
    %
    u_input = cell(1,OutputStruct.N_E);
    y_output = cell(1,OutputStruct.N_E);
    for ii = 1:OutputStruct.N_E
        u_input{ii} = OutputStruct.sets_N_on{ii}./OutputStruct.n_app;
        y_output{ii} = OutputStruct.aggregate_power{ii}.'/...
            OutputStruct.rated_power{ii};
    end
    T_period = OutputStruct.period;
    
    % create multiexperiment data
    %
    OutputStruct.data = iddata(y_output(1:OutputStruct.N_E/2),...
        u_input(1:OutputStruct.N_E/2),T_period);
    OutputStruct.validation = iddata(y_output(OutputStruct.N_E/2+1:end),...
        u_input(OutputStruct.N_E/2+1:end),T_period);
    % detrend to remove non-zero mean
    OutputStruct.data = detrend(OutputStruct.data);
    % estimate if there is a delay in the system
    temp_delay = delayest(OutputStruct.data);
    
    % construct ARX model
    %
    % estimate the orders for the ARX system
    NN = struc(1:3,1:3,temp_delay);
    V_arx = arxstruc(OutputStruct.data,OutputStruct.validation,NN);
    OutputStruct.order_arx = selstruc(V_arx,0);
    % overwrite the orders to (2,2,0), because the fit is almost as good,
    % and complexity is reduced
    OutputStruct.order_arx = [2,2,temp_delay];
    
    % ARX model (with options)
    %
    opt = arxOptions('Focus','stability');
    OutputStruct.arx_sys = arx(OutputStruct.data,OutputStruct.order_arx,...
        opt);
    
    [y_est,fit,~] = compare(OutputStruct.validation,OutputStruct.arx_sys);
    OutputStruct.y_est = y_est;
    OutputStruct.fit = sum(cell2mat(fit.'))/(OutputStruct.N_E/2);
    disp(['Fit estimate ',num2str(OutputStruct.fit)])
    
    % estimate the transfer function
    %
    num = flip(OutputStruct.arx_sys.b);
    den = flip(OutputStruct.arx_sys.a);
    Ts = OutputStruct.period;
    
    OutputStruct.G_est = tf(num,den,Ts);
    % [pp,zz] = pzmap(OutputStruct.G_est);
    
    % save the (modified) output struct
    save('VAGGELIS/saved_data/baseline/result/data_sets_N_ON.mat',...
        'OutputStruct');
end
end

function [mmean,ddevs] = mean_and_deviations(my_matrix)

% compute mean
mmean = mean(my_matrix);

% calculate covariance matrix and deviations for the data
%
Sigma = cov(my_matrix);
vec_dev = sqrt(diag(Sigma));
ddevs.dev = vec_dev.';

ddevs.upper_bound = mmean + ddevs.dev;
ddevs.lower_bound = mmean - ddevs.dev;

end

function plot_stdevs(OutputStruct,sCase)
switch sCase
    case 0
        upper = OutputStruct.P_agg.dev_upper;
        lower = OutputStruct.P_agg.dev_lower;
        mmean = OutputStruct.P_agg.mean;
        t_plot_handles = 'P_{agg} [kW]';
        color_linestyle = '-b';
    case 1
        upper = OutputStruct.N_on.dev_upper;
        lower = OutputStruct.N_on.dev_lower;
        mmean = OutputStruct.N_on.mean;
        t_plot_handles = 'N_{on} [devices]';
        color_linestyle = '*b';
    case 2
        upper = OutputStruct.N_off.dev_upper;
        lower = OutputStruct.N_off.dev_lower;
        mmean = OutputStruct.N_off.mean;
        t_plot_handles = 'N_{off} [devices]';
        color_linestyle = '*b';
end

% plot the coveriance (standard devs.)
figure;
hold on;
% plot shaded area around Results_comparison.Prec
fill([OutputStruct.vTime./3600,fliplr(OutputStruct.vTime./3600)],...
    [upper,fliplr(lower)],[0.9,0.9,0.9],...
    'linestyle','none');
plot(OutputStruct.vTime./3600,mmean,color_linestyle);
hold off;
grid on;
xlabel('Time [h]'); 
ylabel(t_plot_handles);
legend(['Standard dev.'],['mean ',t_plot_handles],'Location','Best');
end