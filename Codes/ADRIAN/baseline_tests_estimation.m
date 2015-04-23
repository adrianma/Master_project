%% ===================================================================
%  function plot_responses
%  by Adrian Martinez Gomez
%  November 2014
%
%  Inputs:
%	* n_app: number of appliances (loads;EWHs)
%   * number_experiments: # of times all the day is simulated
%      (useful for averaging later on)
%	* bEstimation: boolean, set to 1 if the files were already generated		
%		
%  Purpose:
%   * generate 'number_esperiments' files for later analysis
%   * save them for future use and analysis
%  ===================================================================

function baseline_tests_estimation(n_app,number_experiments,bEstimation)

if(~bEstimation)
    % loop to generate N_E files for further analysis
    for ii = 1:number_experiments
        %% Set up simulation
        Params = build_population(n_app);
        % precomputes the preliminary system model
        PrelimModel = precompute_system_model(Params);
        WaterDrawScenarioReal = build_normalized_draw_scenario(Params);
        % water draws required
        Params.bWaterDraw = 1;
        Results = simulate_population(Params,PrelimModel,...
            WaterDrawScenarioReal,WaterDrawScenarioReal,'NoControl');
        
        % save results in file for future analysis
        save(['VAGGELIS/saved_data/baseline/sim_',...
            num2str(n_app),'EWHs_runNumber',num2str(ii),'.mat']);
    end
end
	
if(bEstimation)
	% get the file names
	folder_files = dir(['archived_data/baseline/1000EWH/res/*.mat']);
	Nruns = length(folder_files);

	% iterate to get the aggregate power vectors
	mP_agg = nan(Nruns,8641);
	for ii = 1:Nruns
    
	    % load data-set
	    load(['archived_data/baseline/1000EWH/res/',folder_files(ii).name],'Results_comparison');
	    disp(folder_files(ii).name);
    
	    % get the aggregate power
	    mP_agg(ii,:) = Results_comparison.Prec;
	end

	% compute averaged P_agg
	mean_P_agg = mean(mP_agg);

	% plot the averaged power trajectory mean_P_agg
    Params.t_init = 0;Params.t_sample = 10;Params.t_sim = 24*3600;
	vTime = Params.t_init:Params.t_sample:Params.t_sim;
	
    figure;
	plot(vTime./3600,mean_P_agg);
	grid on;    
	xlabel('Time [h]');
	ylabel('P_{agg} [kW]');
    title(['Mean P_{agg} for N_E = ',num2str(number_experiments),...
        ' experiments']);
    
    % calculate covariance matrix and deviations
    Sigma_P = cov(mP_agg);
    vec_deviations = sqrt(diag(Sigma_P));
    vec_deviations = vec_deviations.';
    
    mean_P_upper = mean_P_agg + vec_deviations;
    mean_P_lower = mean_P_agg - vec_deviations;
    
    % plot the covariance (standard devs.)
    figure;
    hold on;
    % plot shaded area around Results_comparison.Prec
    fill([vTime./3600,fliplr(vTime./3600)],...
        [mean_P_upper,fliplr(mean_P_lower)],[0.9,0.9,0.9],...
        'linestyle','none');
    plot(vTime./3600,mean_P_agg);
    hold off;
    grid on;
    xlabel('Time [h]'); ylabel('P_{agg} [kW]');
    legend(['Standard dev.'],'mean P_{agg}','Location','Best');
    
end

end