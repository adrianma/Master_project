%% ===================================================================
%  Script runs_P_plot
%  by Adrian Martinez Gomez
%  October 2014
%
%  Purpose:
%   *Do nSim simulations for X EWHs of the Vrettos's set of scripts.
%    (remember to change n_app value for X in the script 'build_population'
%   *Store those values for later analysis.
%   *Display the mean of the aggregate power for all simulations for 1 day
%  ===================================================================
function [avP_set] = runs_P_plot(bCalculate,bPlots);
nSim = 10;

% calculate the nSim simulations for 1 day
if(bCalculate)
    for count = 1:nSim
        MainForAdrian(100,0,count,Method);
        
        if(count == 6)
            timesteps = size(Results.urec,2);
        end
    end
end

if(~exist('Params.n_app','var'))
    Params.n_app = 100;
end

% load the power
load(['VAGGELIS/saved_data/PS/sim_',num2str(Params.n_app),...
        'EWHs_ProbSwitching_1_event2_eta_-09.mat']);
mPrec(1,:) = Results_comparison.Prec;
load(['VAGGELIS/saved_data/PS/sim_',num2str(Params.n_app),...
        'EWHs_ProbSwitching_1_event2_eta_-01.mat']);
mPrec(2,:) = Results_comparison.Prec;
load(['VAGGELIS/saved_data/PS/sim_',num2str(Params.n_app),...
        'EWHs_ProbSwitching_1_event2_eta_01.mat']);
mPrec(3,:) = Results_comparison.Prec;

% Calculate average power trajectory for all times
avP_set = mean(mPrec);

vTime = 0:Params.t_sample:Params.t_sim;

%% Plot avg P_agg
if(bPlots)
    figure;
    plot(vTime./3600,avP_set,'LineWidth',3); grid on;
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 14);
    xlim([0,24]);
    xlabel('Time [h]');
    ylabel('Power [kW]');
end

end