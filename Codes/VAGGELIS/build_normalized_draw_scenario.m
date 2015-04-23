function WaterDrawScenario = build_normalized_draw_scenario(Params)
%% ========================================================================
%  Build a water draw scenario based on predefined profile for 1 day
%  Stephan Koch
%  EEH - Power Systems Laborator, ETH Zurich
%  September 14, 2010
%
%  Modified by Evangelos Vrettos
%  October 2014
%  ========================================================================

% input: t_sample in s, t_sim in s
% usage: t_sample in h, t_sim in h --> conversion:
t_sample = Params.t_sample / 3600;
t_sim = Params.t_sim /3600;

WaterDrawProfiles.values = [0.0129 0.0091 0.0091 0.0098 0.0144 0.0263 0.0491 0.054 0.054 0.0538 0.0515 0.0492 0.047 0.043 0.042 0.043 0.0492 0.061 0.0695 0.0685 0.0623 0.0528 0.0435 0.025]'; % from LBNL
% WaterDrawProfiles.values = [0.00875 0.00875 0.00875 0.00875 0.00875 0.01 0.075 0.075 0.065 0.065 0.065 0.046 0.046 0.037 0.037 0.037 0.037 0.063 0.063 0.063 0.063 0.051 0.051 0.00875]';      % from ASRHE

n_hh                                    = Params.n_app;
daily_profile                           = WaterDrawProfiles.values;
n_draw_dist                             = Params.n_draw_dist;
n_draw_dist_param                       = Params.n_draw_dist_param;
shower_dist                             = Params.shower_dist;
shower_dist_param                       = Params.shower_dist_param;
mediumLoad_dist                         = Params.mediumLoad_dist;
mediumLoad_dist_param                   = Params.mediumLoad_dist_param;
shortLoad_dist                          = Params.shortLoad_dist;
shortLoad_dist_param                    = Params.shortLoad_dist_param;
flow_rate_shower_dist                   = Params.flow_rate_shower_dist;
flow_rate_shower_dist_param             = Params.flow_rate_shower_dist_param;
flow_rate_mediumLoad_dist               = Params.flow_rate_mediumLoad_dist;
flow_rate_mediumLoad_dist_param         = Params.flow_rate_mediumLoad_dist_param;
flow_rate_shortLoad_dist                = Params.flow_rate_shortLoad_dist;
flow_rate_shortLoad_dist_param          = Params.flow_rate_shortLoad_dist_param;

WaterDrawScenario = construct_water_draw_scenario(Params,n_hh, ...
    daily_profile,n_draw_dist,n_draw_dist_param,shower_dist,...
    shower_dist_param,mediumLoad_dist,mediumLoad_dist_param,...
    shortLoad_dist,shortLoad_dist_param,flow_rate_shower_dist,...
    flow_rate_shower_dist_param,flow_rate_mediumLoad_dist,...
    flow_rate_mediumLoad_dist_param,flow_rate_shortLoad_dist,...
    flow_rate_shortLoad_dist_param,t_sample,t_sim);

%% Recycle Bin
% n_hh                         = Params.n_app;
% daily_consumption_dist       = 'normal';
% daily_consumption_dist_param = [200,20]%[1, 0];
% daily_profile                = WaterDrawProfiles.values;
% n_draw_dist                  = Params.n_draw_dist;
% n_draw_dist_param            = Params.n_draw_dist_param;
% long_draw_dist               = Params.long_draw_dist;
% long_draw_dist_param         = Params.long_draw_dist_param;
% short_draw_dist              = Params.short_draw_dist;
% short_draw_dist_param        = Params.short_draw_dist_param;
% flow_rate_dist               = Params.flow_rate_dist; 
% flow_rate_dist_param         = Params.flow_rate_dist_param;

% WaterDrawScenario = construct_water_draw_scenario_2(n_hh, daily_consumption_dist, daily_consumption_dist_param, daily_profile, n_draw_dist, n_draw_dist_param, long_draw_dist, long_draw_dist_param, short_draw_dist, short_draw_dist_param, flow_rate_dist, flow_rate_dist_param, t_sample, t_sim);

% Minimum value of flow rate is zero
% WaterDrawScenario.flow_rates = max(0, WaterDrawScenario.flow_rates);

% Total daily water draw
% WaterDrawScenario.total_draw = random(Params.total_draw_dist,Params.total_draw_mean,Params.total_draw_std,Params.n_app,1);


