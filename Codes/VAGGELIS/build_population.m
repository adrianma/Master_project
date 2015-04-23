function Params = build_population(n_app)
%% ===================================================================
%  Build parameter structure and store it in workspace
%  Stephan Koch
%  Oct 17, 2010
%
%  Modified by Evangelos Vrettos
%  April 2012
% 
%  Modified by Evangelos Vrettos
%  October 2014
%  ===================================================================


%% Define how to vary the parameters
% Parameter variation

% Example of reasoning

% 4 Types: Mass: 50, 100, 200, 300 l

% 50 l: 2, 2.5, 3 kW,  draws:  30 ...  70 l
% 100 l: 3, 3.5, 4 kW, draws:  50 ... 100 l
% 200 l: 4, 4.5, 5 kW, draws: 100 ... 220 l
% 300 l: 4, 5, 6 kW,   draws: 150 ... 350 l
% Heating Element: 3, 4, 5, 6 kW
% Total water draws: 50 l ... 200 l

% Alpha: 0.2 ... 1 W/m^2/K

% Deadband: T_set = 55 ... 65 ?C
%           T_dead = 5 ... 15 ?C

% performance rating: outlet temp = setpoint +/- 5?C

%% Implement the parameter variation
% Number of appliances
Params.n_app = n_app;

m_categories        = [100, 200, 300, 400];                                     % [50, 100, 200, 300]; % L
P_rated_cat         = {[3, 3.5, 4], [3.5, 4, 4.5], [4, 4.5, 5], [4.5, 5, 6]};   % {[2, 2.5, 3], [3, 3.5, 4], [4, 4.5, 5], [4, 5, 6]};   %kW
total_draws_cat     = {[30, 70], [50 100], [100 220], [150, 350]};              % not used
m_cat_percentages   = [0.2, 0.3, 0.3, 0.2];                                     %[0.1, 0.25, 0.25, 0.4];
m_cat_cum_perc      = [cumsum(m_cat_percentages)-m_cat_percentages(1) 1];

rand_num        = random('unif', 0, 1, 1, n_app);
cat_idx         = cell(4,1);
m_vec           = zeros(1, 1, n_app);
P_rated_vec     = zeros(1, 1, n_app);
total_draw_vec  = zeros(1, 1, n_app);

% loop that assigns each boiler a tank size vector (m_vec), a power rating vector (P_rated_vec), and a total
% draw vector (total_draw_vec) based on randomly distributed numbers
for i = 1:length(m_categories),   %loop for each category; length() produces longest dimension of an array
    cat_idx{i} = find(rand_num > m_cat_cum_perc(i) & rand_num <= m_cat_cum_perc(i+1));
    m_vec(1, 1, cat_idx{i}) = m_categories(i);
    P_rated_vec(1, 1, cat_idx{i}) = P_rated_cat{i}(random('unid', 3, 1, length(cat_idx{i})));
    total_draw_vec(1, 1, cat_idx{i}) = random('unif', total_draws_cat{i}(1), total_draws_cat{i}(2), 1, length(cat_idx{i}));
end

% Deadband: (T_dead for the whole span, not just one direction!)
Params.T_dead = random('unif', 5, 15, 1, 1, n_app);
Params.T_set = random('unif', 55, 65, 1, 1, n_app);

% Option of setting fixed deadband for all heaters
Params.T_min1 = Params.T_set - 1/2*Params.T_dead;
Params.T_max1 = Params.T_set + 1/2*Params.T_dead;

Params.T_min2 = Params.T_set - 1/2*Params.T_dead;
Params.T_max2 = Params.T_set + 1/2*Params.T_dead;

Params.alpha = random('unif', 0.2, 1, 1, 1, n_app);  % hull heat transfer coeff. [W/m^2/K]

Params.m_categories = m_categories;

%% Set basic parameters
Params.m = m_vec;                                                           % mass [kg]
Params.V = Params.m/1000;                                                   % volume [m^3]
Params.h = 1.5*nthroot(Params.V, 3);                                        % height [m]

Params.A_cross_i = Params.V./Params.h;                                      % cross-section area [m^2]
Params.diam = sqrt(Params.A_cross_i/pi) * 2;                                % diameter [m]

Params.n = 10;                                                              % number of control elements (layers)
Params.V_i = Params.V/Params.n;                                             % volume per element
Params.d = Params.h/Params.n;                                               % distance between elements

Params.A_outer_i = Params.diam .* pi .* Params.d;                           % outer hull surface per slice

Params.A_total = Params.diam .* pi .* Params.h + 2.*Params.A_cross_i;


%% heating elements
% rated electric power
Params.P1_el = P_rated_vec;
Params.P2_el = P_rated_vec;

% Present or not? (Binary variable)
Params.element1_present = 1;        % lower heating element
Params.element2_present = 0;        % upper heating element

Params.element1_location = 3;       % the element location is 2, but we write 3 due to the fake first disk
Params.element2_location = 7;

Params.sensor1_location = 4;        % the sensor location is 3, but we write 4 due to the fake first disk
Params.sensor2_location = 8;

Params.mixing_zone1_size = 2;       % the number of layers that are directly affected by the heating element
Params.mixing_zone2_size = 1;

Params.mixing_zone1_ratio = 0.9;    % the percentage of heat that affects directly the layer where the heating element is     
Params.mixing_zone2_ratio = 1;

Params.eta = 0.95;

%% Heating element state: threshold temperatures, non-operating hours
Params.dominant_element = 2;
Params.blocking_hours1 = [0 0];
Params.blocking_hours2 = [0 0];

%% Water draw parameters
Params.n_draw_dist = 'unif';                                                % number of draws per day   
Params.n_draw_dist_param = [10 20; 20 40; 30 60; 30 60];                    % (lower limit & upper limit)
%                           100lt   200lt  300lt  400lt                        different parameters for different tank volumes
Params.shower_dist = 'normal';                                              % long water draw duration
Params.shower_dist_param = [5/60, 1/60];                                    % (mean, std) in hours
Params.mediumLoad_dist = 'normal';                                          % short water draw duration
Params.mediumLoad_dist_param = [1/60, 1/600];                               % (mean, std) in hours
Params.shortLoad_dist = 'normal';                                           % short water draw duration
Params.shortLoad_dist_param = [1/60, 1/600];                                % (mean, std) in hours
Params.flow_rate_shower_dist = 'normal';                                    % water flow rate distribution
Params.flow_rate_shower_dist_param = [8, 1]*60;                             % kg/h (mean, std) , [8, 2]*60;
Params.flow_rate_mediumLoad_dist = 'normal';                                % water flow rate distribution
Params.flow_rate_mediumLoad_dist_param = [6, 1]*60;                         % kg/h (mean, std) , [6, 2]*60;
Params.flow_rate_shortLoad_dist = 'normal';                                 % water flow rate distribution
Params.flow_rate_shortLoad_dist_param = [1, 2]*60;                          % kg/h (mean, std)

%% thermodynamic parameters
Params.m_i = Params.m/Params.n;     % mass in one element [kg]
Params.k = 0.6;                     % water heat conduction [W/m/K]
Params.c = 4185.5;                  % water heat capacity [J/kg/K]
Params.rho = 1000;                  % water density [kg/m^3]

Params.m_dot = 0.0;

%% Parameters for buoyancy effect
Params.K = 0.41;                    % von Kaman constant [-]
Params.g = 9.81;                    % [kg m/s^2]
Params.beta = 207*1e-6;             % [1/K]
Params.delta_l = Params.diam;
Params.turbMixing = 'No';           % 'Yes' cosiders turbulent mixing, 'No' doesn't

% formula: eps = (K*delta_l)^2 * sqrt(g*beta) * sqrt(d theta / d h)
%              = buoyancy_prefac * sqrt(d_theta / d h)
Params.buoyancy_prefac = (Params.K*Params.delta_l).^2 * sqrt(Params.g * Params.beta);

% testing
Params.buoyancy_prefac = 0;

Params.b = 1;       % the friction coefficient
Params.r = 0.7;     % the mixing coefficient

%% Sim parameters
Params.t_sample = 10;       % in seconds
Params.t_sim = 24*3600;

% initial condition (can be changed for starting and stopping the same sim)
% Params.x_init = repmat([14.4; 60*ones(Params.n,1); 19.7], [1,1,n_app]);
% ...(change)
Params.x_init = repmat([10; 60*ones(Params.n,1); 20], [1,1,n_app]);

Params.u_init = NaN(2,1,n_app);
for i = 1:n_app
    Params.u_init(:,:,i) = [round(rand); 0];            % random initial on/off state for lower heater, off state for upper heater
end

Params.Umi_init = repmat(zeros(Params.n,1), [1,1,n_app]);

Params.t_init = 0;

% Params.q_init is needed for the IAPSF method
for i = 1:n_app
    if Params.x_init(Params.sensor1_location,1,i) > Params.T_max1(:,:,i)
        Params.q_init(:,:,i) = -1;
    elseif Params.x_init(Params.sensor1_location,1,i) < Params.T_min1(:,:,i)
        Params.q_init(:,:,i) = 1;
    else
        Params.q_init(:,:,i) = 0;
    end
end

% original initial condition, not to be changed for simulation reset
Params.t_init_orig = Params.t_init;
Params.x_init_orig = Params.x_init;
Params.u_init_orig = Params.u_init;

Params.x_current = Params.x_init;

%% Circular flow when the element is on
Params.v_roundflow = 0.005;                                                 % [m/s]
Params.m_dot_roundflow = Params.m_i ./ Params.d .* Params.v_roundflow;      % [kg/s]

%% compute initial gradient and buoyancy in the tank
temp_grad = diff(Params.x_init(2:end-1,:,:)).*repmat(-1./Params.d, [Params.n-1, 1, 1]);
temp_grad(temp_grad < 0) = 0;
Params.buoyancy_factor = Params.buoyancy_prefac * sqrt(temp_grad);

%% Control parameters
Params.SOCdef = 2;                  % (1) SOC is defined based only on T_thermostat; (2) SOC is defined based on the temperatures of all layers
                                    % and no negative energies are allowed for each layer; (3) SOC is defined based on the temperatures of all layers and negative energies are allowed
Params.reserve_percent = 0.2;       % percentage of the baseline that is offered as reserves

end

%% Recycle Bin
% Water draw parameters
% Params.total_draw_dist = 'normal';              % daily water usage
% Params.total_draw_mean = 200;                   % in lt/day
% Params.total_draw_std = 20;                     % in lt/day
% Params.n_draw_dist = 'unif';                    % number of draws per day   
% Params.n_draw_dist_param = [10, 30]%[30, 50];            % (lower limit, upper limit)
% Params.long_draw_dist = 'normal';               % long water draw duration
% Params.long_draw_dist_param = [1/6, 1/60];      % (mean, std) in hours
% Params.short_draw_dist = 'normal';              % short water draw duration
% Params.short_draw_dist_param = [1/60, 1/600];   % (mean, std) in hours
% Params.flow_rate_dist = 'normal';               % water flow rate distribution
% Params.flow_rate_dist_param = [8, 2]*60;%[5, 2]*60;        % kg/h (mean, std)