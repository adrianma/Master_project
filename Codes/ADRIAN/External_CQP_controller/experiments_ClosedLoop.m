%% ========================================================================
%  script experiments_ClosedLoop.m
%  by Adrian Martinez Gomez
%  April 2014
%
%  Purpose:
%   * This script calculates closed-loop examples similar to the ones one
%   can find on the thesis in Chapter 5.3
%   * Of course it can be changed to test othes state-space systems or
%   varying the time-horizon, the terminal quadratic cost, etc
%  ========================================================================

%% (0) Load the necessary data, models, and define input structs
%
% load the necessary the data (experimental and validation)
load('archived_data/SPV/1000EWH/Data_same/15032015_same_SPV.mat');

% load the ARMAX/state-space models
load('archived_data/SPV/1000EWH/Models_same/15032015_models_ss_same.mat');

% define struct for running the simulations
%
% prediction time horizon for the MPC
rRefining.N_horizon = 30;
% quadratic terminal cost
rRefining.p_terminal = 20;
% ignore this bHalf value and keep it at zero
rRefining.bHalf = 0;
% upper boundary value ub (u_max) for optimization
rRefining.u_max = 0.95;
% lower boundary value lb (u_min) for optimization
rRefining.u_min = 0.0;
% select to 1 to check "artificial" (user-defined) reference signals
rRefining.bArtificial = 1;
% select 1 to test with the real system
rRefining.bRealSystem = 1;
% many choices for the cost function. as it is written it will be tracking
rRefining.CostFunctionCoeff = [1,0,0,0];
% we don't need to save the figures/plots and data
rRefining.bSave = 0;

%% (1) Perform the examples
%
% loop through the 4 possible "artificial" reference signals
for ii = 1:4
    
    % select the artificial signal type
    rRefining.ArtificialType = ii;
    
    % =================================================================== %
    rRefining.idx = 170;
    disp(['   Hour of the simulation',num2str(idxToHour(rRefining.idx))]);
    % do the run
    main_CL_simulation(ss_16_pos,zData_same_2,rRefining);
    
    % relax the boundaries
    rRefining.u_min = -0.95;
    disp('   relax the boundaries');
    % do the run
    main_CL_simulation(ss_16_pos,zData_same_2,rRefining);
    
    % =================================================================== %
    
    rRefining.idx = 210;
    disp(['   Hour of the simulation',num2str(idxToHour(rRefining.idx))]);
    % select the artificial signal type
    % relax the boundaries
    rRefining.u_min = 0;
    % do the run
    main_CL_simulation(ss_20_pos,zData_same_2,rRefining);
    
    % relax the boundaries
    rRefining.u_min = -0.95;
    disp('   relax the boundaries');
    % do the run
    main_CL_simulation(ss_20_pos,zData_same_2,rRefining);
    
end