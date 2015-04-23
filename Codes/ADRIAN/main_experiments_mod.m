%% ========================================================================
%  function main_experiments
%  by Adrian Martinez Gomez
%  November 2014
%
%  Purpose:
%   * do the experiment for multiple amplitudes as well as different files
%   (= different baseline trajectories)
%  ========================================================================
function main_experiments_mod(Method,hour,n_app);

% Method =                        'ProbSwitching';

% choose the folder for the method (PS:prob. switching,
% SPV:set-point variation)
if(strcmp(Method,'ProbSwitching'))
    tFolderMethod = 'PS';
elseif(strcmp(Method,'SetPointVariation'))
    tFolderMethod = 'SPV';
end
folder_files = dir(['VAGGELIS/saved_data/',tFolderMethod,'/0/*.mat']);

% hour =                          12;
input_type =                    0;
N_E =                           length(folder_files);
bDifferentFiles =               1;
% choose amplitude of the pulse
Amplitudes.eta_pos_low =        0.05;
Amplitudes.eta_pos_high =       0.3;
Amplitudes.eta_neg_low =        -0.05;
Amplitudes.eta_neg_high =       -0.7;
% change duration of the pulse [minutes]
Amplitudes.min_duration =       5;
Amplitudes.max_duration =       15;
NAnalysis =                     1;

% changes in step (how do variations go)
%
interval_change = [0,0];
count = 1;
while(count <= 5)
    
    disp(['%------------------- Amplitudes ',num2str(count+1),...
        ' -------------------%']);
    disp(['   eta_+ ~ Unif(',num2str(Amplitudes.eta_pos_low),',',...
        num2str(Amplitudes.eta_pos_high),')']);
    disp(['   eta_- ~ Unif(',num2str(Amplitudes.eta_neg_low),',',...
        num2str(Amplitudes.eta_neg_high),')']);
    
    % simulate the system
    tests_ext_controller_OL_hourly(Method,hour,input_type,N_E,...
        bDifferentFiles,Amplitudes,NAnalysis,n_app);
    
    % change the amplitude range
    Amplitudes.eta_pos_low = Amplitudes.eta_pos_low + interval_change(1);
    Amplitudes.eta_pos_high = Amplitudes.eta_pos_high + interval_change(1);
    Amplitudes.eta_neg_low = Amplitudes.eta_neg_low - interval_change(2);
    Amplitudes.eta_neg_high = Amplitudes.eta_neg_high - interval_change(2);
    
    % also, remember to update the folder number and the counter
    NAnalysis = NAnalysis + 1;
    count = count + 1;
    
end
end

%% ===================================================================
%  function tests_set_point_variation
%  by Adrian Martinez Gomez
%  October 2014
%   
%  Inputs:
%   *Method: string containing the method to be tested (ProbSwitching,
%   SetPointVariation, etc)
%   *bSimulate: boolean, set to 1 if want to simulate
%   *N_E: number of experiments
%   *NAnalysis: giving the new folder name for the simulations to be
%   stored. Set to "0" if no new simulation set needed.
%
%  Purpose:
%   *do an hour-long experiment varying the input signal
%  ===================================================================
function tests_ext_controller_OL_hourly(Method,hour,input_type,N_E,...
    bDifferentFiles,Amplitudes,NAnalysis,n_app);

if(NAnalysis==0)
    error('Choose a value for NAnalysis different from 0');
end

global eta_in Events event_hour bInput;
% choose 1 if input is not a pulse
bInput = input_type;

% define the simulated time (=1 day)
Params.t_sim = 24*3600;
Params.t_sample = 10.0;

event_hour = hour*3600;

% construct the N_E input signals
% N_E: number of experiments

% choose the folder for the method (PS:prob. switching,
% SPV:set-point variation)
if(strcmp(Method,'ProbSwitching'))
    tFolderMethod = 'PS'; 
elseif(strcmp(Method,'SetPointVariation'))
    tFolderMethod = 'SPV';
end

% eta_in gives the amplitude of the pulse.
%   uniformly between -0.3 and 0.2
temp_neg = random('unif',Amplitudes.eta_neg_high,Amplitudes.eta_neg_low,...
    1,floor(N_E/2));
temp_pos = random('unif',Amplitudes.eta_pos_low,Amplitudes.eta_pos_high,...
    1,ceil(N_E/2));
temp = [temp_pos,temp_neg];
vIdx = randperm(N_E);
eta_in = temp(vIdx);

% define duration (t_duration) time for the pulses
%   duration: uniformly between 5 and 15 minutes
%
% NOTE: to be sure that the event is in [step] units select the
% duration as done below.
% [min]->[sec] would be factor 60
Params.start = event_hour*ones(1,N_E);
Params.duration = randi([Amplitudes.min_duration*60/Params.t_sample,...
    Amplitudes.max_duration*60/Params.t_sample],1,N_E).*Params.t_sample;

Events = cell(1,N_E);
for ii = 1:N_E
    Events{ii} = Params.start(ii):1:(Params.start(ii)+Params.duration(ii));
    
    if(bInput)
        Events{ii} = Params.start(ii):1:(Params.start(ii)+(3600/10 + 1));
    end
end

%% Simulate
tests_main_for_Adrian_hourly(n_app,Method,tFolderMethod,N_E,...
    bDifferentFiles,Amplitudes,NAnalysis);

end