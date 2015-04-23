function LFC_signal = create_LFC_signal(baseline_power, Params)
%% ========================================================================
%  Creates the LFC signal that is sent to the EWH population. It is based
%  on real LFC signal measured in Switzerland.
%
%  Evangelos Vrettos
%  June 2012
%  ========================================================================

%% Load the real LFC signal
if(0)
if Params.t_sample == 10
    load LFCsignal                      % signal from 2009
    LFC_Y = LFCsignal(18300:26940);
elseif Params.t_sample > 10
    load LFCsignal
    check = mod(Params.t_sample,10);
    if check ~= 0
        error('The simulation time step must be a multiple of 10 seconds')
    end
    
    tmp = LFCsignal(18300:26940);
    nSteps = Params.t_sample/10;
    num = floor(length(tmp)/nSteps);
    LFC_Y = NaN(num,1);
    
    i = 1;
    for k = 1:num
        LFC_Y(k) = mean(tmp(i:i+nSteps-1));
        i = i + nSteps;
    end
    LFC_Y(end+1) = LFC_Y(end);
elseif 1 || Params.t_sample < 10
    load LfcSigDec2012                  % signal from 2012
    if 0 && Params.t_sample ~= 1
        error('The simulation time step must be exactly 1 second or larger than 10 seconds')
    end
    LFC_Y = LfcSigDec2012;
end
end

%% Done by Adrian
%
% October 2014
load('VAGGELIS/LfcSigDec2012.mat');
LFC_Y = downsample(LfcSigDec2012,10);

%% Construct the LFC signal to be sent to the EWH central controller
Base_part = baseline_power;

Pband = Params.reserve_percent * baseline_power;

Variable_part = (LFC_Y/100) .* Pband;

LFC_signal = Base_part + Variable_part;
