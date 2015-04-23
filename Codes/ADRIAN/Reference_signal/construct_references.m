%% ========================================================================
%  function construct_reference
%  by Adrian Martinez Gomez
%  December 2014
%
%  Output: y_out is a sturct, with the variable part of the aggregate power
%    (given as a cell array), and the different frequencies and amplitudes
%
%  Purpose:
%   * construct the reference signal for the variable part of the aggregate
%   power (aka, what we would like to track)
%  ========================================================================
function y_out = construct_references(bNormalize,Ts);
tic;
%% 1) initialize the signals
% adjust these 4 vectors to obtain the desired signals; they are crucial
% for the calculation of the (C)QP problem

% for 2) (define standard dev. and M)
standard_devs = (1:1:10).*10^(-3);
vMs = 100.*(2:2:20);

% for 3)
frequencies = 0.0002:0.0002:0.0016;
periods = 1./frequencies;
amplitudes = 0.1:0.1:1.0;

fprintf('\n');
disp(['The period T ranges from ',num2str(periods(1)),...
    ' [seconds] to ',num2str(periods(end)),' [seconds]']);
disp(['The frequency f ranges from ',num2str(frequencies(1)),...
    ' [Hz] to ',num2str(frequencies(end)),' [Hz]']);

% define the number of signals (aka, how many reference signals are
% goung to be generated)
number_signals_normal = length(standard_devs)*length(vMs);
number_signals_sinusoids = length(frequencies)*length(amplitudes);
number_signals = number_signals_normal + number_signals_sinusoids;

fprintf('\n');
disp(['We are going to make ',num2str(number_signals),...
    ' different ref signals']);
disp([' => Doing cumsum of normal distributed random variables:  ',...
    num2str(number_signals_normal),' different ones']);
disp([' => Doing sinusoids with different freq. and amplitudes:  ',...
    num2str(number_signals_sinusoids),' different ones']);

% pre-allocate the output cell-array
y = cell(1,number_signals);
for ii = 1:number_signals
    y{ii} = nan(24*3600+1,1);
end
% pre-allocate the output struct
y_out.signal = cell(1,number_signals);
y_out.frequencies = cell(1,number_signals);
y_out.periods = cell(1,number_signals);
y_out.amplitudes = cell(1,number_signals);

%% 2) generate some of them as sum of normal signals
%  (ala Callaway's paper, where he defines the input)

% iterate through: number of added rand. numb. M ("frequency"), standard
% dev. for the rand. num. standard_dev ("amplitude").
%
count = 1;

for ii = 1:length(standard_devs)
    v = standard_devs(ii).*randn(2*24*3600,1);
    u = nan(24*3600+1,1);
    for jj = 1:length(vMs)
        
        for kk = 1:length(u)
            u(kk) = sum(v(kk:(kk + vMs(jj))));
        end
        % save the result into the
        y{count} = u;
        count = count + 1;
    end
end

%% 3) generate some of them as basic sinusoids
x = 1:(24*3600+1);

for ii = 1:length(frequencies)
    for jj = 1:length(amplitudes)
        % store the sinusoid
        y{count} = amplitudes(jj)*sin(2*pi*frequencies(ii)*x);
        
        % store the frequency information
        y_out.frequencies{number_signals_normal+count} = frequencies(ii);
        y_out.periods{number_signals_normal+count} = periods(ii);
        y_out.amplitudes{number_signals_normal+count} = amplitudes(jj);
        
        count = count + 1;
    end
end

%% 4) Post-process the signals
for ii = 1:number_signals
    % downsample the signal according to Ts
    y_out.signal{ii} = downsample(y{ii},Ts);
    % limit the signal to be between -1 and 1
    y_out.signal{ii} = max(min(y_out.signal{ii},1),-1);
end

fprintf('\n');
disp([' => After the downsampling of ',num2str(Ts),' [seconds]:']);
disp(['The period T ranges from ',num2str(periods(1)/Ts),...
    ' [steps] to ',num2str(periods(end)/Ts),' [steps]']);
disp(['The frequency f ranges from ',num2str(frequencies(1)*Ts),...
    ' [1/steps] to ',num2str(frequencies(end)*Ts),' [1/steps]']);

% if needed, normalize the signal
%
% NOTE: this doesn't seem practical at all for the problem at hand
%  The signals resulting from this can't be used.
if(bNormalize)
    for ii = 1:number_signals
        y_out.signal{ii} = y_out.signal{ii}./norm(y{ii},2);
    end
end

% show the running time for the script
sim_time = toc;
fprintf('\n');
disp([' Simulation time for the script: ',num2str(sim_time),' seconds']);

%% 5) How to call this script (kept here as reference)
%
% y_out = construct_references(0,10);
% ref_signals_1000EWH = y_out;
% save('archived_data/reference_signals/1000EWH/ref_signals_1000EWH.mat',...
%     'ref_signals_1000EWH');

end