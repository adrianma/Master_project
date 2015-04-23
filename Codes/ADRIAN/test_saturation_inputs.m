%% ===================================================================
%  function test_saturation_inputs
%  by Adrian Martinez Gomez
%  November 2014
%
%  Purpose:
%   *
%   *
%  ===================================================================
function UU = test_saturation_inputs(N_length)

%% define first input and its energy
standard_dev = 5*10^-3;
my_standard_dev = 0.5*0.0125; %0.25*5*10^2
vv = random('norm',0,my_standard_dev,1,10000);

uu_1 = nan(1,N_length);
for ii = 1:length(uu_1)
    uu_1(ii) = sum(vv(ii:ii+25));
end

% display the figure
if(0)
    figure;
    plot(uu_1);
end

% calculate energy of the signal
Energy_1 = sum(abs(uu_1).^2);

UU.input_1 = uu_1;
UU.energy_1 = Energy_1;

if(0)
    %% define 2nd input and try to find its saturation(s)
    etas_pos = [0:0.01:0.5].';
    etas_neg = [0:-0.01:-0.5].';
    out_magn_pos = zeros(length(etas_pos),1);
    out_magn_neg = zeros(length(etas_neg),1);
    
    % define time vector
    vTime = Params.t_init:Params.t_sample:Params.t_sim;
    % define frequency of the sinusoids
    omega = 0.1;
    
    % construct the input signal in this loop
    for ii = 1:length(etas_pos)
        % select the input signal as a sinusoid
        uu_2 = etas_pos(ii)*sin(omega*vTime);
    end
    
    % simulate and do analysis in this other loop
    for ii = 1:length(etas_pos)
        uu_2_temp = uu_2
        
        % SIMULATE THE SYSTEM!
        
        %  select only the second half (i.e. allow the tranisent to die down).
        %  and correlate with sin/cos signals.
        Ic = (2/N2)*cos(fsine2a*t2(N2/2:N2))'*yy([N2/2:N2],1);
        Is = (2/N2)*sin(fsine2a*t2(N2/2:N2))'*yy([N2/2:N2],1);
        out_magn_pos(ii) = sqrt(Ic^2 + Is^2)/2;
        
        fprintf('Channel:  |u| =  %g,  |y| = %g, gain = %g\n',...
            etas_pos(ii),out_magn_pos(ii), out_magn_pos(ii)/etas_pos(ii));
        
    end
end
end

%% Recycle bin
%
% low_idx = event_hour/Params.t_sample + 1;
% high_idx = (event_hour + 1*3600)/Params.t_sample + 1;
% low_event_idx = Events{ii}(1)/Params.t_sample + 1 - low_idx;
% high_event_idx = Events{ii}(end)/Params.t_sample + 1 - low_idx;
% for ii = 1:20
% figure;
% hold on;
% plot(Results{ii}.Prec - Results_comparison.Prec(low_idx:high_idx));
% line([0,361],[0,0],'Color','k','LineWidth',3)
% end