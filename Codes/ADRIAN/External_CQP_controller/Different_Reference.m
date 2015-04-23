%% ========================================================================
%  function Different_Reference
%  by Adrian Martinez Gomez
%  March 2015
%
%  Inputs:
%   *opt_ref: choose 1,2,... for different reference signals to test with
%   the controller
%
%  Purpose:
%   Generate an alternative reference signal.
%   Ideally it should be piecewise constant.
%   Try sinusoids as well (Maryam's suggestion)
%
%  Problems/Questions:
%   My models are time varying. How to generate an "artificial" signal for
%   these?
%  ========================================================================
function [RefSignal] = Different_Reference(opt_ref);

% number of filling zeros
N_zeros = 20;
% pre-allocate the reference signal as nan
RefSignal = nan(1,460);
% first N_zeros steps are zero
RefSignal(1:N_zeros) = zeros(1,20);

switch opt_ref
    case 1        
        % this is a trapezoid, close to a pulse, but different
        %
        % do 30 steps of ramp going "up"
        x_ramp_up = 1:1:30;
        RefSignal((N_zeros+1):50) = 0.4/30.*x_ramp_up;
        % do 30 steps with a constant value
        RefSignal(51:80) = 0.4;
        % do 70 steps of ramp going "down"
        x_ramp_down = 1:1:70;
        RefSignal(81:150) = 0.4 - 0.4/70.*x_ramp_down;
        % recover to zero reference signal
        RefSignal(151:end) = zeros(1,length(RefSignal)-150);        
    case 2
        % do 200 steps of ramp going "up"
        x_ramp_up = 1:1:200;
        RefSignal((N_zeros+1):220) = 0.2/200.*x_ramp_up;
        % do remaining steps with a constant value
        RefSignal(221:end) = 0.2;
    case 3
        % generate a quadratic that starts at jump at 0.2, and after one
        % oscillation remains at a constant value of 0.2 for more steps
        x = -1:0.05:1;
        RefSignal((N_zeros+1):61) = 0.3.*ones(size(x)) - 0.12.*x.^2;
        RefSignal(62:130) = 0.2.*ones(1,130-62+1);
        RefSignal(131:end) = zeros(1,length(RefSignal)-130);
    case 4
        % best hour (=idx) for this experiment are 16,12
        %
        x_sinusoid = 1:1:440;
        % generate 2 sinusoidal oscillations
        RefSignal((N_zeros+1):end) = 0.5.*abs(sin(2*pi*x_sinusoid./300));
        % fill with zeros after the 2 oscillations
        RefSignal((N_zeros+300):end) = zeros(size(RefSignal(...
            (N_zeros+300):end)));
    case 5
        
end