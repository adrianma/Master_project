%% ===================================================================
%  function define_ref_signal
%  by Adrian Martinez Gomez
%  November 2014
%
%  Purpose:
%   * compute the reference signal (\eta based) for the step examples
%   performed.
%   * this will allow us to compute transfer function estimate in other
%   functions/scripts.
%  ===================================================================
function ref_signal = define_ref_signal(Params,Results);

vTime = Params.t_init:Params.t_sample:Params.t_sim;
N = length(vTime);
low = find(vTime == Results.Event(1));
high = find(vTime == Results.Event(end));

ref_signal = [zeros(1,(low-1)),Results.etas(1).*ones(1,high-low+1),...
    zeros(1,N-high)];

end