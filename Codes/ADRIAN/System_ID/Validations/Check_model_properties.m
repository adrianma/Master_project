%% ========================================================================
%  function Check_model_properties
%  by Adrian Martinez Gomez
%  February 2015
%
%  Purpose:
%   * checks stability, controllability and observability for a state-space
%   system ss, and shows on the console the properties.
%  ========================================================================
function Check_model_properties(ss)

fprintf('\n');

% stability 
if(isstable(ss))
    disp(' ==> The discrete system is stable or asymptotically stable');
else
    disp(' ==> The discrete system is NOT stable');
end

% observability
bObs = rank(obsv(ss));
if(bObs >= rank(ss.a))
    disp(' ==> Model observable');
else
    disp(' ==> Model NOT observable');
end

% controllability
bCon = rank(ctrb(ss));
if(bCon >= rank(ss.a))
    disp(' ==> Model controllable');
else
    disp(' ==> Model NOT controllable');
end

% eigenvalues of the ss system
disp(' ==> Eigenvalues of the state-space system: |lambda| < 1')
abs(eig(ss.a).')

fprintf('\n');

end