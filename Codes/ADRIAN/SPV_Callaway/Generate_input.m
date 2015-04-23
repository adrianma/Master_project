%% ========================================================================
%  function Assess_hourly
%  by Adrian Martinez Gomez
%  March 2015
%
%  Purpose:
%  Generate an external input which is a moving sum of a white noise 
%  process.
%  This is done to emulate Callaway's paper and experimental setup.
%
%  ========================================================================
function [Input] = Generate_input(sigma,M);
% sigma = 5*10^-2;
% M = 25;
mu = 0;
vv = random('Normal',mu,sigma,1,1000);

u = nan(1,length(vv)/2);
for ii = 1:length(vv)/2
    u(ii) = sum(vv((ii+1):(ii+M)));
    
    % keep the input bounded between -1 and +1
    u(ii) = min(u(ii),1);
    u(ii) = max(u(ii),-1);
    
end

% store the input u 
Input = u;

% plot the input
figure;
title([' \sigma = ',num2str(sigma),' and M = ',num2str(M)]);
plot(u(1:361));
grid on;
xlabel('Step k');

end