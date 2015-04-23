function [uu] = Test_inputs(Sigma_idx);
%% Definitions
%
% I tested the sigmas, and the cases with sigma = {1,0.001} were useless. 
% (the first one because of too big values (jumping btw 1 and -1), the
% second one because of a too low SNR)
if(0)
    sigmas = [1,0.1,0.01,0.001];
else
   sigmas = 0.01:0.01:0.08;
end
Ms = [5,15,25,35,55,75,100,200];

%% Calculate the inputs
uu = cell(length(sigmas),length(Ms),10);
for ii = 1:length(sigmas)
    for jj = 1:length(Ms)
        for kk = 1:10
            uu{ii,jj,kk} = Generate_input(sigmas(ii),Ms(jj));
            close;
        end
    end
end

%% Generate the plots
%
cmap = hsv(10);
for ii = 1:length(Ms)
    
    figure;
    
    hold on;
    for kk = 1:10
        plot(uu{Sigma_idx,ii,kk},'Color',cmap(kk,:));
    end
    hold off;
    
    title([' \sigma=',num2str(sigmas(Sigma_idx)),' and M=',...
        num2str(num2str(Ms(ii)))]);
end
end