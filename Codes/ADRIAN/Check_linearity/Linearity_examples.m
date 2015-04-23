%% ========================================================================
%  function Linearity_examples
%  by Adrian Martinez Gomez
%  February 2015
%
%  Purpose:
%  Plot some examples to check the linearity of the system regarding
%  Input-Output: g(u_1+u_2) = g(u_1) + g(u_2) = y_1 + y_2
%
%  It seems to work best for positive inputs: u_1,u_2 > 0 !
%  ========================================================================
function Linearity_examples();
load('archived_data/SPV/1000EWH/Data/SUM_AVERAGED_DATA_1000_EWH_18022015.mat','zdata');

% indeces for hour 7 for different amplitudes (positive inputs)
idx_03_pos_7 = 134;
idx_06_pos_7 = 541;
idx_09_pos_7 = 961;
idx_03_neg_7 = [128,133,140];

% indeces for hour 22 for different amplitudes (negative inputs)
idx_03_neg_22 = 82;
idx_06_neg_22 = 517;
idx_09_neg_22 = 927;

%%
% 0.3 + 0.3 ?= 0.6 (@ hour 7)
idx{1} = idx_03_pos_7; idx{2} = idx_03_pos_7; idx{3} = idx_06_pos_7;
Linearity_examples_plotting(zdata,idx);
close;

% 0.3 + 0.6 ?= 0.9 (@ hour 7)
idx{1} = idx_03_pos_7; idx{2} = idx_06_pos_7; idx{3} = idx_09_pos_7;
Linearity_examples_plotting(zdata,idx);
close;

% -0.3 + 0.6 ?= 0.3 (@ hour 7)
idx{1} = idx_03_neg_7(3); idx{2} = idx_06_pos_7; idx{3} = idx_03_pos_7;
Linearity_examples_plotting(zdata,idx);

% -0.3 + -0.3 ?= -0.6 (@ hour 22)
idx{1} = idx_03_neg_22; idx{2} = idx_03_neg_22; idx{3} = idx_06_neg_22;
Linearity_examples_plotting(zdata,idx);
close;

% -0.3 + -0.6 ?= -0.9 (@ hour 22)
idx{1} = idx_03_neg_22; idx{2} = idx_06_neg_22; idx{3} = idx_09_neg_22;
Linearity_examples_plotting(zdata,idx);
close;

end

function Linearity_examples_plotting(Data,idx);

figure;
subplot(2,1,1);hold on;
title(['Starting @ Hour ',idxToHour(idx{1}),' of the day']);
plot(Data.u{idx{1}}(:,1) + Data.u{idx{2}}(:,1));
plot(Data.u{idx{3}}(:,1),'r');
legend('u_{1} + u_{2}','u_{3}','Location','Best');
hold off;

subplot(2,1,2);hold on;
plot(Data.y{idx{1}} + Data.y{idx{2}});
plot(Data.y{idx{3}},'r');
legend('y(u_1_{1}) + y(u_{2})','y(u_{3})','Location','Best');
hold off;

end