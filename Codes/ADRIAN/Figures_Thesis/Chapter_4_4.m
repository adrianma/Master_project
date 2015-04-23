%% ========================================================================
%  function Chapter_4_4
%  by Adrian Martinez Gomez
%  March 2015
%
%  Purpose:
%  Generate Figures for Chapter's 4.4 "Selected Examples"
%  ========================================================================
function Chapter_4_4();
%% 0) Load the data and the armax models
%
load('archived_data/SPV/1000EWH/Data_same/15032015_same_SPV.mat');
load('archived_data/SPV/1000EWH/Models_same/15032015_models_armax_same.mat');

% 0.1) split the data into hours (validation data)
[val_4,val_7,val_12,val_16,val_18,val_20,val_22] = ...
    Split_hourly_data(zData_same_2);

% 0.2) split the hourly data into positive and negative data 
[val_4_pos,val_4_neg] = PosNeg_system_ID(val_4);
[val_7_pos,val_7_neg] = PosNeg_system_ID(val_7);
[val_12_pos,val_12_neg] = PosNeg_system_ID(val_12);
[val_16_pos,val_16_neg] = PosNeg_system_ID(val_16);
[val_18_pos,val_18_neg] = PosNeg_system_ID(val_18);
[val_20_pos,val_20_neg] = PosNeg_system_ID(val_20);
[val_22_pos,val_22_neg] = PosNeg_system_ID(val_22);

% the options for the 'compare' command
opt_comp = compareOptions('InitialCondition','z');

%% 1) Do the compare command
%
[yy_4_neg,ff_4_neg,~] = compare(val_4_neg,My_Models.my_mod_4_neg,opt_comp);
[yy_4_pos,ff_4_pos,~] = compare(val_4_pos,My_Models.my_mod_4_pos,opt_comp);
[yy_7_neg,ff_7_neg,~] = compare(val_7_neg,My_Models.my_mod_7_neg,opt_comp);
[yy_7_pos,ff_7_pos,~] = compare(val_7_pos,My_Models.my_mod_7_pos,opt_comp);
[yy_12_neg,ff_12_neg,~] = compare(val_12_neg,My_Models.my_mod_12_neg,opt_comp);
[yy_12_pos,ff_12_pos,~] = compare(val_12_pos,My_Models.my_mod_12_pos,opt_comp);
[yy_16_neg,ff_16_neg,~] = compare(val_16_neg,My_Models.my_mod_16_neg,opt_comp);
[yy_16_pos,ff_16_pos,~] = compare(val_16_pos,My_Models.my_mod_16_pos,opt_comp);
[yy_18_neg,ff_18_neg,~] = compare(val_18_neg,My_Models.my_mod_18_neg,opt_comp);
[yy_18_pos,ff_18_pos,~] = compare(val_18_pos,My_Models.my_mod_18_pos,opt_comp);
[yy_20_neg,ff_20_neg,~] = compare(val_20_neg,My_Models.my_mod_20_neg,opt_comp);
[yy_20_pos,ff_20_pos,~] = compare(val_20_pos,My_Models.my_mod_20_pos,opt_comp);
[yy_22_neg,ff_22_neg,~] = compare(val_22_neg,My_Models.my_mod_22_neg,opt_comp);
[yy_22_pos,ff_22_pos,~] = compare(val_22_pos,My_Models.my_mod_22_pos,opt_comp);

%% 2) Select the examples and plot them
%
figure; 
sub = cell(1,4);
sub{1} = plotting_my_comparison(val_20_pos,yy_20_pos,ff_20_pos,49,20*3600,1);
ylabel('Normalized output','FontSize',12);
sub{2} = plotting_my_comparison(val_12_neg,yy_12_neg,ff_12_neg,53,12*3600,2);
sub{3} = plotting_my_comparison(val_18_pos,yy_18_pos,ff_18_pos,48,18*3600,3);
ylabel('Normalized output','FontSize',12);
xlabel('Time [h]','FontSize',12);
sub{4} = plotting_my_comparison(val_16_neg,yy_16_neg,ff_16_neg,39,16*3600,4);
xlabel('Time [h]','FontSize',12);

fprintf('  => Done!\n');

end

function [my_subplot] = plotting_my_comparison(val,yy,ff,idx,event_hour,...
    sub_idx)

% number of zeros
N_zeros = 20;
% time vector
xx = linspace(event_hour/3600 - N_zeros*10/3600,...
    event_hour/3600 + 1,length(yy{idx}.y));

% plot both responses (including the fit-% in the tile)
my_subplot = subplot(2,2,sub_idx);
hold on;
plot(xx,val.y{idx});
plot(xx,yy{idx}.y,'r');
hold off;
grid on;
xlim([xx(1),xx(end)]);
legend('Real systems response','ARMAX model response','Location','Best');
title(['Fit = ' num2str(ff{idx}),'%'],'FontSize',12);

end