load('VAGGELIS/analysis/SPV/analyzed_1');
vIdx = find(DataAnalysis.vTime > 12*3600);

% =================== EXPERIMENT 1 =================== %
load('/VAGGELIS/saved_data/SPV/refined/sim_100EWHs_SetPointVariation_1_event2_eta_0005.mat');
Results0005 = Results;
Results0005_comparison = Results_comparison;
load('/VAGGELIS/saved_data/SPV/sim_100EWHs_SetPointVariation_1_event2_eta_005.mat');
Results005 = Results;
Results005_comparison = Results_comparison;

% compute the residuals
residuals_exp(1,:) = abs(Results005.Prec(vIdx) - Results0005.Prec(vIdx));
residuals_exp_comparison(1,:) = ...
    abs(Results005_comparison.Prec(vIdx) - Results0005_comparison.Prec(vIdx));

% compute the means AFTER THE EVENT (hour 12)
mresiduals_exp(1) = mean(residuals_exp(1,:));
mresiduals_exp_comparison(1) = mean(residuals_exp_comparison(1,:));

% =================== EXPERIMENT 2 =================== %
load('/VAGGELIS/saved_data/SPV/refined/sim_100EWHs_SetPointVariation_1_event2_eta_0010.mat');
Results001 = Results;
Results001_comparison = Results_comparison;
load('/VAGGELIS/saved_data/SPV/sim_100EWHs_SetPointVariation_1_event2_eta_01.mat');
Results01 = Results;
Results01_comparison = Results_comparison;

% compute the residuals
residuals_exp(2,:) = abs(Results01.Prec(vIdx) - Results001.Prec(vIdx));
residuals_exp_comparison(2,:) = ...
    abs(Results01_comparison.Prec(vIdx) - Results001_comparison.Prec(vIdx));

% compute the means AFTER THE EVENT (hour 12)
mresiduals_exp(2) = mean(residuals_exp(2,:));
mresiduals_exp_comparison(2) = mean(residuals_exp_comparison(2,:));

% =================== EXPERIMENT 3 =================== %
load('/VAGGELIS/saved_data/PS/sim_100EWHs_ProbSwitching_1_event2_eta_01.mat');
Results01 = Results;
Results01_comparison = Results_comparison;
load('/VAGGELIS/saved_data/PS/sim_100EWHs_ProbSwitching_1_event2_eta_05.mat');
Results05 = Results;
Results05_comparison = Results_comparison;
load('/VAGGELIS/saved_data/PS/sim_100EWHs_ProbSwitching_1_event2_eta_07.mat');
Results07 = Results;
Results07_comparison = Results_comparison;
load('/VAGGELIS/saved_data/PS/sim_100EWHs_ProbSwitching_1_event2_eta_09.mat');
Results09 = Results;
Results09_comparison = Results_comparison;

% compute the residuals
residuals_exp(3,:) = abs(Results05.Prec(vIdx) - Results01.Prec(vIdx));
residuals_exp_comparison(3,:) = ...
    abs(Results05_comparison.Prec(vIdx) - Results01_comparison.Prec(vIdx));

residuals_exp(4,:) = abs(Results07.Prec(vIdx) - Results01.Prec(vIdx));
residuals_exp_comparison(4,:) = ...
    abs(Results07_comparison.Prec(vIdx) - Results01_comparison.Prec(vIdx));

residuals_exp(5,:) = abs(Results09.Prec(vIdx) - Results01.Prec(vIdx));
residuals_exp_comparison(5,:) = ...
    abs(Results09_comparison.Prec(vIdx) - Results01_comparison.Prec(vIdx));

% compute the means AFTER THE EVENT (hour 12)
mresiduals_exp(3) = mean(residuals_exp(3,:));
mresiduals_exp_comparison(3) = mean(residuals_exp_comparison(3,:));

mresiduals_exp(4) = mean(residuals_exp(4,:));
mresiduals_exp_comparison(4) = mean(residuals_exp_comparison(4,:));

mresiduals_exp(5) = mean(residuals_exp(5,:));
mresiduals_exp_comparison(5) = mean(residuals_exp_comparison(5,:));

