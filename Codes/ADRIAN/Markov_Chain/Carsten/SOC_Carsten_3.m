%% ========================================================================
%  function SOC_Carsten_3
%  by Adrian Martinez Gomez
%  April 2015
%
%  Purpose:
%   Calculate the SOC definition 3 from Carsten Heinrich's Semester project
%   report.
%   The formula has the form:
%       SOC = min{1,max{0,{sum{T_i - T_min}}/{N*(T_max - T_min)}}}
%
%  ========================================================================
function [SOC] = SOC_Carsten_3(idx_file);
%% 1) Load the file
%
path_file = ['archived_data/baseline/1000EWH/res_same/sim_1000EWHs_runNumber',...
        num2str(idx_file),'.mat'];
load(path_file);

%% 2) Compute the SOC
%
% length of time vector simulation
N_sim = length(Results_comparison.Prec);
% pre-allocate the SOC
SOC = cell(1,N_sim);
for ii = 1:N_sim
    % this two commands correspond to the denominator of the fraction
    numerator = squeeze(Results_comparison.xrec(2:end-1,ii,:)).';
    numerator = sum((numerator - repmat(squeeze(Params.T_min1),1,10)).');
    
    % this calculates the denominator of the equation
    denominator = Params.n.*squeeze(Params.T_max1 - Params.T_min1);
    denominator = denominator.';
    
    % finally, the total definition of the 
    SOC{ii} = min(1,max(0,numerator./denominator));
end
% convert from cell array to matrix
SOC = cell2mat(SOC.');

end