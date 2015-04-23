function Data_new = Adapt_experiment_number(Data);
% 0) change the experiment names for 'Data'
for ii = 1:size(Data.ExperimentName,1)
    Data.ExperimentName{ii} = strcat('Exp',...
        num2str(size(Data.ExperimentName,1) + ii));
end

% pass the output argument
Data_new = Data;

end