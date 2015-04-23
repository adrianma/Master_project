%% ===================================================================
%  function Compare_bad_fits
%  by Adrian Martinez Gomez
%  February 2015
%
%  Inputs:
%
%  Purpose:
%   * Compare a vector of "bad" fits (exp_num_bad) between model (ss) and 
%   real system data (Data).
%  ===================================================================
function Compare_bad_fits(ss,Data,Validation,Fit_vector);
%% 1) Convert the "validation" bad idx to "data" bad idx

% get the "bad" indeces of validation
val_idx_bad = find(Fit_vector < 0);
N_bad = length(val_idx_bad);

disp(['There are ',num2str(N_bad),' experiments where Fit goes badly']);

% retrieve and convert to "bad" indeces of data
data_idx_bad = nan(1,N_bad);
for ii = 1:N_bad
    data_idx_bad(ii) = str2num(Validation.ExperimentName{...
        val_idx_bad(ii)}(4:end));
end

%% 2) Generate the comparison data and get the plots
event_hour = nan(1,N_bad);
for ii = 1:N_bad
    
    event_hour(ii) = idxToHour(data_idx_bad(ii));
    u_value = Data.u{data_idx_bad(ii)}(25);
    
    [yy,~,~] = compare(getexp(Data,data_idx_bad(ii)),ss);
    
    figure;
    subplot(2,1,1);
    plot(Data.u{data_idx_bad(ii)});
    title(['Event hour = ',num2str(event_hour),'h and |u|=',...
        num2str(u_value)]);
    legend('u','Location','Best');
    xlabel('Step Index k');
    ylabel('Temperature set-point variation');
    grid on;
    
    subplot(2,1,2);
    hold on;
    plot(Data.y{data_idx_bad(ii)});
    plot(yy.OutputData,'r');
    legend('y_{REF}','y_{LM}','Location','Best');
    xlabel('Step Index k');
    ylabel('Normalized output - aggregate power - averaged baseline');
    grid on;
    
    hold off;
end

%% 3) Get the number of ocurrences of bad idx on the terminal printed out
vHours = [4,7,12,16,18,20,22];
[vOcurrences,~] = hist(event_hour,vHours);

disp('The issues happen at hours: ');
for ii = 1:length(vHours)
    disp(['   ',num2str(vOcurrences(ii)),' @hour ',num2str(vHours(ii))]);
end
fprintf('\n');

end