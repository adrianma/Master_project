%% Calculate the N_{ON}/N_{OFF} sets and the temperatures
folder_files = dir(['archived_data/baseline/1000EWH/res/*.mat']);
N_E = length(folder_files);

hours = [4,7,12,16,18,20,22];
N_hours = length(hours);

% path of the stored file
path_stored = 'archived_data/baseline/1000EWH/Npop_Temperatures.mat';

if(exist(path_stored,'file') == 2)
    
    N_ON = cell(N_E,1);
    N_OFF = cell(N_E,1);
    
    mTemperature = cell(N_E,1);
    stdTemperature = cell(N_E,1);
    
    Potential_Pow = cell(N_E,1);
    
    Water_Draws = cell(N_E,1);
    
    % gather the ON/OFF sets before the event happens
    for ii = 1:N_E
        
        path_file = strcat('archived_data/baseline/1000EWH/res/',folder_files(ii).name);
        load(path_file);
        
        % get the populations of the EWHs (hourly)
        for jj = 1:N_hours
            %         N_ON{ii,jj} = sum(Results_comparison.urec(1,hours(jj)*3600/10 + 1,:));
            %         N_OFF{ii,jj} = 1000 - N_ON{ii,jj};
            N_ON{ii} = sum(squeeze(Results_comparison.urec(1,:,:)).');
            N_OFF{ii} = 1000 - N_ON{ii};
        end
        
        % get the mean and std temperature of all the layers (aggregated)
        mTemperature{ii} = squeeze(mean(permute(Results_comparison.xrec(2:11,:,:),[3,2,1])));
        stdTemperature{ii} = squeeze(std(permute(Results_comparison.xrec(2:11,:,:),[3,2,1])));
        
        % matrix with the switches ON/OFF for the entire population for the
        % entire day; and vector with the power of all the population
        AA = squeeze(Results_comparison.urec(1,:,:));
        BB = squeeze(Params.P1_el).';
        
        Potential_Pow{ii} = zeros(size(BB));
        for jj = 1:8641
            Potential_Pow{ii}(jj) = sum(BB(1,find(AA(jj,:) == 0)));
        end
        
        % get water draws
        Water_Draws{ii} = sum(squeeze(Results_comparison.mdotrec).');
        
    end
    
else
    % if already computed, don't do it again, load it from file !
    load(path_stored);
end

%% Plot the power reservoir
%
% time vector
xx = linspace(0,24,8641).';

cmap = hsv(N_E);
figure;
hold on;
for ii = 1:N_E
    plot(xx,Potential_Pow{ii},'Color',cmap(ii,:));
end
hold off; grid on;
xlabel('Time [h]');
ylabel('Power reservoir [kW]');

%

figure;
hold on;
temp = nan(1,N_hours);
for ii = 1:N_E
    for jj = 1:N_hours
        temp(jj) = Potential_Pow{ii}(1,hours(jj)*3600/10 + 1);
    end
    plot(hours,temp,'o');
    
end
hold off;
grid on;
xlim([0,24]);
xlabel('Time [h]');
ylabel('Power reservoir [kW]');

%% Plot the N_{ON} set for the N_E experiments for the hours of the events

figure;
hold on;
temp = nan(1,N_hours);
for ii = 1:N_E
    for jj = 1:N_hours
        temp(jj) = N_ON{ii}(1,hours(jj)*3600/10 + 1);
    end
    plot(hours,temp,'o');
    
end
hold off;
grid on;
xlim([0,24]);
ylim([0,230]);
xlabel('Time [h]');
ylabel('Number of loads with the switch ON');
legend('N_{ON}','Location','Best');

%% Do hourly plots for N_{ON} and Potential_Pow on the same figure forall 
%  experiments

temp_N_ON = cell2mat(N_ON);
temp_Potential_Pow = cell2mat(Potential_Pow);

for ii = 1:N_hours
    figure;
    
    % plot the load population N_{ON}
    subplot(2,1,1);
    bar(temp_N_ON(:,hours(ii)*3600/10 + 1));
    grid on;
    xlabel('Experiment number');
    ylabel('Number of loads with the switch ON');
    legend(['Hour = ',num2str(hours(ii))],'Location','bestoutside');
    
    % plot the power reservoir
    subplot(2,1,2);
    bar(temp_Potential_Pow(:,hours(ii)*3600/10 + 1));
    grid on;
    xlabel('Experiment number');
    ylabel('Power reservoir [kW]');
    ylim([3000,5000]);
    
end

%% Plot the water draws

AA = cell2mat(Water_Draws);

mAA = mean(AA);
stdAA = std(AA);
lowAA = mAA - stdAA;
upperAA = mAA + stdAA;
% convert to time steps
mAA = mAA.*10;
lowAA = lowAA.*10;
upperAA = upperAA.*10;

xx = linspace(0,24,8641);

figure;
hold on;
% plot shaded area for the standard deviation around temperatures
fill([xx,fliplr(xx)],[upperAA,fliplr(lowAA)],...
    [0.7,0.7,0.7],'linestyle','none');
plot(xx,mAA);
hold off;
grid on;
ylim([0,65]);
xlim([0,24]);
xlabel('Time [h]');
ylabel('Water draw [liters/step]');
legend('Standard deviation','Mean water draw','Location','Best');

%% Plot the temperatures of the 2nd EWH disk (the one with the heating el.)
%
% time vector
xx = linspace(0,24,8641).';
% corresponds to T_{1+idx_disk} because T_{1} is the water inlet
idx_disk = 2;

% mTemperature = cell2mat(mTemperature);
% stdTemperature = cell2mat(stdTemperature);

for ii = 1:N_E
    figure;
    hold on;
    
    matrixTemps = mTemperature{ii};
    
    Temp_upper = matrixTemps + stdTemperature{ii};
    Temp_lower = matrixTemps - stdTemperature{ii};
    
    % plot shaded area for the standard deviation around temperatures
    fill([xx,fliplr(xx)],[Temp_upper(:,idx_disk),fliplr(Temp_lower(:,idx_disk))],...
        [0.9,0.9,0.9],'linestyle','none');
    plot(xx,matrixTemps(:,idx_disk));
    hold off;
end