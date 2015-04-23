%% ========================================================================
%  function Baseline_Nloads
%  by Adrian Martinez Gomez
%  February 2015
%
%  Modified: April 2015
%
%  Input: bCalc; binary, 1 if it has to be calculated, 0 otherwise.
%
%  Purpose:
%  Generate Figure 2.3 of the Thesis (dynamics of 1 EWH for 1 day)
%  ========================================================================
function Baseline_Nloads(bCalc);
if(bCalc)
    % get the folders with all the baseline experiments
    folder_dir{1} = dir('archived_data/baseline/1EWH/res_same/*.mat');
    folder_dir{2} = dir('archived_data/baseline/5EWH/res_same/*.mat');
    folder_dir{3} = dir('archived_data/baseline/10EWH/res_same/*.mat');
    folder_dir{4} = dir('archived_data/baseline/50EWH/res_same/*.mat');
    folder_dir{5} = dir('archived_data/baseline/100EWH/res_same/*.mat');
    folder_dir{6} = dir('archived_data/baseline/500EWH/res_same/*.mat');
    folder_dir{7} = dir('archived_data/baseline/1000EWH/res_same/*.mat');
    folder_dir{8} = dir('archived_data/baseline/5000EWH/res_same/*.mat');
    
    % number of experiments N_E
    N_E = 20;
    % number of power for loads 10^(N_loads-1)
    N_loads = 7;
    % all the population sizes
    N_populations = [1,5,10,50,100,500,1000,5000];
    
    mPow = cell(1,N_E);
    meanPow = cell(1,N_loads);
    stdPow = cell(1,N_loads);
    
    mPow_norm = cell(1,N_loads);
    meanPow_norm = cell(1,N_loads);
    stdPow_norm = cell(1,N_loads);
    
    mN_ON = cell(N_loads,N_E);
    for ii = 1:N_loads
        if(ii==1)
            path_file = 'archived_data/baseline/1EWH/res_same/';
        elseif(ii==2)
            path_file = 'archived_data/baseline/5EWH/res_same/';
        elseif(ii==3)
            path_file = 'archived_data/baseline/10EWH/res_same/';
        elseif(ii==4)
            path_file = 'archived_data/baseline/50EWH/res_same/';
        elseif(ii==5)
            path_file = 'archived_data/baseline/100EWH/res_same/';
        elseif(ii==6)
            path_file = 'archived_data/baseline/500EWH/res_same/';
        elseif(ii==7)
            path_file = 'archived_data/baseline/1000EWH/res_same/';
        elseif(ii==8)
            path_file = 'archived_data/baseline/5000EWH/res_same/';
        end
        
        if(mod(ii,2))
            disp(['... Loading files with 10^',num2str(ii-1),' EWHs']);
        else
            disp(['... Loading files with 5*10^',num2str(ii-1),' EWHs']);
        end
        
        for jj = 1:length(folder_dir{ii})
            load(strcat(path_file,folder_dir{ii}(jj).name));
            
            mPow{jj} = Results_comparison.Prec;
            mPow_norm{jj} = mPow{jj}./sum(Params.P1_el);
            
            mN_ON{ii,jj} = sum(squeeze(Results_comparison.urec(1,:,:)).');
        end
        
        meanPow{ii} = mean(cell2mat(mPow.'));
        stdPow{ii} = std(cell2mat(mPow.'));
        
        meanPow_norm{ii} = mean(cell2mat(mPow_norm.'));
        stdPow_norm{ii} = std(cell2mat(mPow_norm.'));
    end
    
    save('ADRIAN/Baseline_Nloads/1to1000.mat');
else
    load('ADRIAN/Baseline_Nloads/1to1000.mat');
end

%% Plots
%
indeces = [1,3,5,7];
xx = linspace(0,24,8641);

% 1)
figure;
count = 1;
for ii = indeces
    mm = meanPow_norm{ii};
    ss = stdPow_norm{ii};
    low = mm - ss;
    upp = mm + ss;
    
    subplot(2,2,count);
    hold on;
    % plot shaded area for the standard deviation around temperatures
    fill([xx,fliplr(xx)],[upp,fliplr(low)],...
        [0.7,0.7,0.7],'linestyle','none');
    plot(xx,mm);
    hold off;
    grid on;
    xlim([0,24]);
    title(['N_{app} = ',num2str(N_populations(ii))],'FontSize',12);
    legend('Standard deviation','Mean','Location','Best');
    
    count = count + 1;
end
subplot(2,2,1); 
ylabel('Normalized aggregate power','FontSize',12); 
% adjust the limits of the y-axis to make comparison more understandable
subplot(2,2,2); ylim([-0.5,1]);
subplot(2,2,3); 
ylabel('Normalized aggregate power','FontSize',12);
% adjust the limits of the y-axis to make comparison more understandable
subplot(2,2,4); ylim([-0.1,0.3]);

% 2)
figure
hold on;
cmap = hsv(N_loads - 4);
for ii = 5:N_loads
    plot(xx,stdPow_norm{ii},'Color',cmap(ii-4,:));
end
hold off;
grid on;
ylabel('Standard deviation \sigma','FontSize',12);
xlabel('Time [h]','FontSize',12);
xlim([0,24]);
legend(['N_{app}=',num2str(10)],['N_{app}=',num2str(100)],...
    ['N_{app}=',num2str(1000)],'Location','Best');

%% BIN
if(0)
    clear;
    n_app = 5000;
    Params = build_population(n_app);
    % precomputes the preliminary system model
    PrelimModel = precompute_system_model(Params);
    for ii = 1:20
        WaterDrawScenarioReal = build_normalized_draw_scenario(Params);
        % water draws required
        Params.bWaterDraw = 1;
        Results_comparison = simulate_population(Params,PrelimModel,...
            WaterDrawScenarioReal,WaterDrawScenarioReal,'NoControl');
        % save results in file for future analysis
        save(['archived_data/baseline/5000EWH/res_same/sim_',...
            num2str(n_app),'EWHs_runNumber',num2str(ii),'.mat'],...
            'Params','Results_comparison','n_app','PrelimModel','WaterDrawScenarioReal');
    end
    
end

end