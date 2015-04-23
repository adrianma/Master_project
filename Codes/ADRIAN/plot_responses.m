%% ===================================================================
%  function plot_responses
%  by Adrian Martinez Gomez
%  October 2014
%
%  Purpose:
%   * plot the responses for all the given Events and \etas
%   * save them for future use and analysis
%  ===================================================================

function plot_responses(Method,NFolder,bSaving);
%% Read the files names for later use

% choose the folder for the method (PS:prob. switching,
% SPV:set-point variation)
if(strcmp(Method,'ProbSwitching'))
    tFolderMethod = 'PS';
elseif(strcmp(Method,'SetPointVariation'))
    tFolderMethod = 'SPV';
end
% get the subfolder for the simulation case
tFolder = strcat('/',num2str(NFolder));

% get the file names
folder_files = dir(['VAGGELIS/saved_data/',tFolderMethod,tFolder,...
    '/*.mat']);
load(folder_files(1).name);

vTime = Params.t_init:Params.t_sample:Params.t_sim;

%% Iterate through all the files and plot the response
for ii = 2:length(folder_files)
    
    % load data-set
    load(folder_files(ii).name);
    disp(folder_files(ii).name);
    
    % store the time of the event for the lower subplot
    t_event_st = Results.Event(1);
    t_event_end = Results.Event(end);
    
    % plot and zoom-plot of the response
    handle = figure;
    subplot(2,1,1);
    hold on;
    plot(vTime./3600,Results.Prec);
    plot(vTime./3600,Results_comparison.Prec,'r');
    hold off;
    grid on;
    xlabel('Time [h]'); ylabel('P_{agg}');
    title(['Event at hour ',num2str(t_event_st/3600),' until ',...
        num2str(t_event_end/3600),' and \eta = ',num2str(Results.etas(1))]);
    legend(Method,'Only internal control','Location','Best');
    
    % time index vector to plot only after the event in the second subplot
    vIdx = find(vTime > t_event_st);
    
    % define confidence intervall for convergence
    confidence = 0.05;
    Prec_comparison_u = Results_comparison.Prec(vIdx).*(1 + confidence);
    Prec_comparison_l = Results_comparison.Prec(vIdx).*(1 - confidence);
    
    subplot(2,1,2);
    hold on;
    % plot shaded area around Results_comparison.Prec
    fill([vTime(vIdx)./3600,fliplr(vTime(vIdx)./3600)],...
        [Prec_comparison_u,fliplr(Prec_comparison_l)],[0.9,0.9,0.9],...
        'linestyle','none');
    plot(vTime(vIdx)./3600,Results.Prec(vIdx));
    plot(vTime(vIdx)./3600,Results_comparison.Prec(vIdx),'r');
    hold off;
    grid on;
    xlabel('Time [h]'); ylabel('P_{agg} [kW]');
    title(['Event at hour ',num2str(t_event_st/3600),' until ',...
        num2str(t_event_end/3600),' and \eta = ',num2str(Results.etas(1))]);
    legend([num2str(confidence),'% confidence from red line'],Method,...
        'Only internal control','Location','Best');
    
    % save the figure as png and fig (in case some manipulation is needed)
    %
    if(bSaving)
        saveas(handle,['VAGGELIS/plots/',tFolderMethod,'/',...
            num2str(NFolder),'/sim_',tFolderMethod,'_',...
            num2str(Params.n_app),'EWH_eta',...
            strrep(num2str(Results.etas(1)), '.', '')],'png');
        saveas(handle,['VAGGELIS/plots/',tFolderMethod,'/',...
            num2str(NFolder),'/sim_',tFolderMethod,'_',...
            num2str(Params.n_app),'EWH_eta',...
            strrep(num2str(Results.etas(1)), '.', '')],'fig');
    end
end