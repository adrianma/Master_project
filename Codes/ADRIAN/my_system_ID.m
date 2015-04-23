%% ===================================================================
%  function my_system_ID
%  by Adrian Martinez Gomez
%  October 2014
%
%  Inputs:
%   * Method: string with the method selected (prob. switching or set-point
%   variation.
%   * bDifference: boolean. Select 1 if you want to analyze the variable
%   part of the aggregate power.
%   * NAnalysis: number of the folder where the files are at
%       1 file contains the experiments and the other the baseline case
%       (without external controller).
%
%  Purpose:
%   * compute the ARX/ARMAX/Nonlinear models and compare their performance.
%  ===================================================================

function Output = my_system_ID(Method,bInput,sHour,bWaterDraw,NAnalysis);

% define the sampling time (according to all experiments, it is defined as
% 10 seconds)
Params.t_sample = 10;

% choose the folder for the method (PS:prob. switching,
% SPV:set-point variation)
if(strcmp(Method,'ProbSwitching'))
    tFolderMethod = 'PS';
elseif(strcmp(Method,'SetPointVariation'))
    tFolderMethod = 'SPV';
end
tFolder = strcat('/',num2str(NAnalysis));

% get the file names
folder_files = dir(['VAGGELIS/saved_data/',tFolderMethod,tFolder,...
    '/*.mat']);
% load the baseline case and the responses for the given times
load(folder_files(1).name);
load(folder_files(2).name);

%% (1) Define time vector and confidence interval
Output.vTime = Params.t_init:Params.t_sample:Params.t_sim;

%% (2) Iterate through the experiments to construct the input and the output
%     signals.
Output.u = cell(1,N_E);
Output.y = cell(1,N_E);
Output.y_diff = cell(1,N_E);
Output.period = cell(1,N_E);
for ii = 1:N_E
    % resolve the indices
    low_idx = event_hour/Params.t_sample + 1;
    high_idx = (event_hour + 1*3600)/Params.t_sample + 1;
    low_event_idx = Events{ii}(1)/Params.t_sample + 1 - low_idx;
    high_event_idx = Events{ii}(end)/Params.t_sample + 1 - low_idx;
    
    % output y
    Output.y{ii} = Results{ii}.Prec;
    Output.y_diff{ii} = Results{ii}.Prec - ...
        Results_comparison.Prec(low_idx:high_idx);
    
    % input u
    if(bInput)
        Output.u{ii} = Results{ii}.etas;
    else
        Output.u{ii} = [Results{ii}.etas(1).*ones(1,high_event_idx - ...
            low_event_idx + 1),zeros(1,(high_idx-low_idx)- high_event_idx)];
        % are zeros needed before the Event takes place?
        %   (add them if necessary!)
        if(low_event_idx-1 ~= 0)
            Output.u{ii} = [zeros(1,low_event_idx-1),Output.u{ii}];
        end
        
        % catch error (to ensure that u and y are of the same size)
        if(length(Output.y{ii}) ~= length(Output.u{ii}))
            fprintf('Error in sizes of u and y\n');
            Output.u{ii} = [Output.u{ii},...
                zeros(1,length(Output.y{ii})-length(Output.y{ii}))];
        end
    end
    
    % period Ts
    Output.period{ii} = Params.t_sample;
    
    % Every entry of the 1xN_E cell-array has to be of size N_samplesx1
    % (not 1xN_samples)
    Output.y{ii} = Output.y{ii}.';
    Output.y_diff{ii} = Output.y_diff{ii}.';
    
    % also, remember to normalyze the power demand
    if(1)
        Output.y{ii} = Output.y{ii}./sum(Params.P1_el);
        Output.y_diff{ii} = Output.y_diff{ii}./sum(Params.P1_el);
    end
    
    Output.u{ii} = Output.u{ii}.';
    
    % exogenous input to add the water draw probability into the model
    %
    % 1st) take the empirical probability
    % 2nd) take the SUM of the water draws at each time instance
    if(bWaterDraw(2))
        temp = extend_water_draw(Params,sHour);
        u_exogenous = cell(1,N_E);
        for jj = 1:N_E
            u_exogenous{jj} = temp;
        end
    else
        u_exogenous = cell(1,N_E);
        for jj = 1:N_E
            temp = squeeze(Results{jj}.mdotrec);
            temp = sum(temp.').';
            u_exogenous{jj} = temp;
        end
    end
end

%% (3) Get the data ready (& preprocessing)

%---- 3.0) create multiexperiment data AND validation data ----%
if(bWaterDraw(1))
    fprintf('\n Water draw considered into the model!\n');
    
    % reassemble the input, considering the exogenous input (water draw)
    for ii = 1:N_E
        Output.u{ii} = [Output.u{ii},u_exogenous{ii}];
    end
else
    fprintf('\n No water draw considered into the model!\n');
end
% create the iddata object for the experimental- and validation data
Output.data = iddata(Output.y(1:floor(N_E/2)),Output.u(1:floor(N_E/2)),...
    Output.period(1:floor(N_E/2)));
Output.validation = iddata(Output.y((floor(N_E/2)+1):end),...
    Output.u((floor(N_E/2)+1):end),Output.period((floor(N_E/2)+1):end));
Output.data_diff = iddata(Output.y_diff(1:floor(N_E/2)),...
    Output.u(1:floor(N_E/2)),Output.period(1:floor(N_E/2)));
Output.validation_diff = iddata(Output.y_diff((floor(N_E/2)+1):end),...
    Output.u((floor(N_E/2)+1):end),Output.period((floor(N_E/2)+1):end));

% remove the steady-state (=internal controller-only) part of the curve by
% using the difference
%
Output.data = Output.data_diff;
Output.validation = Output.validation_diff;

%---- 3.1) get advice from command ----%
if(0)
    advice(Output.data);
end

%---- 3.2) handling missing data ----%
dat1 = misdata(Output.data);
% Check how the missing data was estimated on a time plot
if(0)
    figure; plot(Output.data,dat1); grid on;
end

Output.data = dat1;

%---- 3.3) does the system have a time delay? ----%
temp_delays = delayest(Output.data);
disp(['The system has n_k=',num2str(temp_delays),' delays']);

%---- 3.4) detrend the data ----%
if(0)
    T = getTrend(Output.data);
    
    T.InputOffset = I_value;
    T.OutputOffset = O_value;
end
Output.data = detrend(Output.data);

%% (4) Get a model estimate

%---- 4.1) a linear one with ARX/ARMAX structure ----%
if(0)
    % generate model-order combinations
    NN = struc(1:10,1:10,temp_delays);
    % estimate an ARX model for each model order.
    %   1st argument is the experiment data-set
    %   2nd argument is the validation data-set
    V_iv = ivstruc(Output.data,Output.validation,NN);
    V_arx = arxstruc(Output.data,Output.validation,NN);
    % select a model order.
    order_iv = selstruc(V_iv,0);
    order_arx = selstruc(V_arx,0);
    
    % for ARX
    n_a = max(order_iv(1),order_arx(1));
    n_b = max(order_iv(2),order_arx(2));
    n_k = max(order_iv(3),order_arx(3));
    % for ARMAX
    na = n_a;
    nb = n_b;
    nk = n_k;
    
    % find the best coefficient nc
    % (given that na and nb are estimated from the upper code snippet)
    ncs = 1:1:10;
    xxx_1 = armax(Output.data,[na,nb,ncs(1),nk]);
    xxx_2 = armax(Output.data,[na,nb,ncs(2),nk]);
    xxx_3 = armax(Output.data,[na,nb,ncs(3),nk]);
    xxx_4 = armax(Output.data,[na,nb,ncs(4),nk]);
    xxx_5 = armax(Output.data,[na,nb,ncs(5),nk]);
    xxx_6 = armax(Output.data,[na,nb,ncs(6),nk]);
    xxx_7 = armax(Output.data,[na,nb,ncs(3),nk]);
    xxx_8 = armax(Output.data,[na,nb,ncs(4),nk]);
    xxx_9 = armax(Output.data,[na,nb,ncs(5),nk]);
    xxx_10 = armax(Output.data,[na,nb,ncs(6),nk]);
    
    [~,fit,~] = compare(Output.validation,xxx_1,xxx_2,xxx_3,xxx_4,...
        xxx_5,xxx_6,xxx_7,xxx_8,xxx_9,xxx_10);
    
    % which nc from the vector ncs provides the best fit?
    vfit = sum(cell2mat(fit.'));
    [max_value,max_idx] = max(vfit);
    
    disp(['   The best order for n_c is ',num2str(ncs(max_idx))]);
    disp(['   And its average fit for validation is ',...
        num2str(max_value/(N_E/2))]);
    
    % set the best possible value for nc
    nc = ncs(max_idx);
end

%% find best coefficients one and for all! (brute force 3xloop FTW!)
if(1)
    if(bWaterDraw(1))
        nas = 1:1:10;
        nbs_1 = 1:1:10;
        nbs_2 = 1:1:10;
        ncs = 1:1:10;
        % there is no delay to be seen from inspection, therefore use nk=0
        nk = [0,0];
        
        max_value = 0; na = 0; nb_1 = 0; nb_2 = 0; nc = 0;
        for ii = 1:length(nas)
            for jj = 1:length(nbs_1)
                for kk = 1:length(nbs_2)
                    
                    xxx = cell(1,length(ncs));
                    for ll = 1:length(ncs)
                        xxx{ll} = armax(Output.data,[nas(ii),...
                            [nbs_1(jj),nbs_2(kk)],ncs(ll),nk]);
                    end
                    
                    [~,fit,~] = compare(Output.validation,xxx{1},xxx{2},...
                        xxx{3},xxx{4},xxx{5},xxx{6},xxx{7},xxx{8},...
                        xxx{9},xxx{10});
                    
                    % which nc from the vector ncs provides the best fit?
                    vfit = sum(cell2mat(fit.'));
                    [temp_max_value,max_idx] = max(vfit);
                    
                    if(temp_max_value > max_value)
                        max_value = temp_max_value;
                        na = nas(ii);
                        nb_1 = nbs_1(jj);
                        nb_2 = nbs_2(jj);
                        nc = ncs(max_idx);
                    end
                end
            end
        end
        
        disp(['%%%%% Hour number ',num2str(sHour),' %%%%%']);
        disp(['  na =',num2str(na)]);
        disp(['  nb_1 =',num2str(nb_1)]);
        disp(['  nb_2 =',num2str(nb_2)]);
        disp(['  nc =',num2str(nc)]);
        disp(['  Fit =',num2str(max_value)]);        
    else
        nas = 1:1:10;
        nbs = 1:1:10;
        ncs = 1:1:10;
        % there is no delay to be seen from inspection, therefore use nk=0
        nk = 0;
        
        max_value = 0; na = 0; nb = 0; nc = 0;
        for ii = 1:length(nas)
            for jj = 1:length(nbs)
                xxx = cell(1,length(ncs));
                for ll = 1:length(ncs)
                    xxx{ll} = armax(Output.data,[nas(ii),nbs(jj),...
                        ncs(ll),nk]);
                end
                
                [~,fit,~] = compare(Output.validation,xxx{1},xxx{2},xxx{3},...
                    xxx{4},xxx{5},xxx{6},xxx{7},xxx{8},xxx{9},xxx{10});
                
                % which nc from the vector ncs provides the best fit?
                vfit = sum(cell2mat(fit.'));
                [temp_max_value,max_idx] = max(vfit);
                
                if(temp_max_value > max_value)
                    max_value = temp_max_value;
                    na = nas(ii);
                    nb = nbs(jj);
                    nc = ncs(max_idx);
                end
            end
        end
        
        disp(['%%%%% Hour number ',num2str(sHour),' %%%%%']);
        disp(['  na =',num2str(na)]);
        disp(['  nb =',num2str(nb)]);
        disp(['  nc =',num2str(nc)]);
        disp(['  Fit =',num2str(max_value)]);
    end
else
    % values that were obtained running the above loops, and gave the
    % best/highest fit (-> for Prob. Switching)
    if(event_hour == 18*3600 && strcmp(Method,'ProbSwitching'))
        na = 10;
        nb = 1;
        nc = 6;
    elseif(event_hour == 12*3600 && strcmp(Method,'ProbSwitching'))
        na = 8;
        nb = 1;
        nc = 10;
    end
    nk = 0;
    
    [yy,fit,~] = compare(Output.validation,armax(Output.data,[na,nb,nc,nk]));
    Output.armax_best_fit = sum(cell2mat(fit.'))/(N_E/2);
    
    % try calculating the MSE for the estimation
    errors = cell(N_E/2,1);
    MSEs = nan(N_E/2,1);
    for ii = 1:N_E/2
        errors{ii} = Output.y_diff{ii} - yy{ii}.OutputData;
        MSEs(ii) = mean(errors{ii}.^2);
    end
    % plot the MSE
    figure;
    bar(MSEs);
    grid on;
    xlabel('experiment');
    ylabel('MSE');
    
    disp(['   The best orders are (n_a,n_b,n_c) = (',num2str(na),',',...
        num2str(nb),',',num2str(nc),')']);
    disp(['   And its fit for the validation data is ',...
        num2str(Output.armax_best_fit),'%']);
    disp(['   Also, the MSEs are as follows (remember they are normed!) ',...
        num2str(mean(MSEs)),'%']);
end

%% (4) Model selection (linear version)
if(bWaterDraw(1))
    % correct call for the MISO system (with exogenous input)
    Output.armax_model = armax(Output.data,[3,[3,3],3,[0 0]]);
else
    % when NO exogenous input (water draw profile) regarded
    %
    Output.arx_model = arx(Output.data,[n_a,n_b,n_k]);
    Output.armax_model = armax(Output.data,[na,nb,nc,nk]);
    Output.iv_model = iv4(Output.data,[n_a,n_b,n_k]);
    % Output.bj_model = bj(Output.data,[]);
end

%---- 4.2) a nonlinear ARX model estimate ----%

advice_orders = [4,4,1];
disp('Follow the advice function to get the orders of the system![4,4,1]');
disp('  Also, use treepartition nonlinearity structure');
nonlinearity = treepartition('num',100);

Output.nl_arx = nlarx(Output.data,advice_orders,nonlinearity);

%% (5) Model validation
if(0)
    %---- 5.1) compute the residuals ----%
    res_arx = resid(Output.arx_model,Output.validation,'corr');
    Output.res{1} = res_arx;
    res_armax = resid(Output.armax_model,Output.validation,'corr');
    Output.res{2} = res_armax;
    res_nl_arx = resid(Output.nl_arx,Output.validation,'corr');
    Output.res{3} = res_nl_arx;
    
    % plot the residuals for the 3 models
    if(0)
        figure;
        hold on;
        plot(Output.res{1});
        plot(Output.res{2},'r');
        plot(Output.res{3},'k');
        hold off;
        grid on;
    end
    
    %---- 5.2) display model info (including estimated uncertainty) ----%
    if(0)
        present(Output.arx_model);
        present(Output.armax_model);
        present(Output.nl_arx);
    end
    
    %---- 5.3) comparison between the models for a new data set (val) ----%
    [Output.comp_y,Output.comp_fit,Output.comp_x0] = compare(...
        Output.validation,Output.arx_sys,Output.armax_sys,Output.nl_arx);
    % plot the validation(s)
    compare(Output.validation,Output.arx_model,Output.armax_model,...
        Output.nl_arx);
    
    %---- 5.4) compute and plot model prediction errors ----%
    Output.err_arx = pe(Output.arx_model,Output.data);
    Output.err_armax = pe(Output.armax_model,Output.data);
    Output.err_nl_arx = pe(Output.nl_arx,Output.data);
    
    %---- 5.5) plot the noise spectrum ----%
    v = y_measured - y_simulated;
    Output.noise_spectrum = spa(v);
    % plot figure
    if(0)
        spectrum(Output.noise_spectrum,m);
    else
        figure;
        hold on;
        spectrumplot(Output.arx_sys,'b*--');
        spectrumplot(Output.armax_sys,'r*--');
        spectrumplot(Output.nl_arx,'k*--');
        hold off;
        legend('ARX','ARMAX','Non-linear ARX','Location','Best');
    end
    
    %% (6) Extras
    
    %---- 6.1) calculate the ETFE ----%
    g = etfe(Output.data);
    g_diff = etfe(Output.data_diff);
    
    figure; bode(g);
    figure; bode(g_diff);
    
    %---- 6.1.5) calculate ETFE the hard way ----%
    %
    % from Roy Smith's slides
    if(0)
        U = fft(u); % calculate N point FFTs
        Y = fft(y);
        N = length(Y);
        omega = (2*pi/N)*[0:N-1]'; % frequency grid
        idx = find(omega > 0 & omega < pi); % positive frequencies
        loglog(omega(idx),abs(U(idx)));
        loglog(omega(idx),abs(Y(idx)));
        Gest = Y./U; % ETFE estimate
        Gfresp = squeeze(freqresp(G,omega)); % "true" system response
        loglog(omega(idx),abs(Gest(idx)));
        semilogx(omega(idx),angle(Gest(idx)));
        Err = Gest - Gfresp; % calculate error
        loglog(omega(idx),abs(Err(idx)));
    end
    
    %---- 6.2) estimate the frequency response ----%
    Ge = spa(ze);
    Ge_diff = spa(ze);
    
    figure; bode(Ge);
    figure; bode(Ge_diff);
end

%% BIN
% %% Gather all the ARX model systems
% Output.eta = nan(1,length(folder_files)-1);
% Output.u = nan(length(folder_files)-1,length(Output.vTime));
% Output.y = nan(length(folder_files)-1,length(Output.vTime));
% Output.system_arx = cell(1,length(folder_files)-1);
% Output.system_tfest = cell(1,length(folder_files)-1);
% for ii = 2:length(folder_files)
%
%     % load data-set
%     load(folder_files(ii).name);
%     disp(folder_files(ii).name);
%
%     % store the eta (for later sorting)
%     Output.eta(ii-1) = Results.etas(1);
%
%     % gather u and y
%     sIdx = 1 + (sum(sum(Results.urec(2,:,:))) ~= 0);
% %     Output.u(ii-1,:) = sum(squeeze(Results.urec(sIdx,:,:)).') - ...
% %         sum(squeeze(Results_comparison.urec(sIdx,:,:)).');
%     Output.u(ii-1,:) = define_ref_signal(Params,Results);
%     Output.y(ii-1,:) = Results.Prec - Results_comparison.Prec;
%
%     temp = iddata(Output.y(ii-1,:).',Output.u(ii-1,:).',Output.period);
%
%     % 1) use the ARX model
%     %
%     % TODO: como escoger n_a y n_b basado en los experimentos??
%     %   (mirar SysID script para ideas)
%     n_a = 5;
%     n_b = 5;
%     n_k = 0;
%     Output.system_arx{ii-1} = arx(temp,[n_a n_b n_k]);
%
%     % 2) use Matlab's built-in 'tfest'
%     n_p = 5;
%     n_z = 5;
%     Output.system_tfest{ii-1} = tfest(temp,n_p,n_z,'Ts',Params.t_sample);
%
% end
%
% % for convenience, sort the results
% [Output.eta,vIdx] = sort(Output.eta);
% Output.system_arx(vIdx);
% Output.system_tfest(vIdx);
%
%
% %% plot
% if(0)
%     figure; hold on;
%     cmap = hsv(length(folder_files)-1);
%
%     for ii = 1:length(folder_files)-1
%         bode(Output.system_tfest{ii});
%     end
% end
%
%
%
% % ARX: compute the NoiseVariance term needed
% temp_arx = polyest(Output.data,[n_a,n_b,1,1,1,n_k]);
% Output.arx_sys = idpoly(Output.arx_model.A,Output.arx_model.B,...
%     Output.arx_model.C,Output.arx_model.D,Output.arx_model.F,...
%     temp_arx.NoiseVariance,Output.period{1});
%
% % ARMAX: compute the NoiseVariance term needed
% temp_armax = polyest(Output.data,[na,nb,nc,1,1,nk]);
% Output.armax_sys = idpoly(Output.armax_model.A,Output.armax_model.B,...
%     Output.armax_model.C,Output.armax_model.D,Output.armax_model.F,...
%     temp_armax.NoiseVariance,Output.period{1});
%
% % IV: compute the NoiseVariance term needed
% temp_iv = polyest(Output.data,[na,nb,1,1,1,nk]);
% Output.iv_sys = idpoly(Output.armax_model.A,Output.armax_model.B,...
%     Output.armax_model.C,Output.armax_model.D,Output.armax_model.F,...
%     temp_iv.NoiseVariance,Output.period{1});
end

function u_exogenous = extend_water_draw(Params,sHour)

% define and plot the water draw probability
prob_vector = [0.0129 0.0091 0.0091 0.0098 0.0144 0.0263 0.0491 0.054 0.054 0.0538 0.0515 0.0492 0.047 0.043 0.042 0.043 0.0492 0.061 0.0695 0.0685 0.0623 0.0528 0.0435 0.025]';

u_exogenous = [prob_vector(sHour).*ones(1,1*3600/Params.t_sample),0].';

end