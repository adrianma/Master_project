%% ========================================================================
%  function errorsAppendix_A
%  by Adrian Martinez Gomez
%  March 2015
%
%  Inputs:
%   *opt_signs: 1, when both external inputs are positive (folder PosPos)
%               2, when different signs (folder PosNeg)
%               3, when both external inputs are negative (folder NegNeg)
%
%   *opt_error: is a binary number (0 for relative error, 1 for RMSE)
%       almost always select the RMSE!
%       The other one is prone to errors and should be debugged (if time
%       allows)
%
%   *bAddHom: number, giving the file to load. 
%    Choose "1" for Additivity and "0" for Homogeneity and "2" for Maryams
%    examples.
%
%  Purpose:
%  Evaluate the examples calculated for the additivity property for the
%  system.
%       do the call with: errorsAppendix_A(1,1,1);
%  ========================================================================
function errorsAppendix_A(opt_signs,opt_error,bAddHom);
myBool7 = 1;
%% 1)
%
if(opt_signs == 1)
    topt_signs = 'PosPos';
elseif(opt_signs == 2)
    topt_signs = 'PosNeg';
elseif(opt_signs == 3)
    topt_signs = 'NegNeg';
end

% load the stored file
if(bAddHom == 1)
    load(['archived_data/Linearity_SPV/',...
        num2str(topt_signs),'/my_sim_additivity_4_12_20.mat']);
    if(myBool7)
        load(['archived_data/Linearity_SPV/',...
            num2str(topt_signs),'/my_sim_additivity_7.mat']);
        % gather the data from both experiments
        u_1_inputs = [u_1_inputs;u_1_inputs_7];
        u_2_inputs = [u_2_inputs;u_2_inputs_7];
        u_3_inputs = [u_3_inputs;u_3_inputs_7];
        myHours = [myHours,myHours_7];
        %
        uu = [uu;uu_7];
        yy = [yy;yy_7];
    end
elseif(bAddHom == 0)
    load(['archived_data/Linearity_SPV/',...
        num2str(topt_signs),'/my_sim_homogeneity.mat']);
elseif(bAddHom == 2)
    load(['archived_data/Linearity_SPV/',...
        num2str(topt_signs),'/my_sim_zero_Maryam.mat']);
end

%% 2) Gather the errors
%
% number of samples
N_samples = size(yy,2);
% number of hours
N_hours = length(myHours);

% for this get the residuals for comparison with the case u_1+u_2=0
if(bAddHom == 2)
    load('archived_data/baseline/1000EWH/same_baseline_data/power_baseline_same.mat');
    idx = 1;
    y_diffs = P_no_controller_norm(idx,:) - y_baseline_norm;
    
    t_sample = 10;
    N_zeros = 20;
    myDiffs = cell(1,N_hours);
    for ii = 1:N_hours
        % resolve the indices
        low_idx = (myHours(ii)*3600)/t_sample + 1;
        high_idx = ((myHours(ii)*3600) + 1*3600)/t_sample + 1;
        
        myDiffs{ii} = y_diffs(low_idx-N_zeros:high_idx).';
    end
end

% pre-allocation
errRMSE = cell(N_hours,N_samples);
errRel = cell(N_hours,N_samples);
sum_abs = cell(N_hours,N_samples);
% loop through all samples and hours
for ii = 1:N_samples
    for jj = 1:N_hours
        
        % Do Additivity
        if(bAddHom == 1)
            % calculate the RMSE
            errRMSE{jj,ii} = norm(abs(yy{jj,ii}(:,3) - (yy{jj,ii}(:,1) + ...
                yy{jj,ii}(:,2))),2)./sqrt(length(yy{jj,ii}(:,3)));
            
            % calculate the relative error
            errRel{jj,ii} = abs(yy{jj,ii}(:,3) - ...
                (yy{jj,ii}(:,1)+yy{jj,ii}(:,2)))./abs(yy{jj,ii}(:,3));
            % Do Homogeneity
        elseif(bAddHom == 0)
            % calculate the RMSE
            errRMSE{jj,ii} = norm(abs(yy{jj,ii}(:,2) - (alphas(ii).*...
                yy{jj,ii}(:,1))),2)./sqrt(length(yy{jj,ii}(:,3)));
            
            % calculate the relative error
            errRel{jj,ii} = abs(yy{jj,ii}(:,3) - ...
                (alphas(ii).*yy{jj,ii}(:,1)))./abs(yy{jj,ii}(:,3));
        elseif(bAddHom == 2)
            % calculate the RMSE
            errRMSE{jj,ii} = norm(abs(myDiffs{jj} - (yy{jj,ii}(:,1) + ...
                yy{jj,ii}(:,2))),2)./sqrt(length(myDiffs{jj}));
            
             % calculate the relative error
            errRel{jj,ii} = abs(myDiffs{jj}  - ...
                (yy{jj,ii}(:,1)+yy{jj,ii}(:,2)))./abs(myDiffs{jj});
        end
        
        % take 2nd norm of the vector
        errRel{jj,ii} = norm(errRel{jj,ii},2);
        
        % the absolute value of the sum of inputs (magnitude)
        sum_abs{jj,ii} = abs(uu{jj,ii}(23,3));
        
    end
end
% convert back to matrices
%
% select error and convert to [%]
if(opt_error)
    err = 100.*cell2mat(errRMSE);
    text_ylabel = 'RMSE [%]';
    
    disp('  => Use RMSE as the error metric!');
else
    err = 100.*cell2mat(errRel);
    text_ylabel = 'Relative error [%]';
    
    disp('  => Use Relative error as the error metric!');
end
sum_abs = cell2mat(sum_abs);

%% 3) Plot the errors (RMSE or relative error)
%
% (3.1)
figure;
hold on;
stem(err(1,:),'b');
stem(err(2,:),'r');
stem(err(3,:),'g');
stem(err(4,:),'k');
grid on;
xlabel('Experiment number','FontSize',12);
ylabel(text_ylabel,'FontSize',12);
legend(['Hour = ',num2str(myHours(1))],['Hour = ',num2str(myHours(2))],...
    ['Hour = ',num2str(myHours(3))],['Hour = ',num2str(myHours(4))],...
    'Location','Best');

% (3.2)
figure;
hist(err.');
grid on;
ylabel('Experiment number','FontSize',12);
xlabel(text_ylabel,'FontSize',12);
legend(['Hour = ',num2str(myHours(1))],['Hour = ',num2str(myHours(2))],...
    ['Hour = ',num2str(myHours(3))],['Hour = ',num2str(myHours(4))],...
    'Location','Best');

% Predefine color
myColors_1 = repmat([1,0,0],length(u_1_inputs),1);
myColors_2 = repmat([0,1,0],length(u_2_inputs),1);
myColors_3 = repmat([0,0,1],length(u_3_inputs),1);
myColors = [myColors_1;myColors_2;myColors_3];

% (3.2)
%
% only do this one if hour 7 is not present
if(~myBool7)
    figure;
    my_h = scatter3(repmat(u_1_inputs,1,3),repmat(u_2_inputs,1,3),...
        err(:),20,myColors);
    my_h.MarkerFaceColor = [0 .75 .75];
    grid on;
    xlabel('u_1','FontSize',12);
    ylabel('u_2','FontSize',12);
    zlabel(text_ylabel,'FontSize',12);
end

% (3.3)
if(myBool7 && ~(bAddHom==2))
    figure;
    scatter3([repmat(u_1_inputs(1,:),1,3),u_1_inputs(2,:)],...
        [repmat(u_2_inputs(1,:),1,3),u_2_inputs(2,:)],...
        err(:));
    grid on;
    xlabel('u_1','FontSize',12);
    ylabel('u_2','FontSize',12);
    zlabel(text_ylabel,'FontSize',12);
elseif(~(bAddHom==2))
    figure;
    scatter3(repmat(u_1_inputs,1,3),repmat(u_2_inputs,1,3),...
        err(:),20,myColors);
    grid on;
    xlabel('u_1','FontSize',12);
    ylabel('u_2','FontSize',12);
    zlabel(text_ylabel,'FontSize',12);
end

% (3.4)
figure;
subplot(2,1,1);
scatter(sum_abs(:),err(:));
grid on;
ylabel(text_ylabel,'FontSize',12);

subplot(2,1,2);
hold on;
scatter(sum_abs(1,:),err(1,:),'b');
scatter(sum_abs(2,:),err(2,:),'r');
scatter(sum_abs(3,:),err(3,:),'g');
scatter(sum_abs(4,:),err(4,:),'k');
hold off;
grid on;
xlabel('|u_3|','FontSize',12);
ylabel(text_ylabel,'FontSize',12);
legend(['Hour = ',num2str(myHours(1))],['Hour = ',num2str(myHours(2))],...
    ['Hour = ',num2str(myHours(3))],['Hour = ',num2str(myHours(4))],...
    'Location','Best');

% (3.5)
%
% Result is that error is smaller if one input is significantly bigger than
% the other (u_1 >> u_2). Didn't expect that!
figure;
hold on;
scatter(sum_abs(1,1:100),err(1,1:100),'b');
scatter(sum_abs(1,101:end),err(1,101:end),'b*');
scatter(sum_abs(2,1:100),err(2,1:100),'r');
scatter(sum_abs(2,101:end),err(2,101:end),'r*');
scatter(sum_abs(3,1:100),err(3,1:100),'g');
scatter(sum_abs(3,101:end),err(3,101:end),'g*');
scatter(sum_abs(4,1:100),err(4,1:100),'k');
scatter(sum_abs(4,101:end),err(4,101:end),'k*');
hold off;
grid on;
xlabel('|u_3|','FontSize',12);
ylabel(text_ylabel,'FontSize',12);
legend(['Hour = ',num2str(myHours(1))],['Hour = ',num2str(myHours(1)),'; u_1>>u_2'],...
    ['Hour = ',num2str(myHours(2))],['Hour = ',num2str(myHours(2)),'; u_1>>u_2'],...
    ['Hour = ',num2str(myHours(3))],['Hour = ',num2str(myHours(3)),'; u_1>>u_2'],...
    ['Hour = ',num2str(myHours(4))],['Hour = ',num2str(myHours(4)),'; u_1>>u_2'],...
    'Location','Best');


fprintf(' => End of script\n');

end