%% ========================================================================
%  function Assess_hourly
%  by Adrian Martinez Gomez
%  February 2015
%
%  Inputs:
%   Exp_Data : experimental data
%   Val_Data : validation data
%   rOptions : struct with different options on how to perform the analysis
%       Fields: {tMethod,tFocus,tSearchMethod,bIN,myOrder,myOrderNoise,
%       bFeedThrough}
%
%       EXAMPLE FOR rOptions:
%       rOptions.tMethod = 'ARMAX';rOptions.myOrder = 4; rOptions.myOrderNoise = 1;rOptions.tFocus = 'simulation';
%
%  Purpose:
%  This function does generate ARMAX or State-Space models, and performs
%  the validation to the data.
%  Note that it only works for hours 4,7,12,16,18,20,22; you have to adapt
%  it if you want to test others.
%
%  ========================================================================
function [Models,meanFit,stdFit] = Assess_hourly(Exp_Data,Val_Data,...
    rOptions);
%% 0) Initializations

disp('  ==> Start by splitting the data into the hourly intervals...');

% 0.1) split the data into hours (experimental and validation data)
[exp_4,exp_7,exp_12,exp_16,exp_18,exp_20,exp_22] = ...
    Split_hourly_data(Exp_Data);
[val_4,val_7,val_12,val_16,val_18,val_20,val_22] = ...
    Split_hourly_data(Val_Data);

% 0.2) split the hourly data into positive and negative data 
[exp_4_pos,exp_4_neg] = PosNeg_system_ID(exp_4);
[exp_7_pos,exp_7_neg] = PosNeg_system_ID(exp_7);
[exp_12_pos,exp_12_neg] = PosNeg_system_ID(exp_12);
[exp_16_pos,exp_16_neg] = PosNeg_system_ID(exp_16);
[exp_18_pos,exp_18_neg] = PosNeg_system_ID(exp_18);
[exp_20_pos,exp_20_neg] = PosNeg_system_ID(exp_20);
[exp_22_pos,exp_22_neg] = PosNeg_system_ID(exp_22);

[val_4_pos,val_4_neg] = PosNeg_system_ID(val_4);
[val_7_pos,val_7_neg] = PosNeg_system_ID(val_7);
[val_12_pos,val_12_neg] = PosNeg_system_ID(val_12);
[val_16_pos,val_16_neg] = PosNeg_system_ID(val_16);
[val_18_pos,val_18_neg] = PosNeg_system_ID(val_18);
[val_20_pos,val_20_neg] = PosNeg_system_ID(val_20);
[val_22_pos,val_22_neg] = PosNeg_system_ID(val_22);

% in case of debugging, the options for the 'compare' command
opt_comp = compareOptions('InitialCondition','z');

% set options for the SS models: zero initial condition (x_0 = 0)
% and focus on simulation.
if(isfield(rOptions,'tFocus'))
    tFocus = rOptions.tFocus;
else
    tFocus = 'simulation';
end
if(isfield(rOptions,'tSearchMethod'))
    tSearchMethod = rOptions.tSearchMethod;
else
    tSearchMethod = 'auto';
end

%% Armax
if(strcmp(rOptions.tMethod,'ARMAX'))
    % allow noise integration?
    if(isfield(rOptions,'bIN'))
        bIN = rOptions.bIN;
    else
        bIN = 0;
    end
    opt = armaxOptions('InitialCondition','zero','Focus',tFocus,...
        'SearchMethod',tSearchMethod);
    
    % select the orders for ARMAX (na,bbs,nc,nks)
    if(isfield(rOptions,'myOrder'))
        myOrder = rOptions.myOrder;
    else
        myOrder = 4;
    end
    
    if(isfield(rOptions,'myOrderNoise'))
        myOrderNoise = rOptions.myOrderNoise;
    else
        myOrderNoise = 4;
    end
    
    na = myOrder;
    nc = myOrderNoise;
    % if there are 2 inputs present (SPV + exogenous) need vector for nb,nk
    if(size(exp_4.InputName,1) == 2)
        nbs = [myOrder myOrder] + 3.*ones(1,2);
        nks = [0 0];
    else
        nbs = myOrder;
        nks = 0;
    end
    
    %% Steady-state
elseif(strcmp(rOptions.tMethod,'SS'))
    % select if feedthrough allowed (default is 0)
    if(isfield(rOptions,'bFeedthrough'))
        bFeedthrough = rOptions.bFeedthrough;
    else
        bFeedthrough = 0;
    end
    opt = ssestOptions('InitialState','zero','Focus',tFocus,...
        'SearchMethod',tSearchMethod);
    
    % select the orders for SS (na,bbs,nc,nks)
    if(isfield(rOptions,'myOrder'))
        myOrder = rOptions.myOrder;
    else
        myOrder = 5;
    end
    
    nx = myOrder;
end

%% 2) Analysis for positive inputs (u>0)
if(strcmp(rOptions.tMethod,'ARMAX'))
    my_mod_4_pos = armax(exp_4_pos,[na,nbs,nc,nks],'IntegrateNoise',bIN,opt);
    my_mod_7_pos = armax(exp_7_pos,[na,nbs,nc,nks],'IntegrateNoise',bIN,opt);
    my_mod_12_pos = armax(exp_12_pos,[na,nbs,nc,nks],'IntegrateNoise',bIN,opt);
    my_mod_16_pos = armax(exp_16_pos,[na,nbs,nc,nks],'IntegrateNoise',bIN,opt);
    my_mod_18_pos = armax(exp_18_pos,[na,nbs,nc,nks],'IntegrateNoise',bIN,opt);
    my_mod_20_pos = armax(exp_20_pos,[na,nbs,nc,nks],'IntegrateNoise',bIN,opt);
    my_mod_22_pos = armax(exp_22_pos,[na,nbs,nc,nks],'IntegrateNoise',bIN,opt);
elseif(strcmp(rOptions.tMethod,'SS'))
    my_mod_4_pos = ssest(exp_4_pos,nx,'Feedthrough',bFeedthrough,opt);
    my_mod_7_pos = ssest(exp_7_pos,nx,'Feedthrough',bFeedthrough,opt);
    my_mod_12_pos = ssest(exp_12_pos,nx,'Feedthrough',bFeedthrough,opt);
    my_mod_16_pos = ssest(exp_16_pos,nx,'Feedthrough',bFeedthrough,opt);
    my_mod_18_pos = ssest(exp_18_pos,nx,'Feedthrough',bFeedthrough,opt);
    my_mod_20_pos = ssest(exp_20_pos,nx,'Feedthrough',bFeedthrough,opt);
    my_mod_22_pos = ssest(exp_22_pos,nx,'Feedthrough',bFeedthrough,opt);
end

disp(['======================= POSITIVE INPUT (order ',num2str(myOrder),...
    ')========================']);

f_v_4_pos = Validate_model(val_4_pos,my_mod_4_pos); close;
f_v_7_pos = Validate_model(val_7_pos,my_mod_7_pos); close;
f_v_12_pos = Validate_model(val_12_pos,my_mod_12_pos); close;
f_v_16_pos = Validate_model(val_16_pos,my_mod_16_pos); close;
f_v_18_pos = Validate_model(val_18_pos,my_mod_18_pos); close;
f_v_20_pos = Validate_model(val_20_pos,my_mod_20_pos); close;
f_v_22_pos = Validate_model(val_22_pos,my_mod_22_pos); close;

%% 3) Analysis for negative inputs (u<0)
if(strcmp(rOptions.tMethod,'ARMAX'))
    my_mod_4_neg = armax(exp_4_neg,[na,nbs,nc,nks],'IntegrateNoise',bIN,opt);
    my_mod_7_neg = armax(exp_7_neg,[na,nbs,nc,nks],'IntegrateNoise',bIN,opt);
    my_mod_12_neg = armax(exp_12_neg,[na,nbs,nc,nks],'IntegrateNoise',bIN,opt);
    my_mod_16_neg = armax(exp_16_neg,[na,nbs,nc,nks],'IntegrateNoise',bIN,opt);
    my_mod_18_neg = armax(exp_18_neg,[na,nbs,nc,nks],'IntegrateNoise',bIN,opt);
    my_mod_20_neg = armax(exp_20_neg,[na,nbs,nc,nks],'IntegrateNoise',bIN,opt);
    my_mod_22_neg = armax(exp_22_neg,[na,nbs,nc,nks],'IntegrateNoise',bIN,opt);
elseif(strcmp(rOptions.tMethod,'SS'))
    my_mod_4_neg = ssest(exp_4_neg,nx,'Feedthrough',bFeedthrough,opt);
    my_mod_7_neg = ssest(exp_7_neg,nx,'Feedthrough',bFeedthrough,opt);
    my_mod_12_neg = ssest(exp_12_neg,nx,'Feedthrough',bFeedthrough,opt);
    my_mod_16_neg = ssest(exp_16_neg,nx,'Feedthrough',bFeedthrough,opt);
    my_mod_18_neg = ssest(exp_18_neg,nx,'Feedthrough',bFeedthrough,opt);
    my_mod_20_neg = ssest(exp_20_neg,nx,'Feedthrough',bFeedthrough,opt);
    my_mod_22_neg = ssest(exp_22_neg,nx,'Feedthrough',bFeedthrough,opt);
end

disp(['======================= NEGATIVE INPUT (order ',num2str(myOrder),...
    ')========================']);

f_v_4_neg = Validate_model(val_4_neg,my_mod_4_neg); close;
f_v_7_neg = Validate_model(val_7_neg,my_mod_7_neg); close;
f_v_12_neg = Validate_model(val_12_neg,my_mod_12_neg); close;
f_v_16_neg = Validate_model(val_16_neg,my_mod_16_neg); close;
f_v_18_neg = Validate_model(val_18_neg,my_mod_18_neg); close;
f_v_20_neg = Validate_model(val_20_neg,my_mod_20_neg); close;
f_v_22_neg = Validate_model(val_22_neg,my_mod_22_neg); close;

%% 4) Store for output

% store the fit mean and std values for output (POSITIVE)
meanFit(1,1) = mean(f_v_4_pos); stdFit(1,1) = std(f_v_4_pos);
meanFit(2,1) = mean(f_v_7_pos); stdFit(2,1) = std(f_v_7_pos);
meanFit(3,1) = mean(f_v_12_pos); stdFit(3,1) = std(f_v_12_pos);
meanFit(4,1) = mean(f_v_16_pos); stdFit(4,1) = std(f_v_16_pos);
meanFit(5,1) = mean(f_v_18_pos); stdFit(5,1) = std(f_v_18_pos);
meanFit(6,1) = mean(f_v_20_pos); stdFit(6,1) = std(f_v_20_pos);
meanFit(7,1) = mean(f_v_22_pos); stdFit(7,1) = std(f_v_22_pos);

% store the fit mean and std values for output (NEGATIVE)
meanFit(1,2) = mean(f_v_4_neg); stdFit(1,2) = std(f_v_4_neg);
meanFit(2,2) = mean(f_v_7_neg); stdFit(2,2) = std(f_v_7_neg);
meanFit(3,2) = mean(f_v_12_neg); stdFit(3,2) = std(f_v_12_neg);
meanFit(4,2) = mean(f_v_16_neg); stdFit(4,2) = std(f_v_16_neg);
meanFit(5,2) = mean(f_v_18_neg); stdFit(5,2) = std(f_v_18_neg);
meanFit(6,2) = mean(f_v_20_neg); stdFit(6,2) = std(f_v_20_neg);
meanFit(7,2) = mean(f_v_22_neg); stdFit(7,2) = std(f_v_22_neg);

% store all the models
Models.my_mod_4_pos = my_mod_4_pos; Models.my_mod_4_neg = my_mod_4_neg;
Models.my_mod_7_pos = my_mod_7_pos; Models.my_mod_7_neg = my_mod_7_neg;
Models.my_mod_12_pos = my_mod_12_pos; Models.my_mod_12_neg = my_mod_12_neg;
Models.my_mod_16_pos = my_mod_16_pos; Models.my_mod_16_neg = my_mod_16_neg;
Models.my_mod_18_pos = my_mod_18_pos; Models.my_mod_18_neg = my_mod_18_neg;
Models.my_mod_20_pos = my_mod_20_pos; Models.my_mod_20_neg = my_mod_20_neg;
Models.my_mod_22_pos = my_mod_22_pos; Models.my_mod_22_neg = my_mod_22_neg;

end