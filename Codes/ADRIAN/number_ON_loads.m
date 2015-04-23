function Res = number_ON_loads(Params,Results,event_step)

AAA_1 = squeeze(Results.xrec(Params.sensor1_location,:,:));
AAA_2 = squeeze(Results.xrec(Params.sensor2_location,:,:));

Res.sInDeadband = 0;
for ii = 1:Params.n_app
    
    % see the temperature
    T_meas1 = AAA_1(event_step,ii);
    T_meas2 = AAA_2(event_step,ii);
    
    % decide wether in device is in deadband or not
    bInDeadband1 = (T_meas1 <= Params.T_max1(:,:,ii))&&...
        (T_meas1 >= Params.T_min1(:,:,ii))&&Params.element1_present;
    bInDeadband2 = (T_meas2 <= Params.T_max2(:,:,ii))&&...
        (T_meas2 >= Params.T_min2(:,:,ii))&&Params.element2_present;
    bInDeadband = bInDeadband1 || bInDeadband2;
    
    if(bInDeadband)
        Res.sInDeadband = Res.sInDeadband + 1;
    end
end

Res.N_ON_1 = sum(Results.urec(1,event_step,:) == 1);
Res.N_ON_2 = sum(Results.urec(2,event_step,:) == 1);


end