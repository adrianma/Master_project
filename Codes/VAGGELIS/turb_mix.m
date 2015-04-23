function turb=turb_mix(p_dens,U_coeff,mean_visc,beta_thermal,T,Params,i)

% Calculation based on the following:
% Model by Zurigat et al found at Heat Transfer and Stratification in
% Sensible Heat Storage Systems

if U_coeff(1)>0
    lr1=0.018; % inlet diameter, taken from the previous book page 274
    lr2=Params.h(:,:,i); % Tank height
    p_mean=mean(p_dens);
    b_mean=mean(beta_thermal);
    T0=Params.x_init(1,1,i);
    DT=T-T0;
    mean_DT=mean(DT);

    % Calculation of flow speed in the inlet pipe with the assumption of lr1
    ur=(U_coeff(1)*(Params.diam(:,:,i)/2)^2)/((lr1/2)^2);

    % Calculation of the turbulent diffusivity at the inlet
    Re=p_mean*ur*(lr1/mean_visc);
    Ri=Params.g*b_mean*mean_DT*(lr2/ur^2);
    
    scale=0.9;
%     e_eff_in=scale*4.75*(Re/Ri)^0.522;
    e_eff_in=4.75*(10000)^0.522;

    % Calculation of the turbulent diffusivity at every element of the tank
%     A=(e_eff_in-1)/(1-(1/Params.n));
%     B=e_eff_in-A;
% 
%     x=(1:1:Params.n);
%     e_eff=(A./x)+B;
    
%     x=(1:1:Params.n);
%     a=(1/(Params.n-1))*log(e_eff_in+1);
%     e_eff=(e_eff_in+1)-exp(a*(x-1));
% 
x=(1:1:Params.n);
a=((e_eff_in-1)/(Params.n-1));
e_eff=(a*(Params.n-x))+1;

% x=(0:1/Params.n:1);
% e_eff=e_eff_in .* sqrt(1 - x .^ 2);
% e_eff=e_eff(2:end);
% e_eff(e_eff<1)=1;

    turb=e_eff;
else
    turb=1;
end

% x=0:1:100;
% y=10-sqrt(x);
% figure(2);clf;hold on;
% plot(x,y);

