function mean_visc=viscosity(T)

% Calculation based on Arrhenius equation and the parameter values taken
% from: Dependence of Water Viscosity on Temperature and Pressure
% E. R. Likhachev

n0=2.4055*10^(-5); % Pa*sec
E=4.753*10^3; % J/mol
T0=-139.7; % K
R=8.31441; % J/(mol*K)

visc=n0.*exp(E./(R.*(T+273.15+T0)));

% save up some computational time by computing the mean this way
mean_visc=sum(visc)/length(visc);