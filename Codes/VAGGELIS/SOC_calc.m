function SOC = SOC_calc(Params, i, xrec)
%% =====================================================================
%  Calculation of EWH SOC based on the temperature distribution (taken as
%  input in vector T).
%
%  Evangelos Vrettos
%  April 2012
%
%  Modified by Evangelos Vrettos
%  October 2014
%  =====================================================================

T = xrec(2:end-1,:,i);

if Params.SOCdef == 1
    pos = Params.sensor1_location-1;           % subtract 1 because the fake layers are disregarded
    tmp = (T(pos,:,i)-Params.T_min1(:,:,i))/(Params.T_max1(:,:,i)-Params.T_min1(:,:,i));
    SOC = min([1 max([0 tmp])]);
elseif Params.SOCdef == 2
    % Same weights to all layers
    w = ones(1,10);
    
    % Translate temperatures from Celsius to Kelvin
    T = T + 273.15;
    T_min = Params.T_min1(:,:,i) + 273.15;
    T_max = Params.T_max1(:,:,i) + 273.15;
    
    % Calculate the reference thermal content
    num = Params.n * Params.c * (T_max - T_min) * sum(Params.m_i(:,:,i)*w);
    den = sum(w);
    Q0 = (num/den);
        
    % Initialize
    Q = zeros(Params.n,1);
    
    % Only the useful layers contribute to the EWH SOC
    Useful_layers = find(T > T_min);
    Q(Useful_layers) = Params.m_i(:,:,i) * Params.c .* (T(Useful_layers) - T_min);
   
    num = Params.n * sum(w*Q);
    den = sum(w);
    Q_tot = num/den;
    
    tmp = Q_tot/Q0;
    SOC = min([1 tmp]);
elseif Params.SOCdef == 3
    % Same weights to all layers
    w = ones(1,10);
    
    % Translate temperatures from Celsius to Kelvin
    T = T + 273.15;
    T_min = Params.T_min1(:,:,i) + 273.15;
    T_max = Params.T_max1(:,:,i) + 273.15;
    
    % Calculate the reference thermal content
    num = Params.n * Params.c * (T_max - T_min) * sum(Params.m_i(:,:,i)*w);
    den = sum(w);
    Q0 = (num/den);
           
    % All layers contribute to the EWH SOC
    Q = Params.m_i(:,:,i) * Params.c .* (T - T_min);
   
    num = Params.n * sum(w*Q);
    den = sum(w);
    Q_tot = num/den;
    
    tmp = Q_tot/Q0;
    SOC = min([1 max([0 tmp])]);
end



end
