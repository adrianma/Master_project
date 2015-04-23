function CompleteModel = build_final_system_model(Params,Model,m_dot)
%% ===================================================================
%  Assemble the final A-Matrix (Ad) and B-Matrix (Bd) based on precomputed
%  model.
%  Ad = inv(A1 + ?*?eff*A2 + V*A3) * (B1 + ?*?eff*B2 + V*B3)
%  Bd = inv(A1 + ?*?eff*A2 + V*A3) * C
%  x_dot = Ad * x + Bd * u
%
%  Evangelos Vrettos
%  April, 2012
%
%  Modified by Evangelos Vrettos
%  October 2014
%  ===================================================================

CompleteModel.Ad = nan(size(Model.A1));
CompleteModel.Bd = nan(size(Model.C));

% Calculate conductivity, density and diffusivity for disks 2...N-1 (not for 1 and N)
p_dens = Params.rho * ones(Params.n,1);
a_coeff = Params.k / (Params.c * Params.rho) * ones(Params.n,1);

for i = 1:Params.n_app
    %% Calculate the effect of turbulent mixing on thermal diffusivity
    % Thermal expansion coefficient
    dx = Params.d(:,:,i);
    T = Params.x_current(:,:,i);
    beta_thermal = thermal_coeff(T(2:end-1),p_dens);
    
    % Convection speed for disks 2...N-1
    U_coeff = (m_dot(i,:)./ (pi .* p_dens .* ((Params.diam(:,:,i)/2)^2))); % in m/s
    
    % viscosity
    mean_visc = viscosity(T(2:end-1));
    
    % Turbulent diffusivity calculation
    if strcmp(Params.turbMixing,'Yes')
        turb = turb_mix(p_dens,U_coeff,mean_visc,beta_thermal,T,Params,i);
        turb = turb';
    elseif strcmp(Params.turbMixing,'No')
        turb = 1;                                   % it turns out that turbulent diffusivity causes numerical instability. It is ok only for short time intervals.
    end
    
    % Modify thermal diffusivity in case turbulent diffusivity is considered
    a_coeff_tmp = a_coeff .* turb;
    a_coeff_2 = [0; a_coeff_tmp; 0];                % Adds zeros to the fake disks
    dif_matrix = diag(a_coeff_2);
    
    %%  Build temporary matrices
    A = Model.A1(:,:,i) + Model.A2(:,:,i) * dif_matrix + U_coeff(1) * Model.A3(:,:,i);
    
    tmp = Model.B2(:,:,i) * dif_matrix;
    tmp(2,3) = a_coeff(1)/(dx^2);                   % Correct entries for boundary layers (only needed for B)
    tmp((Params.n+2)-1,(Params.n+2)-2) = a_coeff((Params.n+2)-2)/(dx^2);
    B = Model.B1(:,:,i) + tmp + U_coeff(1) * Model.B3(:,:,i);
    
    %% Build complete model matrices
    Ad = A\B;                             % inv(A) * B;
    Bd = A\(Model.C(:,:,i));              % inv(A) * C;

    CompleteModel.Ad(:,:,i) = Ad;
    CompleteModel.Bd(:,:,i) = Bd;
    
end