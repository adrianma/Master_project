function Model = precompute_system_model(Params)
%% ===================================================================
%  Evangelos Vrettos
%  April 2012
%
%  Modified in October 2014, for Powertech paper simulations
%  ===================================================================

%% Model Structure
% Ad = inv(A1 + á*åeff*A2 + V*A3) * (B1 + á*åeff*B2 + V*B3)
% Bd = inv(A1 + á*åeff*A2 + V*A3) * C
% x_dot = Ad * x + Bd * u
% Purpose of this file: precompute A1, A2, A3 and B1, B2, B3 and C (for each EWH) that are time independent.

%% Preallocate
Model.A1     = nan(Params.n+2, Params.n+2, Params.n_app);
Model.A2     = nan(Params.n+2, Params.n+2, Params.n_app);
Model.A3     = nan(Params.n+2, Params.n+2, Params.n_app);
Model.B1     = nan(Params.n+2, Params.n+2, Params.n_app);
Model.B2     = nan(Params.n+2, Params.n+2, Params.n_app);
Model.B3     = nan(Params.n+2, Params.n+2, Params.n_app);
Model.C      = nan(Params.n+2, 2, Params.n_app);

%% Construct the model
for i = 1:Params.n_app,
    %% Model matrix components that can be precomputed
    dt = Params.t_sample;
    dx = Params.d(:,:,i);
        
    % heat loss coefficients
    A1_loss = Params.diam(:,:,i) * pi * Params.h(:,:,i);
    A2_loss = Params.A_outer_i(:,:,i)+Params.A_cross_i(:,:,i);
    k_coeff1 = A1_loss * Params.alpha(:,:,i) / Params.m(:,:,i) / Params.c;              % heat loss coefficient for disks 3,4...N-2
    k_coeff2 = A2_loss * Params.alpha(:,:,i) / Params.m_i(:,:,i) / Params.c;            % heat loss coefficient for disks 2,N-1
    
    %% Calculate matrices A1, A2, A3 and B1, B2 and B3
    % Correspond to diffusion and ambient losses terms. Turbulent mixing during water draw and forced convection are NOT precomputed.
    % A1: does not depend neither on åeff nor on V. A2: depends only on
    % åeff. A3 depends only on V.
    
    % Boundary layers of Matrix A1
    A1 = zeros(Params.n+2);
    
    A1(1,1) = 1;                                        % first line of A1 corresponds to the artificial disk 1, mapping to inlet water temperature
    A1(Params.n+2,Params.n+2) = 1;                      % last line of A1 corresponds to the artificial disk N, mapping to ambient temperature
    
    % The 2 and the N-1 line of A1 correspond to part of the boundary conditions (disks 2 and N-1)
    A1(2,1) = 0;                                        % A(2,1) corresponds to disk 2. Convection speed for disk 2 is given by U_coeff(1), but not used in A1.
    A1(2,2) = 1/dt;
    A1(2,Params.n+2) = -k_coeff2;
    A1((Params.n+2)-1,(Params.n+2)-1) = 1/dt;
    A1((Params.n+2)-1,Params.n+2) = -k_coeff2;
    
    % Boundary layers of Matrix A2
    A2 = zeros(Params.n+2);
    
    % Boundary layers of Matrix A3
    A3 = zeros(Params.n+2);
    A3(2,1) = -1/dx;                                    % A(2,1) corresponds to disk 2. Convection speed for disk 2 is given by U_coeff(1).
    
    % Boundary layers of Matrix B1
    B1 = zeros(Params.n+2);
    
    % The artificial temperatures at disks 1 and N+2 do not change
    B1(1,1) = 1;
    B1(Params.n+2,Params.n+2) = 1;
    
    % The 2 and the N-1 line of B1 correspond to the boundary conditions (disks 2 and N-1)
    C_coeff1 = (1/dt)-k_coeff2;
    B1(2,2) = C_coeff1;
    
    % The disk N-1 corresponds to the last indice of a_coeff, which is N-2 or  (Params.n+2)-2
    C_coeff2 = (1/dt)-k_coeff2;
    B1((Params.n+2)-1,(Params.n+2)-1) = C_coeff2;
    
    % Boundary layers of Matrix B2
    B2 = zeros(Params.n+2);
    
    % The 2 and the N-1 line of B correspond to the boundary conditions (disks 2 and N-1)
    B2(2,2) = -1/(dx^2);
    B2(2,3) = 1/(dx^2);
    % The disk N-1 corresponds to the last indice of a_coeff, which is N-2 or  (Params.n+2)-2
    B2((Params.n+2)-1,(Params.n+2)-1) = -1/(dx^2);
    B2((Params.n+2)-1,(Params.n+2)-2) = 1/(dx^2);
    
    % Boundary layers of Matrix B3
    B3 = zeros(Params.n+2);
    
    % The 2 and the N-1 line of B correspond to the boundary conditions (disks 2 and N-1)
    B3(2,2) = -1/dx;
    % The disk N-1 corresponds to the last indice of a_coeff, which is N-2 or  (Params.n+2)-2
    B3((Params.n+2)-1,(Params.n+2)-1) = -1/dx;
    B3((Params.n+2)-1,(Params.n+2)-2) = 1/dx;
    
    % All the other lines of A1, A2, A3 have similar structure
    % For speed, in the same loop the entries of Matrices B1, B2, B3 are also calculated
    % First for A1 and B1
    for idx=3:1:((Params.n+2)-2)
        A_coeff1 = 0;
        A_coeff2 = 0;
        A_coeff3 = 0;
        B_coeff1 = 0;
        B_coeff3 = 0;
        Coefficients_A = [((-A_coeff1)/2)-B_coeff1 (1/dt)+A_coeff2+(k_coeff1/2) B_coeff3-(A_coeff3/2)];
        Coefficients_B = [((A_coeff1/2)+B_coeff1) ((1/dt)-A_coeff2-(k_coeff1/2)) ((A_coeff3/2)-B_coeff3)];
        k = 1;
        for j=idx-1:1:((idx-1)+2)
            A1(idx,j) = Coefficients_A(k);
            B1(idx,j) = Coefficients_B(k);
            k = k+1;
        end
        A1(idx,Params.n+2) = -k_coeff1;
    end
    
    % Second for A2 and B2
    for idx=3:1:((Params.n+2)-2)
        A_coeff1 = 1/(dx^2);
        A_coeff2 = 1/(dx^2);
        A_coeff3 = 1/(dx^2);
        B_coeff1 = 0;
        B_coeff3 = 0;
        Coefficients_A = [((-A_coeff1)/2)-B_coeff1 A_coeff2 B_coeff3-(A_coeff3/2)];
        Coefficients_B = [((A_coeff1/2)+B_coeff1) -A_coeff2 ((A_coeff3/2)-B_coeff3)];
        k = 1;
        for j=idx-1:1:((idx-1)+2)
            A2(idx,j) = Coefficients_A(k);
            B2(idx,j) = Coefficients_B(k);
            k = k+1;
        end
    end
    
    % Third for A3 and B3
    for idx=3:1:((Params.n+2)-2)
        A_coeff1 = 0;
        A_coeff2 = 0;
        A_coeff3 = 0;
        B_coeff1 = 1/(4*dx);
        B_coeff3 = 1/(4*dx);
        Coefficients_A = [((-A_coeff1)/2)-B_coeff1 A_coeff2 B_coeff3-(A_coeff3/2)];
        Coefficients_B = [((A_coeff1/2)+B_coeff1) -A_coeff2 ((A_coeff3/2)-B_coeff3)];
        k = 1;
        for j=idx-1:1:((idx-1)+2)
            A3(idx,j) = Coefficients_A(k);
            B3(idx,j) = Coefficients_B(k);
            k = k+1;
        end
    end
    
    %% Calculate matrix C
    % Corresponds to internal heating by heating elements
    C = zeros(Params.n+2,2);
    
    loc1 = Params.element1_location;
    mix1 = Params.mixing_zone1_size;
    ratio1 = Params.mixing_zone1_ratio;
    % distribute mixing zone around the location of the element element. Location should be in the middle of the mixing zone
    % if uneven number, put equal amount above and below; if even number, put one more above than below
    upper_elements1 = floor((mix1-1)/2);
    lower_elements1 = ceil((mix1-1)/2);
    C_mixing_idx1 = [loc1-lower_elements1:loc1+upper_elements1];
    % C_mixing_idx1 = [loc1:Params.n+1];
    loc2 = Params.element2_location;
    mix2 = Params.mixing_zone2_size;
    ratio2 = Params.mixing_zone2_ratio;
    % distribute mixing zone around the location of the element element. Location should be in the middle of the mixing zone
    % if uneven number, put equal amount above and below; if even number, put one more above than below
    upper_elements2 = floor((mix2-1)/2);
    lower_elements2 = ceil((mix2-1)/2);
    C_mixing_idx2 = [loc2-lower_elements2:loc2+upper_elements2];
    
    % Calculate the heating power that affects each of the layers within the mixing zone
    for idx = 1:mix1
        pos = C_mixing_idx1(idx);
        if pos == loc1
            C(pos,1) =  ratio1*(Params.P1_el(:,:,i) * Params.eta /Params.m_i(:,:,i) / Params.c * 1000);
        else
            fac = (1-ratio1)/(mix1-1);
            C(pos,1) = fac*(Params.P1_el(:,:,i) * Params.eta /Params.m_i(:,:,i) / Params.c * 1000);
        end
    end
    
    for idx = 1:mix2
        pos = C_mixing_idx2(idx);
        if pos == loc2
            C(pos,2) =  ratio2*(Params.P2_el(:,:,i) * Params.eta /Params.m_i(:,:,i) / Params.c * 1000);
        else
            fac = (1-ratio2)/(mix2-1);
            C(pos,2) = fac*(Params.P2_el(:,:,i) * Params.eta /Params.m_i(:,:,i) / Params.c * 1000);
        end
    end
    
    %% Store A, B and C matrix parameters in Model. structure
    Model.A1(:,:,i) = A1;
    Model.A2(:,:,i) = A2;
    Model.A3(:,:,i) = A3;
    
    Model.B1(:,:,i) = B1;
    Model.B2(:,:,i) = B2;
    Model.B3(:,:,i) = B3;
    
    Model.C(:,:,i) = C;
end
