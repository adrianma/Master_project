%% ========================================================================
%  function solve_CQP_problem
%  by Adrian Martinez Gomez
%  December 2014
%
%  Modified by Adrian Martinez Gomez
%  March 2015
%
%  Output: Results   - struct where the results are stored
%  Input: ss         - state space system
%         N_horizon  - number of time steps for the horizon
%         bExogenous - boolean, usually select as 0 to handle exogenous
%                      input correctly
%         ref_signal - reference signal (usually constructed)
%
%  Purpose:
%   * calculate the solutions for the GUROBI solver for 1 time step with
%   MPC and (C)QP.
%   * matrices are defined here (this is the most difficult part by far)
%  ========================================================================
function Results = solve_CQP_problem(ss,x_current,N_horizon,bExogenous,...
    rRefining,ref_signal,exogenous_input,global_counter);

% define the time horizon
Results.N_horizon = N_horizon;

% get all necessary dimensions
n_states = size(ss.a,1);
n_inputs = size(ss.b,2);
% remove the exogenous input if it is not selected
if(~bExogenous)
    n_inputs = n_inputs - 1;
end

Results.n_z = n_states + n_inputs;

%% 1) Construct the matrices
%
% define terminal cost matrix (it has to be positive definite)
P_terminal = rRefining.p_terminal.*eye(n_states);

Matrices = construct_matrices(ss,x_current,P_terminal,N_horizon,n_states,...
    n_inputs,bExogenous,ref_signal,exogenous_input,rRefining,...
    global_counter);
Results.Matrices = Matrices;

%% 2) Solve the Gurobi problem
%
% solve the problem with gurobi
Res = solve_gurobi_problem(Matrices,N_horizon,n_states,n_inputs,...
    ref_signal,rRefining,global_counter);
Results.gurobi_res = Res;

% Gurobi didn't find a solution; plot the necessary things and terminate
if(~isfield(Res,'x'))
    disp('  ==> Gurobi didnt find a solution <==');
    
    figure;
    plot(ref_signal);
    error(' terminating...(see figure to know what ref_signal was)')
end

% 2.1) extract the wanted information
xpred = cell(1,N_horizon+1);
for ii = 1:N_horizon+1
    % extract the n_states elements (nth order in this case, 3,4)
    xpred{ii} = Res.x(1+(ii-1)*Results.n_z :(ii-1)*Results.n_z+n_states);
end

upred = cell(1,N_horizon);
for ii = 1:N_horizon
    % extract the n_inputs elements (1 in this case)
    upred{ii} = Res.x(1+(ii-1)*Results.n_z+n_states:(ii-1)*Results.n_z+...
        Results.n_z);
end

% store the predicted states and inputs in the 'Results' struct
Results.xpred = xpred;
Results.upred = upred;

end

%% ------------------------------------------------------------------------
%  function construct_matrices
%  by Adrian Martinez Gomez
%  December 2014
%
%  Output: 'Matrices' is a struct containing all the matrices
%
%  Purpose:
%   * construct the (augmented) matrices for the CQP problem
%  ------------------------------------------------------------------------
function Matrices = construct_matrices(ss,x_current,P_terminal,...
    N_horizon,n_states,n_inputs,bExogenous,ref_signal,exogenous_input,...
    rRefining,global_counter);

% define the augmented number of states n_z
n_z = n_states + n_inputs;

A = ss.a;
% IMPORTANT NOTE:
% we have to take into account whether the exogenous input (be as it may)
% has any influence on the system. Here it is assumed NOT (for now)
if(~bExogenous)
    B = ss.b(:,1:end-1);
else
    B = ss.b;
end

%% (I) matrix for the quadratic cost
% submatrix for Q
Q = [ss.c.'*ss.c,ss.c'*ss.d(1:n_inputs);ss.d(1:n_inputs).'*ss.c,...
    ss.d(1:n_inputs).'*ss.d(1:n_inputs)];

% from 0...N_horizon-1
Matrices.Q_bar = blkdiag(kron(eye(N_horizon),Q),P_terminal);
% scale by the coefficient
Matrices.Q_bar = rRefining.CostFunctionCoeff(1).*Matrices.Q_bar;

% see if this is a Convex QP (CQP) or not
%  (see if matrix Q positive semi-definite by checking evals)
if(sum(svd(Matrices.Q_bar) >= 0) == (n_z*N_horizon + n_states))
    disp(num2str(' => Convex QP'));
else
    disp(num2str(' => Non-Convex QP'));
end

%% (II) matrix (vector) for the linear cost
R_bar = cell(N_horizon,1);
for ii = 1:N_horizon
    % stack the results for the R_bar (linear cost term for horizon N)
    R_bar{ii} = [-2*ref_signal((global_counter-1) + ii).*ss.c,...
        -2*ref_signal((global_counter-1) + ii).*ss.d(1:n_inputs)].';
end

% stack the last zeros, and store it on the final struct
Matrices.R_bar = [cell2mat(R_bar);zeros(n_states,1)];
% scale by the coefficient
Matrices.R_bar = rRefining.CostFunctionCoeff(1).*Matrices.R_bar;

%% (III) matrices for the equality constraints
% pre-allocate
Matrices.A_eq = nan((N_horizon + 1)*n_states,N_horizon*n_z + n_states);
cell_A_eq = cell(1,n_inputs);
for ii = 1:(N_horizon + 1)
    % construct the matrix A_eq column-wise; storing and dimensions
    % explained on paper
    switch ii
        % 1st block-column has no zeros ABOVE the 4 matrices
        case 1
            cell_A_eq{ii} = [eye(n_states),zeros(size(B));A,B;...
                zeros((N_horizon + 1 - 2 - (ii-1))*n_states,n_z)];
            % last block-column has no zeros BELOW  the identity matrix
        case (N_horizon + 1)
            cell_A_eq{ii} = [zeros(N_horizon*n_states,n_states);...
                -eye(n_states)];
            % in this cases there are zeros below and above the 4 matrices
        otherwise
            cell_A_eq{ii} = [zeros((ii-1)*n_states,n_z);-eye(n_states)...
                ,zeros(size(B));A,B;...
                zeros((N_horizon + 1 - 2 - (ii-1))*n_states,n_z)];
    end
    if(size(cell_A_eq{ii},1) ~= (N_horizon + 1)*n_states)
        disp(['Dimension mismatch at iteration number ',num2str(ii)]);
    end
end
% store the matrix 
Matrices.A_eq = cell2mat(cell_A_eq);

% right-hand side
%
% pre-allocate
Matrices.b_eq = nan((N_horizon+1)*n_states,1);

% (take into account the exogenous input)
if(~bExogenous)
    temp_RHS = cell(N_horizon,1);
    for ii = 1:N_horizon
        % indexing of the exogenous_input needs the global_counter variable
        % to update after every global iteration
        temp_RHS{ii} = -ss.b(:,2).*exogenous_input((global_counter-1)+ ii);
    end
    Matrices.b_eq = [x_current;cell2mat(temp_RHS)];
    
else
    Matrices.b_eq = [x_current;repmat(zeros(n_states,1),N_horizon,1)];
end

%% (V) bounds for augmented state z
% %
% % pre-allocate
% % actually, it's not needed
% % ...(change)
% Matrices.lower_bounds = nan(n_z*N_horizon + n_states,1);
% Matrices.upper_bounds = nan(n_z*N_horizon + n_states,1);

% assume the bounds don't change over time, and this is for 1 time step
LB_mat = [-inf(n_states,1);rRefining.u_min];
UB_mat = [inf(n_states,1);rRefining.u_max];

% store the bounds
Matrices.lower_bounds = [repmat(LB_mat,N_horizon,1);-inf(n_states,1)];
Matrices.upper_bounds = [repmat(UB_mat,N_horizon,1);inf(n_states,1)];

end

%% ------------------------------------------------------------------------
%  function solve_gurobi_problem
%  by Adrian Martinez Gomez
%  December 2014
%
%  Output: 'Res' is the result of the GUROBI solver; 'Model' is the
%  computed model for the (C)QP problem
%
%  Purpose:
%   * solve the problem with the GUROBI solver
%  ------------------------------------------------------------------------
function [Res,Model] = solve_gurobi_problem(Matrices,N_horizon,n_states,...
    n_inputs,ref_signal,rRefining,global_counter);

n_z = n_states + n_inputs;

% affine term of the reference output signal tracking
affine_term = sum(ref_signal(1 + (global_counter-1):...
    1 + (global_counter-1) + N_horizon).^2);
% scale by the coefficient
affine_term = rRefining.CostFunctionCoeff(1).*affine_term;

%% Adapt the cost function if needed
% select this snippet to adapt the cost function 
if(rRefining.CostFunctionCoeff(2) ~= 0)
    % ===== 1) penalize input with SUM_k=0...N-1 {(u_k)^2} =====
    sub_SIGMA = [zeros(n_states),zeros(n_states,n_inputs);...
        zeros(n_inputs,n_states),eye(n_inputs)];
    
    SIGMA = blkdiag(kron(eye(N_horizon),sub_SIGMA),zeros(n_states));
    
    Matrices.Q_bar = Matrices.Q_bar + rRefining.CostFunctionCoeff(2).*...
        SIGMA.'*SIGMA;
    
elseif(rRefining.CostFunctionCoeff(3) ~= 0)
    % ===== 2) penalize difference between consecutive inputs =====
    
    % starting and ending block matrices    
    sub_LAMBDA = [zeros(n_states),zeros(n_states,n_inputs);...
        zeros(n_inputs,n_states),eye(n_inputs)];
    
    % the block matrices of: the diagonal (diag), above the diagonal 
    % (upper), and below the diagonal (lower)
    diag_LAMBDA = bsxfun(@times,sub_LAMBDA,...
        permute([1,repmat(2,1,N_horizon-2),1],[3,1,2]));
    upper_LAMBDA = bsxfun(@times,sub_LAMBDA,...
        permute(repmat(-2,1,N_horizon-1),[3,1,2]));
    lower_LAMBDA = zeros(size(upper_LAMBDA));
    
    % use the command 'blocktridiag', which allows to construct such
    % matrices
    %
    % look at: http://bit.ly/1AJzz6Z
    LAMBDA = blktridiag(diag_LAMBDA,lower_LAMBDA,upper_LAMBDA);
    LAMBDA = full(LAMBDA);
    LAMBDA = blkdiag(LAMBDA,zeros(n_states));
    
    Matrices.Q_bar = Matrices.Q_bar + ...
        rRefining.CostFunctionCoeff(3).*LAMBDA;
    
elseif(rRefining.CostFunctionCoeff(4) ~= 0)
    % ===== 3) penalize input with SUM_k=0...N-1 {(u_k - u_max)^2} =====
    
    % first: the quadratic term
    sub_SIGMA = [zeros(n_states),zeros(n_states,n_inputs);...
        zeros(n_inputs,n_states),eye(n_inputs)];
    
    SIGMA = blkdiag(kron(eye(N_horizon),sub_SIGMA),zeros(n_states));
    
    Matrices.Q_bar = Matrices.Q_bar + ...
        rRefining.CostFunctionCoeff(4).*SIGMA.'*SIGMA;
    
    % second: the linear term
    L_mat = [repmat([zeros(1,n_states),-2.*rRefining.CostFunctionCoeff(3)...
        .*rRefining.u_max],1,N_horizon),zeros(1,n_states)].';
    
    Matrices.R_bar = Matrices.R_bar + ...
        rRefining.CostFunctionCoeff(4).*L_mat;
    
    % third: the additive (affine) term
    affine_term = affine_term + ...
        rRefining.CostFunctionCoeff(4).*N_horizon.*rRefining.u_max.^2;
end

%% Define the gurobi matrices 
% model objective
%   minimize           c'*x + x'*Q*x + alpha
%      x
Model.Q = sparse(Matrices.Q_bar);
% in my definition, R_bar := c'
% gurobi expects c as an input, therefore I have to invert R_bar here
Model.obj = Matrices.R_bar;
% this corresponds to the affine part of the objective ("alpha" above)
Model.objcon = affine_term;

% ...(change) not quite sure if this is right!
if(rRefining.bHalf)
    Model.Q = 0.5.*Model.Q;
    Model.obj = 1.*Model.obj;
end

% check if the problem is convex; if not, throw an error
% ...(change) 
%     for rRefining.CostFunctionCoeff(2)~=0 there is a problem with this
%     part. have to investigate it in the future!
% if(rRefining.CostFunctionCoeff(2)~=0)
%     if(sum(eig(Model.Q) < -0.000000000001) ~= 0)
%         error('Model NOT Convex!');
%     end
% end

%   subject to         A*x = b,
%                      l <= x <= u,
%                      some xj integral,
%                      some xk must lie within second order cones,
%                      x'*Qc*x + q'*x <= beta,
%                      some xi in SOS constraints.
%
% linear equalities and inequalities
eq = repmat(sprintf('='), 1, size(Matrices.A_eq,1));

% matrix A and rhs b
A = [Matrices.A_eq];
b = [Matrices.b_eq];
sense = [eq];

Model.A = sparse(A);
Model.rhs = b;
Model.sense = sense;
Model.modelsense = 'min';

% lower and upper bound
Model.lb = Matrices.lower_bounds;
Model.ub = Matrices.upper_bounds;

% variables type (C = continuous, B = binary)
%
% (generate the right number of 'C's)
vtype = repmat(repmat(sprintf('C'),1,n_z),1,N_horizon+1); % here 105
Model.vtype = vtype(1:size(A,2)); % need only 104

% www.gurobi.com/documentation/6.0/reference-manual/parameter_descriptions
% (for info about the parameters)
params.BarConvTol = 1e-8;
params.BarHomogeneous = 1;

% solve model
Res = gurobi(Model,params);
% disp(Res)

end