function [xnew, Uminew] = buoyCor(x,Params,Umi,idx)
%% ===================================================================
%  This function takes as input the temperature distribution within the 
%  tank at current step and calculates the new distribution after buoyancy
%  is considered. Buoyancy is calculated based on the density gradient 
%  of water layers and Newton's second law of motion. 
%
%  IMPORTANT: different methods are used depending on Params.t_sample
%  Method 1: valid only for small time steps (< 10 sec).
%            the displacement of layers depends only on the buoyant force
%            and the viscous force, i.e., we need to keep track of the
%            speeds of different layers
%            
%  Method 2: valid for large time steps (> 10 sec).
%            the displacement of layers depends only on the buoyant force;
%            the viscous force is zero because the time step is large
%            enough to assume that the speed of the layers is zero at the
%            end of each time step (i.e., a speed equilibrium has been reached)
%
%  Evangelos Vrettos
%  October 2014
%
%  Adrian Martinez Gomez
%  December 2014
%  ===================================================================

%% Calculate the acceleration of each disk due to buoyancy forces
% Calculate the density difference of each layer j with the layers above layer j
T = x(2:length(x)-1);
rho_ref = rhoref(T);

% Predefinitions to reduce overhead
m_i = Params.m_i(:,:,idx);
d_t = Params.t_sample;
b = Params.b;
d = Params.d(:,:,idx);
L = Params.h(:,:,idx);
g = Params.g;
rr = Params.r;
n = Params.n;

% Calculate the velocity of each layer
[Uminew, a_ges] = speed(T,Umi,rho_ref,m_i,d_t,b,g);

% Calculate displacement for each layer
dis = distance(Uminew,Umi,d_t,d,L,a_ges);

% Calculate new temperatures after mixing
%
% EDIT: CALL THE MEX FILE!
% note that you need some inversions
if(0)
    xtmp = mixing_mex(T.',dis.',rr,int8(n));
    xnew = [x(1); xtmp.'; x(end)];
end

if(1)
    xtmp = mixing(T,dis,rr,n);
    xnew = [x(1); xtmp; x(end)];
end


end

%% Recycle Bin