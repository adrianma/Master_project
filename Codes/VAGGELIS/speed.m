function [y, a_ges] = speed(T,U,rho_ref,m_i,d_t,b,g)  
%% ===================================================================
%  input : Temperature array T, velocity array of previous timestep U,
%          timegrid d_t,average density array rho_ref, friction parameter b
%
%  output: new speed of cell
%
%  purpose: calculates the speed of a cell 
%
%  Evangelos Vrettos
%  October 2014
%
%  Adrian Martinez Gomez
%  December 2014
%  ===================================================================

% d_t = Params.t_sample;
% b = Params.b;
% m_i = Params.m_i(:,:,idx);

% Parameter alpha is used to control whether the acceleration is computed
% using the viscous force or not.
if d_t<1
    alpha = 1;
elseif (1<=d_t) && (d_t<=10)
    alpha = (-1/9)*d_t + (10/9);
elseif d_t>10
    alpha = 0;
end

a_ges = (buoy(T,rho_ref,m_i,g) - alpha*b*U)./m_i;
y = U+d_t*a_ges;

end
