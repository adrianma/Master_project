function y = buoy(T,rho_ref,m_i,g)  
%% ===================================================================  
%  input : temperature array T,average density above cell array rho_ref, volume of cell v_0
%  
%  output: buoyant force on cell 
%  
%  Purpose: calculate the buoyancy force according to Vrettos-Koch-Andersson
%
%  Evangelos Vrettos
%  October 2014
%
%  Adrian Martinez Gomez
%  December 2014
%  ===================================================================

% m_i = Params.m_i(:,:,idx);
% g = Params.g;
rho_i = density(T);

y = g*m_i.*(rho_ref-rho_i)./(rho_ref+rho_i);
