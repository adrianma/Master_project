function y = distance(Uminew,Umi,d_t,d,L,a_ges)
%% ===================================================================
%  input : speed at time t Umi, speed at time tT1 Uminew, time grid constant d_t
%
%  output: travelled distance in time d_t
%
%  Purpose: calculates the travelled distance of water
% 
%  Evangelos Vrettos
%  October 2014
%
%  Adrian Martinez Gomez
%  December 2014
%  ===================================================================

% d_t = Params.t_sample;
% d = Params.d(:,:,idx);
% L = Params.h(:,:,idx);

if d_t < 10
    y = min(L, max((0.5*(Umi+Uminew)*d_t),0)) / d;
else
    y = min(L, max(0.5*a_ges*(d_t^2),0)) / d;
end

end
