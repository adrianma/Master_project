function y = mixing(T,dis,rr,n)  
%% =================================================================== 
%  input: temperature of cells array T, distance array dis, number of cells
%         n, mixing ratio r
%
%  output:  new temperature array
%
%  purpose: calculates the temperature array of next timestep
% 
%  Evangelos Vrettos
%  October 2014
%
%  Adrian Martinez Gomez
%  December 2014
%  ===================================================================

% rr = Params.r;
% n = Params.n;

for j = 1:n
    i = 1;
    r = rr;
    while i < dis(j)+1 && j+i <= n
        if (dis(j) < i) 
            r = r*(dis(j)-floor(dis(j))); 
        end
        d = T(j+i-1);
        T(j+i-1) = r*T(j+i) + (1-r)*T(j+i-1);
        T(j+i) = r*d + (1-r)*T(j+i);
        i = i + 1;
    end
end

y = T;
end