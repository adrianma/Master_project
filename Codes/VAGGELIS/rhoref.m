function rho_ref = rhoref(T)
%% ===================================================================
%  input : Tempature array T, cell above which rho_ref should be
%  calculated j, total number of cells n
%
%  output: average density above cell j
%
%  Evangelos Vrettos
%  October 2014
%  ===================================================================

a = size(T,1);
rho_ref = zeros(size(T));
dens = density(T);

if(0)
    % for each layer j, rho_ref(j) is the maximum between rho_ref(j+1) and mean(rho_ref(j+1:end))
    for i = 1:a-1
        rho_ref(i)= max(mean(dens(i+1:a)),dens(i+1));
    end
    rho_ref(a) = dens(a);
end

%% this unwraps the loop from above to accelerate this function x6
% unfortunately, this takes out the ability for variable number of layers
% (here only 10)
if(0)
    rho_ref(1) = max(sum(dens(2:a))/(9),dens(2));
    rho_ref(2) = max(sum(dens(3:a))/(8),dens(3));
    rho_ref(3) = max(sum(dens(4:a))/(7),dens(4));
    rho_ref(4) = max(sum(dens(5:a))/(6),dens(5));
    rho_ref(5) = max(sum(dens(6:a))/(5),dens(6));
    rho_ref(6) = max(sum(dens(7:a))/(4),dens(7));
    rho_ref(7) = max(sum(dens(8:a))/(3),dens(8));
    rho_ref(8) = max(sum(dens(9:a))/(2),dens(9));
    rho_ref(9) = max(sum(dens(a:a))/(1),dens(10));
    rho_ref(a) = dens(a);
end

if(1)
    vv = flip(cumsum(flip(dens))./(1:1:a)');
    rho_ref = [max(vv(2:end),dens(2:a));dens(a)];
end


end

    
