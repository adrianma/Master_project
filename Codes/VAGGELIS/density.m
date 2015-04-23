% source: McCutcheon, S.C., Martin, J.L, Barnwell, T.O. Jr. 1993. Water 
% Quality in Maidment, D.R. (Editor). Handbood of Hydrology, 
% McGraw-Hill, New York, NY (p. 11.3 )

function p_dens=density(T)

p_dens = 1000*(1 - (T+288.9414)./(508929.2.*(T+68.12963)).*(T-3.9863).^2);
