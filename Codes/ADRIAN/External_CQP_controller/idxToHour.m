%% ========================================================================
%  function main_CL_simulation
%  by Adrian Martinez Gomez
%  January 2015
%
%  Output: hour of the event
%  Input: reference index for the data
%  ========================================================================
function event_hour = idxToHour(ref_exo_idx);
% sub-optimal, but this maps the experiment number back to the hour of the
% day where the change happens
event_hour = nan;
if(1<=mod(ref_exo_idx,140) && mod(ref_exo_idx,140)<=20)
    event_hour = 12;
elseif(21<=mod(ref_exo_idx,140) && mod(ref_exo_idx,140)<=40)
    event_hour = 16;
elseif(41<=mod(ref_exo_idx,140) && mod(ref_exo_idx,140)<=60)
    event_hour = 18;
elseif(61<=mod(ref_exo_idx,140) && mod(ref_exo_idx,140)<=80)
    event_hour = 20;
elseif(81<=mod(ref_exo_idx,140) && mod(ref_exo_idx,140)<=100)
    event_hour = 22;
elseif(101<=mod(ref_exo_idx,140) && mod(ref_exo_idx,140)<=120)
    event_hour = 4;
elseif(121<=mod(ref_exo_idx,140) && mod(ref_exo_idx,140)<=139 ...
        || mod(ref_exo_idx,140)==0)
    event_hour = 7;
end

if(0)
    
    xx = 1:10:3810;
    [yy,~,~]=compare(getexp(zData_same_2,rRefining.idx),ss_16_pos,opt_comp);
    figure;
    subplot(2,1,1);plot(zData_same_2.u{rRefining.idx}(:,1));grid on;
    legend('External input','Location','Best');
    subplot(2,1,2);hold on;
    plot(xx,zData_same_2.y{rRefining.idx}(:,1));
    plot(yy,'r');
    hold off;grid on;
    legend('Real system','ARMAX','Location','Best');
    
end

end