function [value,isterminal,direction]=spikeevents(t,y)
 global Vrest n r
 %if y>=V_t
    value = y((n/r+n/r.^2+1):end)+64*ones(n+n/r+n/r.^2,1);
% en
isterminal = ones(n+n/r+n/r.^2,1);
direction = ones(n+n/r+n/r.^2,1); 
