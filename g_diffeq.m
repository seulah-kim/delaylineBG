function dydt = g_diffeq(t,y,tau_syn,g_incoming)
prob_syn = 0.8; % probability of successful synaptic transmission
probSyn_gp = rand(length(t_span),n./r)>prob_syn; 
dydt = -y./tau_syn + g_incoming
