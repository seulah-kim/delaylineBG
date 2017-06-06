function dydt = Vm_diffeq(t,y,tau_cell,C,g)

I_syn =  g*(y-E_gpe);
dydt = -1./tau*y+