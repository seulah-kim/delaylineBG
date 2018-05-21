clear all

dt=0.0001; % 0.1ms integration steps

%Runs 5s simulation without external current measure the baseline activity.  
[g_gp2snr] = BGdelayline_setinit2('I_exc_gp',60); % initial conductance of gp to snr synapses 

parfor l = 1:100
%Silent striatum, testing different constant excitatory input to GPe
[Vm_gp,Vm_snr,Vm_str, Igp, Isnr] = BGdelayline('n',100,'stimCellsPer',0,'I_exc_gp',60,'I_exc_snr',40,...
'prob_syn_gp2snr',0.35,'g_gp2snr_i',g_gp2snr,'connectivity','all');

spk_gp{l} = Vm_gp==15; %binary spike array
spk_snr{l}=  Vm_snr==15;
spk_str{l}= Vm_str==15;
end
%%


%%

A=cell2mat(spk_snr);