% excitatory input to GP parameter range
gp_in = 80:20:270;    %80 to 270 correspond to 10 to 40Hz
fr_gp = zeros(length(gp_in),1);

% gp2snr synaptic transmission range
P_gp2snr = 0.1:0.05:0.5;
I_tot_snr_mean = zeros(length(gp_in),length(P_gp2snr));
I_tot_snr_std = zeros(length(gp_in),length(P_gp2snr));
for param1 = 1:length(gp_in)
    for param2 = 1:length(P_gp2snr)
        spk_gp=[];
        spk_snr=[];
        spk_str=[];
        Isyn_snr=[];
        %Runs 20s simulation without any stimulus to measure steady-state values. 
        g_gp2snr = BGdelayline_setinit('I_exc_gp',80,'prob_syn_gp2snr',0.35); % initial conductance of gp to snr synapses 
        
        for l = 1:10
        %%Simulation
        [Vm_gp,Vm_snr,Vm_str, Isyn_gp, Isnr] = BGdelayline('stimCellsPer',100,'I_exc_gp',gp_in(param1),'prob_syn_gp2snr',P_gp2snr(param2),'g_gp2snr_i',g_gp2snr);

        spk_gp = [spk_gp; Vm_gp==15];
        spk_snr = [spk_snr; Vm_snr==15];
        spk_str = [spk_str; Vm_str==15];
        Isyn_snr = [Isyn_snr ; Isnr];
        end

        binWidth_gp = 50; %ms

        t_bar_gp = 1:binWidth_gp:length(spk_gp);
        psth_gp = zeros(1,length(t_bar_gp));

        for psth_i = 1:length(t_bar_gp)-1

            psth_gp(psth_i) = sum(sum(spk_gp(:,t_bar_gp(psth_i):(t_bar_gp(psth_i)+binWidth_gp-1))));

        end

        fr_gp(param1,param2) = mean(psth_gp(find(1<=t_bar_gp*0.001<=3)))/(size(spk_gp,1)*(binWidth_gp*0.001));
        I_tot_snr_mean(param1,param2) = mean(mean(Isyn_snr(:,3000:5000)));
        I_tot_snr_std(param1,param2) = std(mean(Isyn_snr(:,3000:5000)));
     end
    
end