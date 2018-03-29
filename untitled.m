GPe_input = [80]; % pA, inhibitory
for i = 1:length(GPe_input)
%Runs 20s simulation without any stimulus to measure steady-state values. 
[g_gp2snr] = BGdelayline_setinit2('I_exc_gp',GPe_input(i),'prob_syn_gp2snr',0.35); % initial conductance of gp to snr synapses 

spk_gp = [];
spk_snr=[];
dt=0.0001; % 0.1ms integration steps

[Vm_gp, Vm_snr, Isyn_out] = BGdelayline_GPe2SNr('Iexc_snr',80,'Iexc_gp',GPe_input(i),'g_gp2snr_i',g_gp2snr,'n',100); % only SNr layer. 100 neurons
spk_snr = [spk_snr; Vm_snr==15];
spk_gp = [spk_gp ; Vm_gp ==15];

binWidth = 100; % bin size, scale of 0.1ms

t_bar = 1:binWidth:length(spk_snr);
psth_snr = zeros(1,length(t_bar));
psth_gp = zeros(1,length(t_bar));

for psth_i = 1:length(t_bar)-1

    psth_snr(psth_i) = sum(sum(spk_snr(:,t_bar(psth_i):(t_bar(psth_i)+binWidth-1))));
    psth_gp(psth_i) = sum(sum(spk_gp(:,t_bar(psth_i):(t_bar(psth_i)+binWidth-1))));
    
end


figure
subplot(3,1,1)
bar((binWidth/2+t_bar-1)*0.0001,psth_gp/(size(spk_gp,1)*(binWidth*0.0001)))
xlim([0.95 1.05])
vline(1)
ylabel('gp PSTH (spikes/s)')
ylim([0 100])
xlim([0.9 1.2])
xticks([0.9 1 1.1 1.2])
xticklabels({'-100','0','100','200'})
title(sprintf('input to GPe: %d pA, net excitation to SNr = 100pA',GPe_input(i)))

subplot(3,1,2)
plotRaster(spk_snr);
vline(1)
ylabel('SNr cell')
xlim([0.9 1.2])
xticks([0.9 1 1.1 1.2])
xticklabels({'-100','0','100','200'})


subplot(3,1,3)
bar((binWidth/2+t_bar-1)*0.0001,psth_snr/(size(spk_snr,1)*(binWidth*0.0001)))
xlim([0.9 1.2])
xticks([0.9 1 1.1 1.2])
xticklabels({'-100','0','100','200'})
ylim([0 100])
vline(1)
ylabel('snr PSTH (spikes/s)')
end
