% Monte Carlo
spk_gp=[];
spk_snr=[];
spk_str=[];


%Runs 20s simulation without any stimulus to measure steady-state values. 
[g_gp2snr,Isyn_snr] = BGdelayline_setinit2('I_exc_gp',50,'prob_syn_gp2snr',0.35); % initial conductance of gp to snr synapses 

% I_ext = 6 is the number that matches Vthr - Vrest

g_gp2snr
Isyn_snr

for l = 1:10
%%Simulation
[Vm_gp,Vm_snr,Vm_str, Igp, Isnr] = BGdelayline('stimCellsPer',60,'I_exc_gp',50,'prob_syn_gp2snr',0.35,'g_gp2snr_i',g_gp2snr,'Isyn_snr_i',Isyn_snr);

spk_gp = [spk_gp; Vm_gp==15];
spk_snr = [spk_snr; Vm_snr==15];
spk_str = [spk_str; Vm_str==15];
%Isyn_snr = [Isyn_snr ; Isnr];
end

%%Plot
figure
subplot(5,1,1)
plotRaster(spk_str);
title('Spike rater plots (5pA)')
ylabel('Str cells')
subplot(5,1,2)

plotRaster(spk_gp);
ylabel('GP cells')
subplot(5,1,4)

plotRaster(spk_snr);
ylabel('SNr cell')


binWidth_gp = 200; %ms

t_bar_gp = 1:binWidth_gp:length(spk_gp);
psth_gp = zeros(1,length(t_bar_gp));

for psth_i = 1:length(t_bar_gp)-1

    psth_gp(psth_i) = sum(sum(spk_gp(:,t_bar_gp(psth_i):(t_bar_gp(psth_i)+binWidth_gp-1))));

end

%fr_gp = mean(psth_gp(find(1<=t_bar_gp*0.0001<=3)))/(size(spk_gp,1)*(binWidth_gp*0.0001));
subplot(5,1,3)
bar(t_bar_gp*0.0001,psth_gp/(size(spk_gp,1)*(binWidth_gp*0.0001)))
xlim([0 length(spk_gp)*0.0001])
ylabel('gp PSTH (spikes/s)')


binWidth_snr = 200; %ms
t_bar_snr = 1:binWidth_snr:length(spk_snr);
psth_snr = zeros(1,length(t_bar_snr));

for psth_j = 1:length(t_bar_snr)-1
psth_snr(psth_j) = sum(sum(spk_snr(:,t_bar_snr(psth_j):(t_bar_snr(psth_j)+binWidth_snr-1))));
end

subplot(5,1,5)
bar(t_bar_snr*0.0001,psth_snr/(size(spk_snr,1)*(binWidth_snr*0.0001)))
xlim([0 length(spk_gp)*0.0001])
ylabel('snr PSTH (spikes/s)')

