 %% Condition 1 - 3pA
tic
spk_gp=[];
spk_snr=[];
spk_str=[];
Isyn_snr=[];

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
        for l = 1:10
        %%Simulation
        [Vm_gp,Vm_snr,Vm_str, Isyn_gp, Isnr] = BGdelayline('stimCellsPer',100,'I_exc_gp',gp_in(param1),'prob_syn_gp2snr',P_gp2snr(param2));

        spk_gp = [spk_gp; Vm_gp==15];
        spk_snr = [spk_snr; Vm_snr==15];
        spk_str = [spk_str; Vm_str==15];
        Isyn_snr = [Isyn_snr ; Isnr];
        end

        binWidth_gp = 200; %ms

        t_bar_gp = 1:binWidth_gp:length(spk_gp);
        psth_gp = zeros(1,length(t_bar_gp));

        for psth_i = 1:length(t_bar_gp)-1

            psth_gp(psth_i) = sum(sum(spk_gp(:,t_bar_gp(psth_i):(t_bar_gp(psth_i)+binWidth_gp-1))));

        end

        fr_gp(param1,param2) = mean(psth_gp(find(1<=t_bar_gp*0.001<=3)))/(size(spk_gp,1)*(binWidth_gp*0.001));
        I_tot_snr_mean(param1,param2) = mean(mean(Isyn_snr(:,2500:4000)));
        I_tot_snr_std(param1,param2) = std(mean(Isyn_snr(:,2500:4000)));
     end
    
end
toc
% figure;errorbar(P_gp2snr,I_tot_snr_mean(1,:),I_tot_snr_std(1,:))
% xlabel('firing rate (Hz)')
% ylabel('GP-to-SNr net inhibition per cell')

%Surface plot
figure;surfc(fr_gp(:,1),P_gp2snr,I_tot_snr_mean')

%%
%%Plot
figure (1)
subplot(5,1,1)
plotRaster(spk_str);
title('Spike rater plots (3pA)')
ylabel('Str cells')
subplot(5,1,2)

plotRaster(spk_gp);
ylabel('GP cells')
subplot(5,1,4)

plotRaster(spk_snr);
ylabel('SNr cell')


figure (2)
subplot(2,1,1)
plot(Isyn_gp')
legend('gp_inhibition')
title('input to GP = 3pA')
subplot(2,1,2)
plot(Isyn_snr')
legend('snr_inhibition')

%%

binWidth_gp = 200; %ms

t_bar_gp = 1:binWidth_gp:length(spk_gp);
psth_gp = zeros(1,length(t_bar_gp));

for psth_i = 1:length(t_bar_gp)-1
   
    psth_gp(psth_i) = sum(sum(spk_gp(:,t_bar_gp(psth_i):(t_bar_gp(psth_i)+binWidth_gp-1))));
   
end

fr_gp = mean(psth_gp(find(1<=t_bar_gp*0.001<=3)))/(size(spk_gp,1)*(binWidth_gp*0.001));

figure (1)
subplot(5,1,3)
bar(t_bar_gp*0.001,psth_gp/(size(spk_gp,1)*(binWidth_gp*0.001)))
xlim([0 7])
ylabel('gp PSTH (spikes/s)')

binWidth_snr = 200; %ms
t_bar_snr = 1:binWidth_snr:length(spk_snr);
psth_snr = zeros(1,length(t_bar_snr));

for psth_j = 1:length(t_bar_snr)-1
psth_snr(psth_j) = sum(sum(spk_snr(:,t_bar_snr(psth_j):(t_bar_snr(psth_j)+binWidth_snr-1))));
end

subplot(5,1,5)
bar(t_bar_snr*0.001,psth_snr/(size(spk_snr,1)*(binWidth_snr*0.001)))
xlim([0 7])
ylabel('snr PSTH (spikes/s)')

% [N, EDGES] = histcounts(spk_gp(1,:),100);
