spk_gp=[];
spk_snr=[];
spk_str=[];
Isnr_net=[];
dt=0.0001; % 0.1ms integration steps

%Runs 5s simulation without external current measure the baseline activity.  
[g_gp2snr] = BGdelayline_setinit2('I_exc_gp',60); % initial conductance of gp to snr synapses 

for l = 1:100
%Silent striatum, testing different constant excitatory input to GPe
[Vm_gp,Vm_snr,Vm_str, Igp, Isnr] = BGdelayline('n',100,'stimCellsPer',100,'I_exc_gp',70,'I_exc_snr',30,...
'prob_syn_gp2snr',0.35,'g_gp2snr_i',g_gp2snr,'connectivity','all');

%meanIsnr = mean(mean(Isnr(:,end-1/dt:end))); %average net inhibition exerted by GPe to SNr
%Isnr_net = [Isnr_net, meanIsnr]; %makes an array of average inhibition across trials

spk_gp = [spk_gp; Vm_gp==15]; %binary spike array
spk_snr = [spk_snr; Vm_snr==15];
spk_str = [spk_str; Vm_str==15];
end

%%Plot
figure(1)
subplot(3,1,1)
plotRaster(spk_str);
title(sprintf('%d percent str activated',100))
ylabel('Str cells')
vline(1)
xlabel('')
xlim([0.9 1.2])
xticks([0.9 1 1.1 1.2])
xticklabels({'-100','0','100','200'})

subplot(3,1,2)
plotRaster(spk_gp);
ylabel('GP cells')
vline(1)
xlabel('')
xlim([0.9 1.2])
xticks([0.9 1 1.1 1.2])
xticklabels({'-100','0','100','200'})

subplot(3,1,3)
plotRaster(spk_snr);
ylabel('SNr cell')
vline(1)
xlabel('')
xlim([0.9 1.2])
xticks([0.9 1 1.1 1.2])
xticklabels({'-100','0','100','200'})

binWidth = 100; %0.1 ms

t_bar = 1:binWidth:length(spk_snr);
psth_gp = zeros(1,length(t_bar));
psth_snr = zeros(1,length(t_bar));

for psth_i = 1:length(t_bar)-1

    psth_gp(psth_i) = sum(sum(spk_gp(:,t_bar(psth_i):(t_bar(psth_i)+binWidth-1))));
    psth_snr(psth_i) = sum(sum(spk_snr(:,t_bar(psth_i):(t_bar(psth_i)+binWidth-1))));
end

figure(2)
%fr_gp = mean(psth_gp(find(1<=t_bar_gp*0.0001<=3)))/(size(spk_gp,1)*(binWidth_gp*0.0001));
subplot(2,1,1)
bar((binWidth/2+t_bar-1)*0.0001,psth_gp/(size(spk_gp,1)*(binWidth*0.0001)))

ylabel('gp PSTH (spikes/s)')
ylim([0 100])
vline(1)
xlim([0.9 1.2])
xticks([0.9 1 1.1 1.2])
xticklabels({'-100','0','100','200'})
title(sprintf('input to GPe: %d pA, net excitation to SNr = 70pA',60))

subplot(2,1,2)
bar((binWidth/2+t_bar-1)*0.0001,psth_snr/(size(spk_snr,1)*(binWidth*0.0001)))
ylim([0 100])
vline(1)
xlim([0.9 1.2])
xticks([0.9 1 1.1 1.2])
xticklabels({'-100','0','100','200'})
ylabel('snr PSTH (spikes/s)')

[M,I]=max(psth_snr(2:end));
tdelay = t_bar(I)*dt;

figure(3)
subplot(2,1,1)
plot(linspace(0,3,length(Igp)),Igp')
ylabel('pA')
vline(1)
legend('Incoming current to gp')

subplot(2,1,2)
plot(linspace(0,3,length(Igp)),Isnr(1,:))
ylabel('pA')
legend('Incoming current to snr')
vline(1)

