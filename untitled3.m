StimStr = 10:10:100; % phasic activation of striatum amount

for m = 1:length(StimStr)

for sim_i =1:100

%spk_gp=[];
spk_snr=[];
%spk_str=[];
%Isnr_net=[];
dt=0.0001; % 0.1ms integration steps

%Runs 5s simulation without external current measure the baseline activity.  
[g_gp2snr] = BGdelayline_setinit2('I_exc_gp',60,'prob_syn_gp2snr',0.35); % initial conductance of gp to snr synapses 

for l = 1:100
%Silent striatum, testing different constant excitatory input to GPe
[~,Vm_snr,~, ~, ~] = BGdelayline('stimCellsPer',StimStr(m),'I_exc_gp',60,'prob_syn_gp2snr',0.35,'g_gp2snr_i',g_gp2snr);

%meanIsnr = mean(mean(Isnr(:,end-1/dt:end))); %average net inhibition exerted by GPe to SNr
%Isnr_net = [Isnr_net, meanIsnr]; %makes an array of average inhibition across trials

spk_snr = [spk_snr; Vm_snr==15];

end

binWidth_snr = 100; %ms
t_bar = 1:binWidth_snr:length(spk_snr);
psth_snr = zeros(1,length(t_bar));

for psth_j = 1:length(t_bar)-1
psth_snr(psth_j) = sum(sum(spk_snr(:,t_bar(psth_j):(t_bar(psth_j)+binWidth-1))));
end

binloc = (binWidth/2+t_bar-1)*0.0001;
        
[M,I]=max(psth_snr);

if binloc(I)>1 & binloc(I)<1.2
    tdelay(sim_i,i) = binloc(I);
else
    tdelay(sim_i,i) = NaN;
end

end
end
tdelay

figure
errorbar(StimStr,nanmean(tdelay,1),nanstd(tdelay,1)./sqrt(sum(~isnan(tdelay),1)));
xlabel('percent of Str')
ylabel('Time to peak SNr firing (ms)' )
ylim([1 1.1])
yticks(linspace(1,1.1,10))
yticklabels({'0','10','20','30','40','50','60','80','90','100'})
title('response time when initially silent Str is phasically activated ')
