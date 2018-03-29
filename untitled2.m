clear all
StimCell = 70;
for percent_i = 1:length(StimCell)
spk_gp=[];
spk_snr=[];
spk_str=[];
%Isnr_net=[];
dt=0.0001; % 0.1ms integration steps

%Runs 5s simulation without external current measure the baseline activity.  
[g_gp2snr] = BGdelayline_setinit2('I_exc_gp',60,'prob_syn_gp2snr',0.35); % initial conductance of gp to snr synapses 

    for l = 1:100
    %Silent striatum, testing different constant excitatory input to GPe
    [Vm_gp,Vm_snr,Vm_str, ~, ~] = BGdelayline('n',100,'stimCellsPer',StimCell(percent_i),'I_exc_gp',60,'I_exc_snr',40,...
    'prob_syn_gp2snr',0.35,'g_gp2snr_i',g_gp2snr,'connectivity','all');

    %meanIsnr = mean(mean(Isnr(:,end-1/dt:end))); %average net inhibition exerted by GPe to SNr
    %Isnr_net = [Isnr_net, meanIsnr]; %makes an array of average inhibition across trials
    spk_str = [spk_str; Vm_str==15];
    spk_gp = [spk_gp; Vm_gp==15];
    spk_snr = [spk_snr; Vm_snr==15];

    end
    
    binWidth = 100; %ms
    t_bar = 1:binWidth:length(spk_snr);
    psth_snr = zeros(1,length(t_bar));
    psth_gp = zeros(1,length(t_bar));

    for psth_j = 1:length(t_bar)-1
    psth_snr(psth_j) = sum(sum(spk_snr(:,t_bar(psth_j):(t_bar(psth_j)+binWidth-1))));
    psth_gp(psth_j) = sum(sum(spk_gp(:,t_bar(psth_j):(t_bar(psth_j)+binWidth-1))));
    end

    X = psth_snr/(size(spk_snr,1)*(binWidth*0.0001));
%     figure;plot(X)

    histogram(X,'FaceColor',[0 0 1/percent_i])
    hold on
    if percent_i ==1
        Xmean = mean(X(t_bar*0.0001<1));
        Xmean = mean(X(1:end-1));
        Xstd = std(X(1:end-1),1); % normalized by number of observations N 
        CI = 2.58*Xstd;
    vline(Xmean+CI)
    vline(Xmean-CI)
    end
end
    
    

%%
    binloc = (binWidth/2+t_bar-1)*0.0001;

    [M,I]=max(psth_snr);

    if binloc(I)>1 & binloc(I)<1.2
        tdelay(sim_i,m) = binloc(I);
    else
        tdelay(sim_i,m) = NaN;
    end
    %%
    figure; subplot(3,1,1); plotRaster(spk_str)
    subplot(3,1,2);plotRaster(spk_gp)
    subplot(3,1,3);plotRaster(spk_snr)
    