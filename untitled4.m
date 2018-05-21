clear all

dt=0.0001; % 0.1ms integration steps

%Runs 5s simulation without external current measure the baseline activity.  
[g_gp2snr] = BGdelayline_setinit2('I_exc_gp',60); % initial conductance of gp to snr synapses 

for i = 1:2
    switch i
        case 1
            exp = 0;
        case 2
            exp = 100;
    end
        
    clear spk_str spk_snr spk_gp
    
    parfor l = 1:100
    %Silent striatum, testing different constant excitatory input to GPe
    [Vm_gp,Vm_snr,Vm_str, Igp, Isnr] = BGdelayline('n',100,'stimCellsPer',exp,'I_exc_gp',60,'I_exc_snr',40,...
    'prob_syn_gp2snr',0.35,'g_gp2snr_i',g_gp2snr,'connectivity','all');

    spk_gp{l,1} = Vm_gp==15; %binary spike array
    spk_snr{l,1}=  Vm_snr==15;
    spk_str{l,1}= Vm_str==15;
    end

    spk_str = double(cell2mat(spk_str));
    spk_snr = double(cell2mat(spk_snr));
    spk_gp = double(cell2mat(spk_gp));

    binWidth = 50; %0.1 ms

t_bar = 1:binWidth:length(spk_snr);
%psth_gp = zeros(1,length(t_bar));
psth_snr = zeros(1,length(t_bar));

for psth_i = 1:length(t_bar)-1

    %psth_gp(psth_i) = sum(sum(spk_gp(:,t_bar(psth_i):(t_bar(psth_i)+binWidth-1))));
    psth_snr(psth_i) = sum(sum(spk_snr(:,t_bar(psth_i):(t_bar(psth_i)+binWidth-1))));
end

%figure;
X=psth_snr(1:end-1)/(size(spk_snr,1)*(binWidth*0.0001));
%plot((binWidth/2+t_bar(1:end-1)-1)*0.0001,X)
%hold on;
%plot((binWidth/2+t_bar(1:end-1)-1)*0.0001,filter(1/10, [1 1/10-1],X))

L = length(X);
Fs = 1/(binWidth*0.0001);% sampling rate

y = fft(X,1024);
P = y.*conj(y)/(1024*L);
length(y)

figure(1);
f = (0:1024/2-1)*Fs/1024;        % Frequency vector
length(f)
plot(log(f),P(1:1024/2))
xlabel('Freq (Hz)')
xticks(log([0.1 1 10 100 1000]))
xticklabels({'0.1','1','10','100','1000'})

hold on;

%figure;
%histogram(X(1:end-1),'Normalization','probability')

avg=mean(X(1:end-1))
stdev=std(X(1:end-1))


%figure;
%plotRaster(spk_snr)
end
