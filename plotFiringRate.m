function [r_str, r_gp, r_snr]= plotFiringRate(traceM)

%Converts to a binary matrix indicating spike occurence
SpikeM = traceM==15;

binWidth = 200; % in ms
binSpeed = 5; % in ms/step
binEnd = length(traceM);
tt = 1:binSpeed:(binEnd-binWidth+1);
spikeCount = [];

for  i_bin = 1:length(tt)
	spikeCount = [spikeCount round(sum(SpikeM(:,tt(i_bin):(tt(i_bin)+binWidth)),2)/(0.001*binWidth),2)];
end
plot(linspace(0,20,length(tt)),spikeCount')
end
