function plotRaster(SpikeM)

%Converts to a binary matrix indicating spike occurence
%SpikeM = traceM==15;
%x and y coordinates
[neuron,time] = find(SpikeM);

%To avoid indexing problem of "find" function for an array
if size(SpikeM,1)==1
    time = find(SpikeM)';
    neuron = ones(length(time),1);
end

neuron = neuron';
time = time';
lineLength = 1/2;
%Creating x and y arrays for vertical lines
xPoints = [ time ;
            time ;
            NaN(size(time)) ];
yPoints = [ neuron - lineLength ;
            neuron + lineLength ;
            NaN(size(neuron)) ];
xPoints = xPoints(:);
yPoints = yPoints(:);
t_plot = xPoints*0.001;
plot(t_plot,yPoints,'k');
      
set(gca,'YDir','reverse');
set(gca,'TickDir','out') 
xlabel('Time (s)');
ylabel('Neuron');
xlim([0 (size(SpikeM,2)+1)*0.001]);
ylim([0 size(SpikeM,1)+1]);
if size(SpikeM,1) == 1
    set(gca,'YTick', [0 1])                        % don't draw y-axis ticks
    set(gca,'PlotBoxAspectRatio',[1 0.05 1])    % short and wide
    ylim([0.5 1.5])
end
set(gca, 'box', 'off')      % Remove top and right axis lines
end
