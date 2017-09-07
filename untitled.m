testRange = 10*(6:10);

for test_i = 1:length(testRange)
    spk_gp=[];
    spk_snr=[];
    spk_str=[];
    for l = 1:1
    %%Simulation
        [Vm_gp,Vm_snr,Vm_str] = BGdelayline('stimCellsPer',testRange(test_i));

        spk_gp = [spk_gp; Vm_gp==15];
        spk_snr = [spk_snr; Vm_snr==15];
        spk_str = [spk_str; Vm_str==15];
    end

    %%Plot
    figure
    subplot(3,1,1)
    plotRaster(spk_str);
    title(sprintf('Spike rater plots (5pA),%d percent',testRange(test_i)))
    ylabel('Str cells')
    subplot(3,1,2)

    plotRaster(spk_gp);
    ylabel('GP cells')
    subplot(3,1,3)

    plotRaster(spk_snr);
    ylabel('SNr cell')
    
end