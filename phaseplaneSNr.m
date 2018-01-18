% Monte Carlo
spk_gp=[];
spk_snr=[];
spk_str=[];

%Runs 20s simulation without any stimulus to measure steady-state values. 
[g_gp2snr,Isyn_snr] = BGdelayline_setinit2('I_exc_gp',80,'prob_syn_gp2snr',0.6); % initial conductance of gp to snr synapses 

% I_ext = 6 is the number that matches Vthr - Vrest

%g_gp2snr
%Isyn_snr
h = animatedline;
for l = 1:10
%%Simulation
[Vm_gp,Vm_snr,Vm_str, Igp, Isnr] = BGdelayline('stimCellsPer',100,'I_exc_gp',80,'prob_syn_gp2snr',0.6,'g_gp2snr_i',g_gp2snr,'Isyn_snr_i',Isyn_snr);

spk_gp = [spk_gp; Vm_gp==15];
spk_snr = [spk_snr; Vm_snr==15];
spk_str = [spk_str; Vm_str==15];

% for tt = 1:length(Vm_snr)
% addpoints(h,Vm_snr(tt),Isnr(tt));
% title(sprintf('%d s',tt*dt))
% drawnow
% end
% hold on
end


%%
Inull = (-Vm_snr(0.5/dt:1/dt)-70)/100*1000 + 60;

figure
plot(Vm_snr(0.5/dt:0.7/dt),Isnr(0.5/dt:0.7/dt),'b')
hold on 
plot(Vm_snr(1/dt:1.2/dt),Isnr(1/dt:1.2/dt),'r')


%%
figure
plot(Vm_snr,Inull,'r')
hold on

figure
plot(diff(I_snr),Vm_snr)
h = animatedline;
x= diff(Isnr);
y= Vm_snr(3:end);


for tt = 1:length(y)
addpoints(h,x(tt),y(tt));
title(sprintf('%d s',tt*dt))
drawnow
end
hold on
