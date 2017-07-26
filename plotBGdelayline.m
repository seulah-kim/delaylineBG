%% plotBGdelayline.m
% This code shall be run following BGdelayline{2,3}.m to visualize neural activitiy in each layer of the delay line 
figure
subplot(2,1,1)
plot(Isyn_gp_out')
subplot(2,1,2)
plot(Isyn_snr_out')

figure
subplot(3,1,1)
 plot(t_span,Vm_str(5,1:length(t_span)),'b')
hold on
%  plot(t_span,Vm_gp(1,1:length(t_span))','b')
 plot(t_span,Vm_snr(1,1:length(t_span)),'r')
% legend('Vm_{str}','Vm_{gp}','Vm_{snr}')
legend('Vm_{str} #1','Vm_{snr}')
ylabel('membrane Voltage (mV)')
% xlabel('Time (s)')

subplot(3,1,2) 
hold on
plot(t_span,g_gp2snr(:,1:length(t_span)),'g')
% legend('g_{str2gp}','g_{gp2snr}')
legend('g_{str2gp}')
ylabel('conductance (nS)')

subplot(3,1,3)
plot(t_span,Vm_str(5,1:length(t_span)),'g')
hold on
plot(t_span,Vm_gp(:,1:length(t_span)),'b')
%  plot(t_span,Vm_snr(1,1:length(t_span)),'r')
% legend('Vm_{str}','Vm_{gp}','Vm_{snr}')
legend('Vm_{str} #1','Vm_{gp} #1')
ylabel('membrane Voltage (mV)')
xlabel('Time (s)')

%%total synaptic conductance
% g_str=sum(spk_count_str.*probSyn_str      % how to advance through probSyn_str??i
