%% BG delay line model
close all
clear 
clc
r = 10; % convergence ratio between layers
n = 100; % number of str neurons

dt = 0.001; %(s)
t_span = 0:dt:20;

% Define Constants
%%cellular
tau_cell_str = ones(n,1);   % cell decay constant (ms)
tau_cell_gp = ones(n/r,1);   % cell decay constant (ms)
tau_cell_snr = ones(n/r.^2,1);   % cell decay constant (ms)
Cm = 0.1;        % (nF)
Erev_i = -85;   % (mV)
I_const = 0.1;   % (pA)
Vrest = -70;    % (mV)
Vpeak = 15;     % (mV)
V_thres = -64;  % (mV)
k =1000;    % steepness of logistic function
xd = 5;    % x location of midpoint of logistic function
L = 40;     % max value of curve
synNonlin = @(x)[L./(1+exp(-k*(x-xd)))];   % check if this is in the working range

%%synaptic
tau_syn =5; % synaptic decay constant (ms)
prob_syn = 0.9; % probability of successful synaptic transmission
g_uni = 1 ; % conductance of a single synapse
% synSuccess_i = rand(length(Vm_i),1)<prob_syn; % flip a coin: success (1) or failure (0)
% probSyn_gp = rand(n./r)>prob_syn; 

% probSyn_gp = rand(length(t_span),n./r)>prob_syn; 
% sum(spk_count_str.*synSuccess

g_str2gp = 0.3*ones(n/r,1);    % (nS)
g_gp2snr = 0.3*ones(n/r.^2,1);
Vm_str = Vrest+0.2*randn(n,1);
Vm_gp = Vrest+0.2*randn(n/r,1);
Vm_snr = Vrest+0.2*randn(n/r.^2,1);

% Jump Vm_i if it exceeds threshold voltage?
delta_str = zeros(n,1);
delta_gp = zeros(n/r,1);
delta_snr = zeros(n/r.^2,1);

Isyn_gp_out =[];
Isyn_snr_out =[];
for t = 1:length(t_span)

%Striatum
% delta_str = double(Vm_str>V_thres);
% Vm_str(Vm_str>V_thres) = Vrest;   % bring it back to Vrest
Iext_str = 7*I_const*ones(n,1);  % constant input to striatum ( 0.6~0.61Iconst for relatively quiet str cells) (in pA)
dVm_str = (-1*(Vm_str(:,t)-Vrest*ones(n,1))./tau_cell_str + Iext_str/Cm)*dt;   % no nonlinear function here!

%GP
% delta_gp = double(Vm_gp>V_thres);
% Vm_gp(Vm_gp>V_thres) = Vrest;   % bring it back to Vrest
synSuccess_str2gp = double(rand(n,n/r)<prob_syn);
Isyn_gp = g_str2gp(:,t).*(Vm_gp(:,t)-Erev_i*ones(n/r,1));  % (pA)
Iext_gp = 50*I_const*ones(n/r,1);
dg_str2gp = (-g_str2gp(:,t)./tau_syn + transpose(delta_str'*synSuccess_str2gp*g_uni))*dt;
% dVm_gp = (-(Vm_gp(:,t)-Vrest*ones(n/r,1))./tau_cell_gp - synNonlin(Isyn_gp)/Cm+Iext_gp/Cm)*dt;
dVm_gp = (-(Vm_gp(:,t)-Vrest*ones(n/r,1))./tau_cell_gp +Iext_gp/Cm)*dt;
%SNr
% Vm_snr(Vm_snr>Vm_thres) = Vrest;   % bring it back to Vrest
synSuccess_gp2snr = double(rand(n/r,n/r.^2)<prob_syn);
Isyn_snr = -g_gp2snr(:,t).*(Vm_snr(:,t)-Erev_i.*ones(n/r.^2,1));
Iext_snr = Isyn_snr+60*I_const*ones(n/r.^2,1);
dg_gp2snr = (-g_gp2snr(:,t)./tau_syn + transpose(delta_gp'*synSuccess_gp2snr*g_uni))*dt;
dVm_snr = (-(Vm_snr(:,t)-Vrest*ones(n/r.^2,1))./tau_cell_snr + synNonlin(Iext_snr)/Cm)*dt;

%%Update
g_str2gp(:,t+1) = g_str2gp(:,t)+dg_str2gp;
g_gp2snr(:,t+1) = g_gp2snr(:,t)+dg_gp2snr;

% check for spikes in the cells
%update spike patterns
delta_str = double((Vm_str(:,t)+dVm_str)>=V_thres);
delta_gp = double((Vm_gp(:,t)+dVm_gp)>=V_thres);
delta_snr = double((Vm_snr(:,t)+dVm_snr)>=V_thres);

Vm_str(:,t+1) = Vm_str(:,t)+dVm_str;
if sum(delta_str)>=1
    Vm_str((Vm_str(:,t)+dVm_str)>=V_thres,t) = Vpeak;
    Vm_str((Vm_str(:,t)+dVm_str)>=V_thres,t+1) = Vrest;
end

Vm_gp(:,t+1) = Vm_gp(:,t)+dVm_gp;
if sum(delta_gp)>=1
    Vm_gp((Vm_gp(:,t)+dVm_gp)>=V_thres,t) = Vpeak;
    Vm_gp((Vm_gp(:,t)+dVm_gp)>=V_thres,t+1) = Vrest;
end

Vm_snr(:,t+1) = Vm_snr(:,t)+dVm_snr;
if sum(delta_snr)>=1
    Vm_snr((Vm_snr(:,t)+dVm_snr)>=V_thres,t) = Vpeak;
    Vm_snr((Vm_snr(:,t)+dVm_snr)>=V_thres,t+1) = Vrest;
end 

%save variables 
Isyn_gp_out = [Isyn_gp_out,Isyn_gp];
Isyn_snr_out = [Isyn_snr_out,Isyn_snr];
end
%%
figure
subplot(2,1,1)
plot(Isyn_gp_out')
subplot(2,1,2)
plot(Isyn_snr_out')
%%
figure
subplot(2,1,1)
 plot(t_span,g_str2gp(:,1:length(t_span)),'r')
hold on
% plot(t_span,g_gp2snr(:,1:length(t_span)),'g')
legend('g_{str2gp}','g_{gp2snr}')
ylabel('conductance (nS?)')

subplot(2,1,2)
%  plot(t_span,Vm_str(1,1:length(t_span))','k')
hold on
 plot(t_span,Vm_gp(1,1:length(t_span))','b')
%  plot(t_span,Vm_snr(:,1:length(t_span)),'r')
legend('Vm_{str}','Vm_{gp}','Vm_{snr}')
ylabel('membrane Voltage (mV)')
xlabel('Time (s)')
%%total synaptic conductance
% g_str=sum(spk_count_str.*probSyn_str      % how to advance through probSyn_str??