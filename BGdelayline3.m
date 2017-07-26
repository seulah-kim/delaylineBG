%% BGdelayline3.M
%This is the main code for simulating BG delay line model
tic
close all
clear 
clc
%%Network size
r = 10; % convergence ratio between layers
n = 100; % number of str neurons

%%Simulation time
dt = 0.001; %(s)
t_span = 0:dt:20;

%%Define Constants
%Cellular
tau_cell_str = ones(n,1);   % cell decay constant (ms)
tau_cell_gp = ones(n/r,1);   % cell decay constant (ms)
tau_cell_snr = ones(n/r.^2,1);   % cell decay constant (ms)
Cm = 0.1;        % Cell capacitance (nF)
Erev_i = -85;   % Synaptic reversal potential(mV)
I_const = 0.1;   % General scale for constant excitatory input(pA)
Vrest = -70;    % Resting potential(mV)
Vpeak = 15;     % Peak potential(mV)
V_thres = -64;  % Threshold voltage(mV)
%Synaptic
tau_syn =5; % synaptic decay constant (ms)
prob_syn = 0.95; % probability of successful synaptic transmission
g_uni = 1 ; % conductance of a single synapse

%%Nonlinear function to filter synaptic input (sigmoid function)
%values for GP 
k =1000;    % steepness of logistic function
xd = 5;    % x location of midpoint of logistic function
L = 40;     % max value of curve
syms synNonlin(x)
synNonlin(x)= piecewise(x<=(xd-2/k),0,(xd-2/k)<x<=xd, L*k/4*(x-xd)+L/2,x>xd,L./(1+exp(-k*(x-xd))));    % check if this is in the working range
%values for SNr
k =1000;    % steepness of logistic function
xd = 6.5;    % x location of midpoint of logistic function
L = 40;     % max value of curve
%synNonlin_snr = @(x)[L./(1+exp(-1000*(x-xd)))];   % check if this is in the working range
syms synNonlin_snr(x)
synNonlin_snr(x) =piecewise(x<=(xd-2/k),0,(xd-2/k)<x<=xd, L*k/4*(x-xd)+L/2,x>xd,L./(1+exp(-k*(x-xd))));   % check if this is in the working range

%%Input current to Str layer
randNum = round((20-10).*rand(1,1)+10);
randNum = 5;
tStim = [randNum:dt:randNum+0.5]; % (s) 
Iext_str = zeros(n,length(t_span));
% only the Str cell #5 will receive input 
Iext_str(5,ismember(t_span,tStim)) = 50*I_const;  % constant input to striatum ( 0.6~0.61Iconst for relatively quiet str cells) (in pA)

%%Initialize variables
%homogeneous conductance
g_str2gp = 0.05*ones(n/r,1);    % (nS)
g_gp2snr = 0.3*ones(n/r.^2,1);
%heterogeneous Vstart
Vm_str = Vrest+2*randn(n,1);  % (mV)
Vm_gp = Vrest+2*randn(n/r,1);
Vm_snr = Vrest+2*randn(n/r.^2,1);
delta_str = zeros(n,1); % binary
delta_gp = zeros(n/r,1);
delta_snr = zeros(n/r.^2,1);
Isyn_gp_out =[];
Isyn_snr_out =[];

for t = 1:length(t_span)
%Striatum
%Iext_str = I_const*ones(n,1);  % constant input to striatum ( 0.6~0.61Iconst for relatively quiet str cells) (in pA)
dVm_str = (-1*(Vm_str(:,t)-Vrest*ones(n,1))./tau_cell_str + Iext_str(:,t)/Cm)*dt;   % no nonlinear function here!

%GP
synSuccess_str2gp = double(rand(n,n/r)<prob_syn);
Isyn_gp = g_str2gp(:,t).*(Vm_gp(:,t)-Erev_i*ones(n/r,1));  % (pA)
Iext_gp = 50*I_const*ones(n/r,1);
dg_str2gp = (-g_str2gp(:,t)./tau_syn + transpose(delta_str'*synSuccess_str2gp*g_uni))*dt;
dVm_gp = (-(Vm_gp(:,t)-Vrest*ones(n/r,1))./tau_cell_gp - synNonlin(Isyn_gp)/Cm+Iext_gp/Cm)*dt;
% dVm_gp = (-(Vm_gp(:,t)-Vrest*ones(n/r,1))./tau_cell_gp +Iext_gp/Cm)*dt; %

%SNr
synSuccess_gp2snr = double(rand(n/r,n/r.^2)<prob_syn);
Isyn_snr = g_gp2snr(:,t).*(Vm_snr(:,t)-Erev_i.*ones(n/r.^2,1));
Iext_snr = 40*I_const*ones(n/r.^2,1);
dg_gp2snr = (-g_gp2snr(:,t)./tau_syn + transpose(delta_gp'*synSuccess_gp2snr*g_uni))*dt;
dVm_snr = (-(Vm_snr(:,t)-Vrest*ones(n/r.^2,1))./tau_cell_snr - synNonlin_snr(Isyn_snr)/Cm+Iext_snr/Cm)*dt;

%Update conductances
g_str2gp(:,t+1) = g_str2gp(:,t)+dg_str2gp;
g_gp2snr(:,t+1) = g_gp2snr(:,t)+dg_gp2snr;

%Check for spikes in the cells
delta_str = double((Vm_str(:,t)+dVm_str)>=V_thres);
delta_gp = double((Vm_gp(:,t)+dVm_gp)>=V_thres);
delta_snr = double((Vm_snr(:,t)+dVm_snr)>=V_thres);

%Update spike patterns
%Str
Vm_str(:,t+1) = Vm_str(:,t)+dVm_str;
if sum(delta_str)>=1
    Vm_str((Vm_str(:,t)+dVm_str)>=V_thres,t) = Vpeak;
    Vm_str((Vm_str(:,t)+dVm_str)>=V_thres,t+1) = Vrest;
end
%GP
Vm_gp(:,t+1) = Vm_gp(:,t)+dVm_gp;
if sum(delta_gp)>=1
    Vm_gp((Vm_gp(:,t)+dVm_gp)>=V_thres,t) = Vpeak;
    Vm_gp((Vm_gp(:,t)+dVm_gp)>=V_thres,t+1) = Vrest;
end
%SNr
Vm_snr(:,t+1) = Vm_snr(:,t)+dVm_snr;
if sum(delta_snr)>=1
    Vm_snr((Vm_snr(:,t)+dVm_snr)>=V_thres,t) = Vpeak;
    Vm_snr((Vm_snr(:,t)+dVm_snr)>=V_thres,t+1) = Vrest;
end 

%Save variables 
Isyn_gp_out = [Isyn_gp_out,Isyn_gp];
Isyn_snr_out = [Isyn_snr_out,Isyn_snr];
end
toc