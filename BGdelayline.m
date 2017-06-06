%% BG delay line model

r = 10; % convergence ratio between layers
n = 100; % number of str neurons

%%Initial Conditions
Vrest = -70;  % mV
Cm = 1; % Membrane capacitance
g_in = 0.3; % inhibitory synaptic conductance
g_ex = 1; % excitatory synaptic conductance (optional)
tau_cell = 5;   % cell decay constant (ms)
tau_syn = 5;  % synaptic decay constant (ms)
E_gp = -85; % GP-synaptic reversal potential (mV)
E_snr = -85; % SNr-synaptic reversal potential (mV)
spk_count_str = zeros(n,1);
spk_count_gp = zeros(n./r,1);
spk_count_snr = zeros(n./r.^2,1);
g_unitary = 1;  % conductance of each synapse (mS)
prob_syn = 0.8; % probability of successful synaptic transmission

%% Feedforward inhibition
% s = rng;


V_str = Vrest.*ones(n,1);
V_gp = Vrest.*ones(n./r,1);
V_snr = Vrest.*ones(n./r.^2,1);
%Simulation time
t_span = 0:0.001:10;

%Probablistic synapse
rng(1)  % for repeatable simulations-- set random number generator with sequence #1
probSyn_str = rand(length(t_span),n)>prob_syn; % flip a coin: success (1) or failure (0)
probSyn_gp = rand(length(t_span),n./r)>prob_syn; 


[TOUT,YOUT] = ode45(@g_diffeq,t_span,g_str)

%%total synaptic conductance
g_str=sum(spk_count_str.*probSyn_str      % how to advance through probSyn_str??