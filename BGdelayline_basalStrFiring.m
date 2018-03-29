%% BGdelayline.m
%This is the main code for simulating BG delay line model
function [Vm_gp, Vm_snr, Vm_str, Isyn_gp_out, Isyn_snr_out] = BGdelayline_basalStrFiring(varargin);
p = inputParser;	% construct input parser object

%%Network size
p.addParameter('r', 10); % convergence ratio between layers
p.addParameter('n', 100); % number of str neurons

%%Simulation time
dt = 0.0001; %(s) dt is 0.1ms resolution
t_span = 0:dt:3;

%%Define Constants
R = 100;        % Membrane Resistance (MOhm) ** changed from 100 1/8/18 
Erev_i = -85;   % Synaptic reversal potential(mV)
%I_const = 0.1;   % General scale for constant excitatory input(nA)
Vrest = -70;    % Resting potential(mV)
Vpeak = 15;     % Peak potential(mV)
V_thres = -64;  % Threshold voltage(mV)
%Synaptic
p.addParameter('tau_syn',0.005); % synaptic decay constant (s) because RK method will multiply everything by dt
p.addParameter('prob_syn', 0.35); % probability of successful synaptic transmission
p.addParameter('prob_syn_gp2snr',0.35); % probability of successful synaptic transmission for gp2snr
g_uni = 1000; % synaptic conductance (pS) -- do not know which value to use
p.addParameter('g_gp2snr_i',0.0006);


%%Input current to Str layer
%randNum = round((20-10).*rand(1,1)+10);
randNum = 1;	% timing of Stimulation
tStim = (1/dt):((1+0.01)/dt); % (s) 
%IextRatio_gpsnr = 4; %77/8 or 1?
p.addParameter('stimCellsPer',35);	% Percentage of Str cells receiving stimulation
p.addParameter('I_exc_gp',60);		% Total excitatory input to GP cells controls firing rate  60pA (pA) -> equivalent to 6mV with R=100MOhm

%% Parse and validate input arguments
p.parse(varargin{:}); 

% Assign variables from parsed input arguments
r = p.Results.r;
n = p.Results.n;
tau_syn = p.Results.tau_syn;
prob_syn = p.Results.prob_syn;
stimCellsPer = p.Results.stimCellsPer;
I_exc_gp = p.Results.I_exc_gp;
prob_syn_gp2snr = p.Results.prob_syn_gp2snr;
g_gp2snr_i = p.Results.g_gp2snr_i;


%%Initialize variables
%Cellular
tau_cell_str = 0.01*ones(n,1);   % cell decay constant (s)
tau_cell_gp = 0.01*ones(n/r,1);   % cell decay constant (s)
tau_cell_snr = 0.01*ones(n/r.^2,1);   % cell decay constant (s)
%homogeneous conductance
g_str2gp = 0*ones(n/r,length(t_span));    % (nS)	
g_gp2snr = g_gp2snr_i*ones(n/r.^2,length(t_span));	% steady state value changes depending on g_uni (3.5 for Iext=80, 6.5 for Iext=150, 12.6 for Iext=280)
%heterogeneous Vstart
Vm_str = Vrest+5*randn(n,1);  % (mV)
Vm_gp = Vrest+5*randn(n/r,1);
Vm_snr = Vrest+5*randn(n/r.^2,1);
del_str = zeros(n,1);   % binary
del_gp = zeros(n/r,1);

Isyn_gp_out =[];
Isyn_snr_out =[];
%Iext_str = 20.*randn(n,length(t_span));   % noise standard deviation is 20pA 
Iext_str = zeros(n,length(t_span));   % noise standard deviation is 20pA 
stimCells=datasample(1:n,n*stimCellsPer/100,'Replace',false);	% Random sample of str cells for receiving stimulus, defined by percentage 

basalstimCells = datasample(1:n,n*0.1,'Replace',false);

% input to Str
Iext_str(stimCells,tStim) = 1000 ;  % external input to Str (1000 pA) 

%%Simulation
for t = 1:length(t_span)
%Striatum
dVm_str = dt./tau_cell_str.*(-Vm_str(:,t)+Vrest*ones(n,1)+Iext_str(:,t)/1000*R) + 20.*randn(n,1).*sqrt(dt);   % dV/dt = 1/tau*(-V +IR). Here white noise input is described by Langenvin eqn.  

%GP
synSuccess_str2gp = double(rand(n,n/r)<prob_syn);	% flipping coin: n x n/r binary matrix 
Isyn_gp = g_str2gp(:,t).*(Vm_gp(:,t)-Erev_i*ones(n/r,1));  % synaptic (nS x mV = pA)
Iext_gp = I_exc_gp*ones(n/r,1);	% external input (pA) -  original value 60
dg_str2gp = -g_str2gp(:,t)./tau_syn.*dt + transpose(del_str'*synSuccess_str2gp*0.001*g_uni); %nS (g_uni is converted from pS to nS)
dVm_gp =1./tau_cell_gp.* (-(Vm_gp(:,t)-Vrest*ones(n/r,1)) + (Iext_gp-Isyn_gp)/1000*R)*dt+ 20.*randn(n/r,1).*sqrt(dt); %convert pA to nA by dividing by 1000. nA*MOhm = mV

%SNr
synSuccess_gp2snr = double(rand(n/r,n/r.^2)<prob_syn_gp2snr);	% flipping coin: n/r x n/r^2 binary matrix
Isyn_snr = g_gp2snr(:,t).*(Vm_snr(:,t)-Erev_i.*ones(n/r.^2,1));	% synaptic (nS x mV = pA)
Iext_snr = 65*ones(n/r.^2,1);	% external input (pA) ** change this value to basal Isyn_snr input based on gp f.r.
dg_gp2snr = -g_gp2snr(:,t)./tau_syn.*dt + transpose(del_gp'*synSuccess_gp2snr*0.001*g_uni); %nS (g_uni is converted from pS to nS)
dVm_snr = 1./tau_cell_snr.*(-(Vm_snr(:,t)-Vrest*ones(n/r.^2,1))+(Iext_snr-Isyn_snr)/1000*R)*dt + 20.*randn(n/r.^2,1).*sqrt(dt); %convert pA to nA by dividing by 1000. nA*MOhm = mV

%Update conductances
g_str2gp(:,t+1) = g_str2gp(:,t)+dg_str2gp;
g_gp2snr(:,t+1) = g_gp2snr(:,t)+dg_gp2snr;

%Check for spikes in the cells
del_str = (Vm_str(:,t)+dVm_str)>=V_thres;   % keep the array for matrix operation in line 78 (dg_str2gp)
del_gp = (Vm_gp(:,t)+dVm_gp)>=V_thres;

ind_str = find(del_str);    % indices of spikes
ind_gp = find(del_gp);
ind_snr = find((Vm_snr(:,t)+dVm_snr)>=V_thres);

%Update voltages based on spike pattern 
%Str
Vm_str(:,t+1) = Vm_str(:,t)+dVm_str;
if ~isempty(ind_str)
    Vm_str(ind_str,t) = Vpeak;
    Vm_str(ind_str,t+1) = Vrest;
end

%GP
Vm_gp(:,t+1) = Vm_gp(:,t)+dVm_gp;
if ~isempty(ind_gp) 
    Vm_gp(ind_gp,t) = Vpeak;
    Vm_gp(ind_gp,t+1) = Vrest;
end

%SNr
Vm_snr(:,t+1) = Vm_snr(:,t)+dVm_snr;
if ~isempty(ind_snr) 
    Vm_snr(ind_snr,t) = Vpeak;
    Vm_snr(ind_snr,t+1) = Vrest;
end 

%Save variables 
Isyn_gp_out = [Isyn_gp_out,Isyn_gp];
Isyn_snr_out = [Isyn_snr_out,Isyn_snr];
end

end
