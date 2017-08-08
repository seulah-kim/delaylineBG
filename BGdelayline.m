%% BGdelayline.m
%This is the main code for simulating BG delay line model
function [Vm_gp, Vm_snr, Vm_str, Isyn_gp_out, Isyn_snr_out] = BGdelayline(varargin);

tic
close all

p = inputParser;	% construct input parser object

%%Network size
p.addParameter('r', 10); % convergence ratio between layers
p.addParameter('n', 100); % number of str neurons

%%Simulation time
dt = 0.001; %(s)
t_span = 0:dt:20;

%%Define Constants
R = 10;        % Membrane Resistance (MOhm)
Erev_i = -85;   % Synaptic reversal potential(mV)
I_const = 0.1;   % General scale for constant excitatory input(pA)
Vrest = -70;    % Resting potential(mV)
Vpeak = 15;     % Peak potential(mV)
V_thres = -64;  % Threshold voltage(mV)
%Synaptic
p.addParameter('tau_syn',5); % synaptic decay constant (ms)
p.addParameter('prob_syn', 0.35); % probability of successful synaptic transmission
g_uni = 100 ; % conductance of a single synapse

%%Input current to Str layer
randNum = round((20-10).*rand(1,1)+10);
randNum = 5;	% timing of Stimulation
tStim = [randNum:dt:randNum+0.5]; % (s) 
IextRatio_gpsnr = 10/3;
p.addParameter('stimCellsPer',100);	% Percentage of Str cells receiving stimulation

%% Parse and validate input arguments
p.parse(varargin{:}); 

% Assign variables from parsed input arguments
r = p.Results.r;
n = p.Results.n;
tau_syn = p.Results.tau_syn;
prob_syn = p.Results.prob_syn;
stimCellsPer = p.Results.stimCellsPer;

%%Initialize variables
%Cellular
tau_cell_str = ones(n,1);   % cell decay constant (ms)
tau_cell_gp = ones(n/r,1);   % cell decay constant (ms)
tau_cell_snr = ones(n/r.^2,1);   % cell decay constant (ms)
%homogeneous conductance
g_str2gp = 0.05*ones(n/r,1);    % (nS)	
g_gp2snr = 0.3*ones(n/r.^2,1);
%heterogeneous Vstart
Vm_str = Vrest+2*randn(n,1);  % (mV)
Vm_gp = Vrest+2*randn(n/r,1);
Vm_snr = Vrest+2*randn(n/r.^2,1);
del_str = zeros(n,1);   % binary
del_gp = zeros(n/r,1);
Isyn_gp_out =[];
Isyn_snr_out =[];
Iext_str = zeros(n,length(t_span));
stimCells=datasample(1:n,n*stimCellsPer/100,'Replace',false);	% Random sample of str cells for receiving stimulus, defined by percentage 

Iext_str(stimCells,ismember(t_span,tStim)) = 0*I_const;  % constant input to striatum ( 0.6~0.61Iconst for relatively quiet str cells) (in pA)

%%Simulation
for t = 1:length(t_span)
%Striatum
dVm_str = 1./tau_cell_str.*(-1*(Vm_str(:,t)-Vrest*ones(n,1)) + Iext_str(:,t)*R)*dt;   % dV/dt = 1/tau*(-V +IR) 

%GP
synSuccess_str2gp = double(rand(n,n/r)<prob_syn);	% flipping coin: n x n/r binary matrix 
Isyn_gp = g_str2gp(:,t).*(Vm_gp(:,t)-Erev_i*ones(n/r,1));  % synaptic (pA)
Iext_gp = 30*I_const*ones(n/r,1);	% external input (pA)
dg_str2gp = (-g_str2gp(:,t)./tau_syn + transpose(del_str'*synSuccess_str2gp*g_uni/n))*dt;
dVm_gp =1./tau_cell_gp.* (-(Vm_gp(:,t)-Vrest*ones(n/r,1)) - Isyn_gp*R+Iext_gp*R)*dt;

%SNr
synSuccess_gp2snr = double(rand(n/r,n/r.^2)<prob_syn);	% flipping coin: n/r x n/r^2 binary matrix
Isyn_snr = g_gp2snr(:,t).*(Vm_snr(:,t)-Erev_i.*ones(n/r.^2,1));	% synaptic (pA)
Iext_snr = round(30*IextRatio_gpsnr,2)*I_const*ones(n/r.^2,1);	% external input (pA)
dg_gp2snr = (-g_gp2snr(:,t)./tau_syn + transpose(del_gp'*synSuccess_gp2snr*g_uni/(n/r)))*dt;
dVm_snr = 1./tau_cell_snr.*(-(Vm_snr(:,t)-Vrest*ones(n/r.^2,1)) - Isyn_snr*R+Iext_snr*R)*dt;

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
toc
end
