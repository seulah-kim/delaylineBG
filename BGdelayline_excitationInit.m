%% BGdelayline_excitationInit.m
%This is the main code for simulating BG delay line model
function [g_gp2snr_out , Isyn_snr, del_gp_out] = BGdelayline_excitationInit(varargin);
p = inputParser;	% construct input parser object

%%Network size
p.addParameter('r', 10); % convergence ratio between layers
p.addParameter('n', 100); % number of gp neurons

%%Simulation time
dt = 0.0001; %(s) dt is 0.1ms resolution
t_span = 0:dt:3;

%%Define Constants
R = 100;        % Membrane Resistance (MOhm) ** changed from 100 1/8/18 
Erev_i = 0;   % Synaptic reversal potential(mV) for AMPAR
%I_const = 0.1;   % General scale for constant excitatory input(nA)
Vrest = -70;    % Resting potential(mV)
Vpeak = 15;     % Peak potential(mV)
V_thres = -64;  % Threshold voltage(mV)
%Synaptic
p.addParameter('tau_syn',0.005); % synaptic decay constant (s) because RK method will multiply everything by dt
p.addParameter('prob_syn', 0.35); % probability of successful synaptic transmission
p.addParameter('prob_syn_gp2snr',0.35); % probability of successful synaptic transmission for gp2snr
g_uni = 300; % synaptic conductance (pS) -- do not know which value to use
p.addParameter('g_gp2snr_i',0.0006);


%%Input current to Str layer
%randNum = round((20-10).*rand(1,1)+10);
randNum = 1;	% timing of Stimulation
tStim = (1/dt):((1+0.01)/dt); % (s) 10 ms activation
%IextRatio_gpsnr = 4; %77/8 or 1?
p.addParameter('stimCellsPer',35);	% Percentage of Str cells receiving stimulation
p.addParameter('I_exc_gp',60);		% Total excitatory input to GP cells controls firing rate  60pA (pA) -> equivalent to 6mV with R=100MOhm
p.addParameter('I_exc_snr',60);		% Total excitatory input to GP cells controls firing rate  60pA (pA) -> equivalent to 6mV with R=100MOhm

p.addParameter('connectivity','random');   % GPe2SNr connectivity. Default is all cells convering to 1.
%% Parse and validate input arguments
p.parse(varargin{:}); 

% Assign variables from parsed input arguments
r = p.Results.r;
n = p.Results.n;
tau_syn = p.Results.tau_syn;
prob_syn = p.Results.prob_syn;
stimCellsPer = p.Results.stimCellsPer;
I_exc_gp = p.Results.I_exc_gp;
I_exc_snr = p.Results.I_exc_snr;
prob_syn_gp2snr = p.Results.prob_syn_gp2snr;
g_gp2snr_i = p.Results.g_gp2snr_i;
connectivity = p.Results.connectivity;

%%Initialize variables
%Cellular
tau_cell_gp = 0.01*ones(n,1);   % cell decay constant (s)
tau_cell_snr = 0.01*ones(n/r,1);   % cell decay constant (s)
%homogeneous conductance
g_gp2snr = g_gp2snr_i*ones(n/r,length(t_span));	% steady state value changes depending on g_uni (3.5 for Iext=80, 6.5 for Iext=150, 12.6 for Iext=280)
%heterogeneous Vstart
Vm_gp = Vrest+5*randn(n,1);
Vm_snr = Vrest+5*randn(n/r,1);
del_gp = zeros(n,1);
W_gp = zeros(n/r,n/r);

switch connectivity
    case 'all'
    case 'random'   % picks selected number of inputs but non-overlapping connection
        nInputs = r;    % number of inputs per cell
        
        for ii=1:n/r
            connection = randsample([1:n],nInputs);
            gp2snr_connectivity(connection,ii) = 1; %connectivity
            W_gp(:,ii)= connection; %GPe input location to each SNr cell
        end 
            
    case 'segregated'   % picks selected number of inputs and converge in non-overlapping
        nInputs = r;    % number of inputs per cell

        for ii=1:n/r
            connection = ((ii-1)*nInputs+1):ii*nInputs;
            gp2snr_connectivity(connection,ii) = 1; %connectivity
            W_gp(:,ii)= connection; %GPe input location to each SNr cell
        end 
        
end

W_gp2snr = zeros(n,n/r);
Isyn_snr_out =[];

stimCells=datasample(1:n,n*stimCellsPer/100,'Replace',false);	% Random sample of str cells for receiving stimulus, defined by percentage 

% input to Str
Iext_gp = I_exc_gp*ones(n,1);  % external input to Str (1000 pA) 

%%Simulation
for t = 1:length(t_span)
%GPe
dVm_gp = dt./tau_cell_gp.*(-Vm_gp(:,t)+Vrest*ones(n,1)+Iext_gp/1000*R) + 20.*randn(n,1).*sqrt(dt);   % dV/dt = 1/tau*(-V +IR). Here white noise input is described by Langenvin eqn.  

%SNr
switch connectivity
    case 'all'
        W_gp2snr = double(rand(n,n/r)<prob_syn_gp2snr);
    otherwise
        for n_snr = 1:n/r
            synSuccess_gp2snr = double(rand(nInputs,1)<prob_syn_gp2snr);	% flipping coin 10 times
            W_gp2snr(W_gp(:,n_snr),n_snr) = synSuccess_gp2snr;     % Assigning flipping results back to weight matrix
        end
end

Isyn_snr = g_gp2snr(:,t).*(Vm_snr(:,t)-Erev_i.*ones(n/r,1));	% synaptic (nS x mV = pA)
Iext_snr = I_exc_snr*ones(n/r,1);	% external input (pA) ** change this value to basal Isyn_snr input based on gp f.r.
dg_gp2snr = -g_gp2snr(:,t)./tau_syn.*dt + transpose(del_gp'*W_gp2snr*0.001*g_uni*(10/nInputs)); %nS (g_uni is converted from pS to nS)
dVm_snr = 1./tau_cell_snr.*(-(Vm_snr(:,t)-Vrest*ones(n/r.^2,1))+(Iext_snr-Isyn_snr)/1000*R)*dt + 20.*randn(n/r.^2,1).*sqrt(dt); %convert pA to nA by dividing by 1000. nA*MOhm = mV

%Update conductances

g_gp2snr(:,t+1) = g_gp2snr(:,t)+dg_gp2snr;

%Check for spikes in the cells
del_gp = (Vm_gp(:,t)+dVm_gp)>=V_thres;

ind_gp = find(del_gp);
ind_snr = find((Vm_snr(:,t)+dVm_snr)>=V_thres);

%Update voltages based on spike pattern 
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
Isyn_snr_out = [Isyn_snr_out,Isyn_snr];
end

%reset initial g_gp2snr
g_gp2snr_out = mean(mean(g_gp2snr(:,end-1/dt:end)));
%determine net inhibitory synaptic current during last 1 s.
Isyn_snr = mean(mean(Isyn_snr_out(:,end-1/dt:end)));
end
