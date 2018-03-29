%% BGdelayline_SNrdelay.m
%This is the main code for simulating BG delay line model
function [Vm_snr] = BGdelayline_SNrdelay(varargin);

p = inputParser;	% construct input parser object

%%Network size
nSNr = 100; % number of SNr neurons

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
p.addParameter('Iexc_snr',60); % (pA)
p.addParameter('Igaba_snr',10); % (pA)

%% Parse and validate input arguments
p.parse(varargin{:}); 

% Assign variables from parsed input arguments

Iexc_snr = p.Results.Iexc_snr;
Igaba_snr = p.Results.Igaba_snr;

%%Initialize variables
%Cellular
tau_cell_snr = 0.01*ones(nSNr,1);   % cell decay constant (s)
%heterogeneous Vstart
Vm_snr = Vrest+5*randn(nSNr,1);

%%Simulation
for t = 1:length(t_span)
Isyn_snr = Igaba_snr*ones(nSNr,1);	% synaptic (nS x mV = pA)
if t> (1/dt) && t<(1.01/dt) % GPe ceases firing at 1s <= t <= 1.01s
    Isyn_snr = zeros(nSNr,1);
end
Iext_snr = Iexc_snr*ones(nSNr,1);	% external input (pA) 
dVm_snr = 1./tau_cell_snr.*(-(Vm_snr(:,t)-Vrest*ones(nSNr,1))+(Iext_snr-Isyn_snr)/1000*R)*dt + 20.*randn(nSNr,1).*sqrt(dt); %convert pA to nA by dividing by 1000. nA*MOhm = mV

%Check for spikes in the cells
ind_snr = find((Vm_snr(:,t)+dVm_snr)>=V_thres);

%Update voltages based on spike pattern 

%SNr
Vm_snr(:,t+1) = Vm_snr(:,t)+dVm_snr;
if ~isempty(ind_snr) 
    Vm_snr(ind_snr,t) = Vpeak;
    Vm_snr(ind_snr,t+1) = Vrest;
end 

end
end
