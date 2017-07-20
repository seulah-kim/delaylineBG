function dydt = Vm_diffeq(t,y,n,r)
% y = [g_str1, g_str2,... g_gp1,...,g_snr1,
% Variables
% g_i = y(1:(n+n/r+n/r.^2));
g_str2gp = y(1:n/r);
g_gp2snr = y((n/r+1):(n/r+n/r.^2));
% Vm_i = y((n+n/r+n/r.^2+1):length(y));
Vm_str = y((n/r+n/r.^2+1):(n+n/r+n/r.^2));
Vm_gp = y((n+n/r+n/r.^2+1):(n+2*n/r+n/r.^2));
Vm_snr = y((n+2*n/r+n/r.^2+1):length(y));
% Jump Vm_i if it exceeds threshold voltage?




% Define Constants
%%cellular
tau_cell = 1;
Cm = 1;
Erev_i = -85;
I_const = 30; 
Vrest = -70;
Vpeak = 15; 
Vm_thres = -64;
k =0.1;    % steepness of logistic function
xd = 5;    % x location of midpoint of logistic function
L = 10;     % max value of curve
synNonlin = @(x)[L./(1+exp(-k*(x-xd)))];   % check if this is in the working range

%%synaptic
tau_syn =1;
prob_syn = 0.3; % probability of successful synaptic transmission
g_uni = 0.01 ; % conductance of a single synapse 
% synSuccess_i = rand(length(Vm_i),1)<prob_syn; % flip a coin: success (1) or failure (0)
% probSyn_gp = rand(n./r)>prob_syn; 

% probSyn_gp = rand(length(t_span),n./r)>prob_syn; 
% sum(spk_count_str.*synSuccess

%Striatum
delta_str = double(Vm_str>Vm_thres);
Vm_str(Vm_str>Vm_thres) = Vrest;   % bring it back to Vrest
% Vm_str(Vm_str>=Vpeak) = V
% I_syn =  -g_i.*(Vm_i-Erev_i);

Iext_str = 0.1*I_const*ones(length(Vm_str),1);  % constant input to striatum

dVm_str = -(Vm_str-Vrest*ones(length(Vm_str),1))/tau_cell + Iext_str;   % no nonlinear function here!

%GP
delta_gp = double(Vm_gp>Vm_thres);
Vm_gp(Vm_gp>Vm_thres) = Vrest;   % bring it back to Vrest

synSuccess_str2gp = double(rand(length(Vm_str),length(Vm_gp))<prob_syn);
Isyn_gp = -g_str2gp.*(Vm_gp-Erev_i*ones(length(Vm_gp),1));
Iext_gp = Isyn_gp+I_const*ones(length(Vm_gp),1);
dg_str2gp = -g_str2gp./tau_syn + transpose(delta_str'*synSuccess_str2gp*g_uni);

dVm_gp = -(Vm_gp-Vrest*ones(length(Vm_gp),1))/tau_cell + synNonlin(Iext_gp);

%SNr
% delta_snr = Vm_snr>Vm_thres;
Vm_snr(Vm_snr>Vm_thres) = Vrest;   % bring it back to Vrest
synSuccess_gp2snr = double(rand(length(Vm_gp),length(Vm_snr))<prob_syn);
Isyn_snr = -g_gp2snr.*(Vm_snr-Erev_i.*ones(length(Vm_snr),1));
Iext_snr = Isyn_snr+I_const*ones(length(Vm_snr),1);
% dg_gp = -g_snr./tau_syn + sum(delta_snr.*synSuccess_i.*g_uni);
dg_gp2snr = -g_gp2snr./tau_syn + transpose(delta_gp'*synSuccess_gp2snr*g_uni);
dVm_snr = -(Vm_snr-Vrest*ones(length(Vm_snr),1))/tau_cell + synNonlin(Iext_snr);

% dg_i= -g_i./tau_syn + sum(delta_i.*synSuccess_i.*g_uni);
% I_ext = I_syn + I_const*ones(length(Vm_i),1);
% dVm_i = -(Vm_i-Vrest*ones(length(Vm_i),1))/tau_cell + 20./(1+exp(-I_ext));
% dydt = [dg_i;dVm_i];
dydt = [dg_str2gp;dg_gp2snr;dVm_str;dVm_gp;dVm_snr];