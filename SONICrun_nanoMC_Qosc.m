function [TTime] = SONICrun_nanoMC_Qosc(Tsim,MODE,USpstart,USpd,USfreq,USdc,USprf,USisppa,ESpstart,ESpd,...
    ESdc,ESprf,ESisppa,PLOT,model,USibegin,USiend,SearchMode,aBLS,fBLS,varargin)
% SONIC model that  accounts for charge oscillations by fourier series components (see Paper 3)
coder.extrinsic('nakeinterp1');
if nargin < 20
fBLS = '1';
if nargin < 19
aBLS = '32e-9';
if nargin < 18
error('Not enough input arguments');
end
end
end
proteinMode = '0'; threshMode = '0'; modeStr = ''; gateMultip = '1'; NFS = '2';
if length(varargin) > 4 , disp('Warning: extra input parameters will be ignored'); end
if length(varargin) == 4,  proteinMode = varargin{1}; threshMode = varargin{2}; gateMultip = varargin{3}; NFS = varargin{4};
 modeStr = ['-(proteinMode,threshMode,gateMultip,NFS)=(' proteinMode ',' threshMode ',' gateMultip ',' NFS ')']; end
if length(varargin) == 3,  proteinMode = varargin{1}; threshMode = varargin{2}; gateMultip = varargin{3};
 modeStr = ['-(proteinMode,threshMode,gateMultip)=(' proteinMode ',' threshMode ',' gateMultip ')']; end
if length(varargin) == 2, proteinMode = varargin{1}; threshMode = varargin{2};
 modeStr = ['-(proteinMode,threshMode)=(' proteinMode ',' threshMode ')']; end
if length(varargin) == 1, proteinMode = varargin{1}; modeStr = ['-(proteinMode)=(' proteinMode ')']; end
% SONICrun_nanoMC is a MultiCompartmental nanoscale model (see Fig. 3 and
% Fig. 10 in Lemaire et al.). Based on SONICrun, which is a speed-up version 
% of funPES based on multi-scale optimization with effective parameters. Based on:
% SONIC (Lemaire et al., 2018), Piezoelectric solver (see Plaksin et al.,2014; Plaksin et al.,2016;
% Krasovitski et al., 2011, PhD thesis Plaksin, 2016, Tarnaud et al.,2018).

% This function should be used, after tables are made with SONICtabulate.m.
% -- Parallellized function for HPC calculations ---
% 1. Output is saved on the current path
% Possible output: MODE=1 -> Isppa_thresh for 1 AP
%                  MODE=2 -> Array of AP-times
% In MODE==1 USisppa is ignored
% 2. Input: US = ultrasonic stimulus, ES = electrical stimulus
% pstart = Pulse start, pd = pulse duration, freq = frequency, dc = duty
% cycle, isppa = spatial peak pulse average intensity, plot: if 1 make a
% plot -> use 0 for HPC simulations, if 2: no plots but save necessary data
% MODEL = neuron model USibegin and USiend are used when MODE=1 to indicate the
% search-range=[USibegin,USiend] for the intensity threshold 
% SearchMode = 0 indicates that ~any(~search-range) = 0,
% while SearchMode = 1 handles the case ~any(~search-range) = 1.
% ProteinMode: coverage of proteins: (mode 0: full coverage; mode
% 1: partial coverage, leak current has full coverage; mode 2: all gates partial
% coverage including leak current)
% threshMode: (0: neuron excitation during Tsim for threshold
% determination; 1: only neuron excitation during stimulus duration).
% gateMultip is a multiplier, applied to the conductivity gains that are
% reduced if proteinMode ~= 0. E.g. if gains are actually determined at
% fBLS_exp = 0.75, the local gains are 4 times the HH-gains.
% 3. Units: SI-units are used, except for the membrane voltage [mV].
% -> All input variables are strings, because this works more fluently
% -PBS files for HPC simulations:
Tsim = str2double(Tsim); MODE = str2double(MODE); USpstart = str2double(USpstart);
USpd = str2double(USpd); USfreq = str2double(USfreq); USdc= str2double(USdc);
USprf = str2double(USprf); USipa = str2double(USisppa); ESpstart = str2double(ESpstart);
ESpd = str2double(ESpd); ESdc = str2double(ESdc); ESprf = str2double(ESprf);
ESipa = str2double(ESisppa); PLOT = str2double(PLOT); 
USibegin = str2double(USibegin); USiend = str2double(USiend);
SearchMode = str2double(SearchMode); aBLS = str2double(aBLS);
fBLS = str2double(fBLS); proteinMode = str2double(proteinMode); threshMode = str2double(threshMode); 
gateMultip = str2double(gateMultip); NFS = str2double(NFS);

SearchPrecision = 2; % If MODE=1, number of significant digits of the solution
% of the titration algorithm
MODEL = str2double(model);
% Choose model: 1=RS,2=FS,3=LTS,4=TC,5=RE,6=RS_FVX,7=FS_FVX,8=LTS_CAX,
% 9=STN, 10=Th-RT, 11=GPi, 12=GPe, 13=Str-MSN, 14 = HH

% RS,FS,LTS,RS_FVX,FS_FVX and LTS_CAX are Pospischil models
% RE and TC are Destexhe-based models
% Th-RT, GPi, GPe are based on Rubin-Terman
% STN is based on Otsuka et al.
% Str-MSN is based on McCarthy et al.

if any(MODEL == [10,11,12])
fprintf(['\n Warning: SONIC multi-scale model (modelnr %u) makes use of instantaneous gate parameters ' ...
'(x=xinf). Implementation will assume x=xinf(Veff). \n'],MODEL);
% This implementation might not be the most accurate. 
% 1) Note that if taux << 1/USfreq - as is implicit in funPES - then SONIC-averaging 
% over the period of the non-linear membrane current would be most accurate.
% E.g. : I_Na = g_Na m_inf(V)^3 (V-E_Na) -> effective parameters for
% m_inf(Q/Cm)^3*(Q/Cm) and m_inf(Q/Cm)^3 are required. 
% 2) If taux > 1/USfreq, while taux << dtES (motivating the assumption x =
% xinf), then from dx/dt = (x_inf-x)/taux it follows that x = xeff(Q) would
% be a good approximation (assumption: taux == cte) 
end

DISPLAY = 0; DISPLAYbase = 1;
% Display level. Note: higher display level will give more runtime information but will slow the program 
% DISPLAY = 0 -> No information displayed (use this option for HPC simulations)
% DISPLAY = 1 -> Display progress based on update nr. alone

% Discretisation parameters
dt = min(50e-6,0.1*min(USdc/USprf,ESdc/ESprf));            % Discretisation time (s)
atol = 1e-7; rtol = 1e-4; % absolute and relative VSVO-tolerances
maxRate = 1e6;        % (1/s). This is the maximal allowed rate constant of (a,apb). Set to 'inf' if no maxRate
% Physically protein gate opening/closing will have a minimal time delay, irrespective of voltage (~ 1/maxRate). 
% Computationally, very high maxRate will increase the stiffness of the set of equations, without important alterations in the solution set 
% (e.g. for all practical purposes, gate opening from 0->1 in 1 us is equal to infinitely fast opening) 
Tupdate = 500e-6;           % Update of fourier coefficients [s]
pretabulSONIC = 0;

tic;
% 1a. Multicompartmental parameters
deffV = 100*10^(-9);   % The effective depth beneath the membrane area for axial current calculations (m)
rhoi = 1;            % axoplasmatic resistivity (Ohm*m)
% 1b. General parameters
Qthresh = 0;            % Threshold for Q for AP-discrimination [nC/cm^2]
Cm0 = 0.01;				% Rest capacitance (F/m^2)	 
a=aBLS;	        		% radius leaflet boundary (m)
c = 1515;				% Speed of sound surrounding medium (m/s)
rhol = 1028;			% Density surrounding medium (kg/m^3)
Rg = 8.314;             % Universal gas constant (J/(K*mol))
Far = 96485.3329;       % Faraday constant (C/mol) 
%Temp = 309.15; 		    % Surrounding medium temperature (K)	
%Temp = 306.15;
Temp = 33+273.15;
deffCa = 100*10^(-9);   % The effective depth beneath the membrane area for calcium concentration calculations (m)
cCae = 2;               % Extracellular Ca2+ concentration (mol/m^3)
k1 = 2.5*10^22;         % k1 factor for hyperpolarization driven Ca regulation (M^-4*s^-1)
k2 = 0.4;               % k2 factor for hyperpolarization driven Ca regulation (s^-1)
k3 = 100;               % k3 factor for hyperpolarization driven Ca regulation (s^-1)
k4 = 1;                 % k4 factor for hyperpolarization driven Ca regulation (s^-1)
if MODEL==1             % Regular spiking neuron
Gna = 560;				% Maximal conductance of the Na-channel (S/m^2)
Vna = 50;				% Na nernst potential (mV)
Gk = 60;				% Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Gm = 0.75;				% Maximal conductance of the slow non-inactivating K-channel (S/m^2)
Vk = -90;				% Potassium nernst potential (mV)
Gl = 0.205;				% Maximal conductance of the non-voltage-dependent non-specific ion channel (S/m^2)
Vl = -70.3;				% Leak nernst potential (mV)
Vm0 = -71.9;            % rest membrane potential (mV)
VT = -56.2;             % Spike treshold adjustment parameter (mV)
taumax = 608*10^(-3);   % Decay time constant slow non-inactivating K+ (s)
modelName = 'RS';
elseif MODEL == 2       % Fast spiking neuron
Gna = 580;				% Maximal conductance of the Na-channel (S/m^2)
Vna = 50;               % Na nernst potential (mV)
Gk = 39;				% Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Gm = 0.787;				% Maximal conductance of the slow non-inactivating K-channel (S/m^2)
Vk = -90;				% Potassium nernst potential (mV)
Gl = 0.38;				% Maximal conductance of the non-voltage-dependent non-specific ion channel (S/m^2)
Vl = -70.4;				% Leak nernst potential (mV)
Vm0 = -71.4;            % Rest membrane potential (mV)
VT = -57.9;             % Spike treshold adjustment parameter (mV)
taumax = 502*10^(-3);   % Decay time constant slow non-inactivating K+ (s)
modelName = 'FS';
elseif MODEL == 3       % Low threshold spiking neuron
Gna = 500;				% Maximal conductance of the Na-channel (S/m^2)
Vna = 50;               % Na nernst potential (mV)
Gk = 40;				% Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Gm = 0.28;				% Maximal conductance of the slow non-inactivating K-channel (S/m^2)
Vk = -90;				% Potassium nernst potential (mV)
Gl = 0.19;				% Maximal conductance of the non-voltage-dependent non-specific ion channel (S/m^2)
Vl = -50;				% Leak nernst potential (mV)
GT = 4;                 % Maximal conductance of low-threshold Ca2+ channels (S/m^2)
VCa = 120;              % Nernst potential of Ca2+ (mV)
Vm0 = -54;              % Rest membrane potential (mV)
VT = -50;               % Spike treshold adjustment parameter (mV)
taumax = 4;             % Decay time constant slow non-inactivating K+ (s)
Vx = -7;				% Shift Ca2+ voltage (mV)
modelName = 'LTS';
elseif MODEL == 4       % Thalamocortical neuron
Gna = 900;				% Maximal conductance of the Na-channel (S/m^2)
Vna = 50;               % Na nernst potential (mV)
Gk = 100;				% Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Vk = -90;				% Potassium nernst potential (mV)
Gl = 0.1;				% Maximal conductance of the non-voltage-dependent non-specific ion channel (S/m^2)
Vl = -70;				% Leak nernst potential (mV)
GT = 20;                % Maximal conductance of low-threshold Ca2+ channels
GKL = 0.138;            % Maximal conductance of leak potassium currents
Gh = 0.175;             % Maximal conductance of hyperpolarization-activated mixed cationic current
ginc = 2;               % Locked gate relative conductance (-)
Vh = -40;               % Reversal potential of a hyperpolarization-activated mixed cationic urrent
Vm0 = -63.4;            % Rest membrane potential (mV)
tauCa = 5*10^(-3);      % Calcium decay time constant (s)
VT = -52;               % Spike treshold adjustment parameter (mV)
Vm0 = -63.4;            % Rest membrane potential (mV)
Vx = 0;                 % Shift Ca2+ voltage (mV)
modelName = 'TC';
elseif MODEL == 5       % Reticular thalamus neuron
Gna = 2000;				% Maximal conductance of the Na-channel (S/m^2)
Vna = 50;               % Na nernst potential (mV)
Gk = 200;				% Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Vk = -90;				% Potassium nernst potential (mV)
Gl = 0.5;				% Maximal conductance of the non-voltage-dependent non-specific ion channel (S/m^2)
Vl = -90;				% Leak nernst potential (mV)
GT = 30;                % Maximal conductance of low-threshold Ca2+ channels
Vm0 = -89.5;            % Rest membrane potential (mV)
tauCa = 5*10^(-3);      % Calcium decay time constant (s)
VT = -67;               % Spike treshold adjustment parameter (mV)
Vm0 = -89.5;            % Rest membrane potential (mV)
modelName = 'RE';
elseif MODEL == 6       % Regular spiking ferret visual cortex neuron
Gna = 500;				% Maximal conductance of the Na-channel (S/m^2)
Vna = 50;               % Na nernst potential (mV)
Gk = 50;				% Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Vk = -90;				% Potassium nernst potential (mV)
Gm = 0.7;				% Maximal conductance of the slow non-inactivating K-channel (S/m^2
Gl = 1;				    % Maximal conductance of the non-voltage-dependent non-specific ion channel (S/m^2)
Vl = -70;				% Leak nernst potential (mV)
Vm0 = -70.4;            % Rest membrane potential (mV)
VT = -55;               % Spike treshold adjustment parameter (mV)
taumax = 1;             % Decay time constant slow non-inactivating K+ (s)
modelName = 'RS_FVX';
elseif MODEL == 7       % Fast spiking ferret visual cortex neuron
Gna = 500;				% Maximal conductance of the Na-channel (S/m^2)
Vna = 50;               % Na nernst potential (mV)
Gk = 100;				% Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Vk = -90;				% Potassium nernst potential (mV)
Gm = 0;			    	% Maximal conductance of the slow non-inactivating K-channel (S/m^2
Gl = 1.5;				% Maximal conductance of the non-voltage-dependent non-specific ion channel (S/m^2)
Vl = -70;				% Leak nernst potential (mV)
Vm0 = -70;              % Rest membrane potential (mV)
VT = -55;               % Spike treshold adjustment parameter (mV)
taumax = 1;             % Decay time constant slow non-inactivating K+ (s)
modelName = 'FS_FVX';
elseif MODEL == 8       % Low threshold spiking cat association cortex neuron
Gna = 500;				% Maximal conductance of the Na-channel (S/m^2)
Vna = 50;               % Na nernst potential (mV)
Gk = 50;				% Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Vk = -90;				% Potassium nernst potential (mV)
Gm = 0.3;				% Maximal conductance of the slow non-inactivating K-channel (S/m^2
Gl = 0.1;				% Maximal conductance of the non-voltage-dependent non-specific ion channel (S/m^2)
Vl = -85;				% Leak nernst potential (mV)
GT = 4;                 % Maximal conductance of low-threshold Ca2+ channels
VCa = 120;              % Nernst potential of Ca2+ (mV)
Vm0 = -84.6;            % Rest membrane potential (mV)
VT = -55;               % Spike treshold adjustment parameter (mV)
taumax = 1;             % Decay time constant slow non-inactivating K+ (s)
Vx = -2;				% Shift Ca2+ voltage (mV)
modelName = 'LTS_CAX';
elseif MODEL == 9       % Subthalamic nucleus neurons -> Otsuka model
Gna=490;                % Maximal conductance of the Na-channel (S/m^2)
Vna=60;                 % Na nernst potential (mV)
Gk=570;                 % Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Vk=-90;                 % Potassium nernst potential (mV)
Gl=3.5;                 % Maximal conductance of the non-voltage-dependent non-specific ion channel (S/m^2)
Vl=-60;                 % Leak nernst potential (mV)
GT=50;                  % Maximal conductance of T-type low-threshold Ca2+ channels (S/m^2)
GL=150;                 % Maximal conductance of L-type high-threshold Ca2+ channels (S/m^2)
GA =50;                 % Maximal conductance of A-type K channel (S/m^2)
GCa=10;                 % Maximal conductance of Ca2+ activated K-channel (S/m^2)
Vm0=-58;                % Resting potential (mV)
tauCa = 0.5*10^(-3);    % Calcium decay time constant (s) 
rinf = @(kCA) 1./(1+exp(-(kCA-0.17)./0.08));  % rest r-value [-]
d2inf = @(kCA) 1./(1+exp((kCA-0.1)./0.02));        % rest d2-value [-]
taur = @(kCA) (10^(-3))*2;  % Time-constant [s] for r
taud2 = @(kCA) (10^(-3))*130; % Time constant [s] for d2
modelName = 'STN';
elseif MODEL == 10      % Thalamic Rubin-Terman based neuron
Gl=0.5;                 % Maximal conductance of non-specific non-voltage dependent ion channel (S/m^2)  
Vl=-70;                 % Leak nernst potential (mV)
Gna=30;                 % Maximal conductance of the Na-channel (S/m^2)
Vna=50;                 % Na nernst potential (mV)
Gk=50;                  % Maximal conductance of the K-channel (S/m^2)
Vk=-75;                 % K nernst potential (mV)
GT=50;                  % Maximal conductance of T-type low threshold Ca2+ channels (S/m^2)
Vm0=-65;                % Resting potential (mV)
VT=0;                   % Nernst potential of T-type low threshold Ca2+ channel (mV)
minf = @(V) 1./(1+exp(-(V+37)/7));            % rest m-value [-]
pinf = @(V) 1./(1+exp(-(V+60)/6.2));            % rest p-value [-]
modelName = 'Th_RT';
elseif MODEL == 11      % Globus pallidus internus neuron
Gl=1;                   % Maximal conductance of the non-voltage dependent non-specific ion channels (S/m^2)
Vl=-65;                 % Leak nernst potential (mV)
Gna=1200;               % Maximal conductance of the Na-channel (S/m^2)
Vna=55;                 % Na nernst potential (mv)
Gk=300;                 % Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Vk=-80;                 % K nernst potential (mV)
GT=5;                   % Maximal conductance of T-type low threshold Ca2+ channel (S/m^2)
VT=0;                   % Nernst potential of T-type Ca2+ channel (mV) 
% !! DON'T CONFUSE VT WITH SPIKE-THRESHOLD ADJUSTMENT PARAMETER IN OTHER MODELS !!
GCa=1.5;                % Maximal conductance of Ca2+ activated K-channel (S/m^2)
VCa=120;                % Nernst potential of Ca2+ activated K-channel (mV)
Gahp=100;               % Maximal conductance of the afterhyperpolarization K-channel (S/m^2)
Vahp=-80;               % Nernst potential of the AHP-K channel (mV)
Vm0 = -65;              % Resting potential (mV)
VT=0;                   % Nernst potential of T-type low threshold Ca2+ channel (mV)
minf = @(V) 1./(1+exp(-(V+37)./10));            % rest m-value [-]
ainf = @(V) 1./(1+exp(-(V+57)./2));            % rest a-value [-]
sinf = @(V) 1./(1+exp(-(V+35)./2));             % rest s-value [-]
modelName = 'GPi';
elseif MODEL == 12      % Globus pallidus externus neuron
Gl=1;                   % Maximal conductance of the non-voltage dependent non-specific ion channels (S/m^2)
Vl=-65;                 % Leak nernst potential (mV)
Gna=1200;               % Maximal conductance of the Na-channel (S/m^2)
Vna=55;                 % Na nernst potential (mv)
Gk=300;                 % Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Vk=-80;                 % K nernst potential (mV)
GT=5;                   % Maximal conductance of T-type low threshold Ca2+ channel (S/m^2)
VT=0;                   % Nernst potential of T-type Ca2+ channel (mV) 
GCa=1.5;                % Maximal conductance of Ca2+ activated K-channel (S/m^2)
VCa=120;                % Nernst potential of Ca2+ activated K-channel (mV)
Gahp=100;               % Maximal conductance of the afterhyperpolarization K-channel (S/m^2)
Vahp=-80;               % Nernst potential of the AHP-K channel (mV)
Vm0 = -65;              % Resting potential (mV)
VT=0;                   % Nernst potential of T-type low threshold Ca2+ channel (mV)
minf = @(V) 1./(1+exp(-(V+37)./10));            % rest m-value [-]
ainf = @(V) 1./(1+exp(-(V+57)./2));            % rest a-value [-]
sinf = @(V) 1./(1+exp(-(V+35)./2));             % rest s-value [-]
modelName = 'GPe';
elseif MODEL == 13      % Striatum medium spiny neuron
Gl=1;                   % Maximal conductance of the non-voltage dependent non-specific ion channel (S/m^2)
Vl=-67;                 % Leak nernst potential (mV)
Gna=1000;               % Maximal conductance of the Na-channel (S/m^2)
Vna=50;                 % Na nernst potential (mV)
Gk=800;                 % Maximal conductance of the delayed rectifier K-channel
Vk=-100;                % K delayed rectifier nernst potential (mV)
Vm=-100;                % K non-inactivating nernst potential (mV)
Vm0=-87;                % Resting potential (mV)
modelName = 'Str';
elseif MODEL == 14      % Hodgking-Huxley model 
Gna = 1200;             % Maximal conductance of the Na-channel (S/m^2)
Gk = 360;               % Maximal conductance of the K-channel (S/m^2)
Gl = 3;                 % Maximal conductance of the non-voltage dependent non-specific ion channels (S/m^2)
Vna = 45;               % Na nernst potential (mV)
Vk = -82;               % K nernst potential (mV)
Vl = -59.4;             % Leak nernst potential (mV)
Vm0 = -70;              % Resting potential (mV)
Qam = 3; Qbm = 3; Qan = 3; Qbn = 3; Qah = 3; Qbh = 3;
Qap = 3; Qbp = 3;       % Temperature factors (-)
Temp0 = 6.3+273.15;     % Default temperature (K)
modelName = 'HH';
end

f1rtV = struct;         % 1 dimensional rates in V
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 4 || MODEL == 5 || MODEL == 6 ||...
        MODEL == 7 || MODEL == 8
f1rtV.('a_m') =  @(V) -1000*0.32*((-4)*double((V-VT-13)==0)+...
    ((V-VT-13)./(double((V-VT-13)~=0).*exp(-((V-VT-13)/4))-1)));% Rate constant alpha_m [1/s]
f1rtV.('apb_m') = @(V) f1rtV.('a_m')(V)+1000*0.28*(5*double((V-VT-40)==0)+...
    (V-VT-40)./(double((V-VT-40)~=0).*exp(((V-VT-40)/5))-1)); % Rate constant beta_m [1/s]
f1rtV.('a_n') = @(V) -1000*0.032*((-5)*double((V-VT-15)==0)+...
    (V-VT-15)./(double((V-VT-15)~=0).*exp(-((V-VT-15)/5))-1));% Rate constant alpha_n [1/s]
f1rtV.('apb_n') = @(V) f1rtV.('a_n')(V)+1000*0.5*exp(-(V-VT-10)/40);                    % Rate constant beta_n [1/s]
f1rtV.('a_h') = @(V) 1000*0.128*exp(-((V-VT-17)/18));                % Rate constant alpha_h [1/s]
f1rtV.('apb_h') = @(V) f1rtV.('a_h')(V)+(1000*4)./(1+exp(-((V-VT-40)/5)));               % Rate constant beta_h [1/s]
end
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 6 || MODEL == 7 || MODEL == 8 
pinf = @(V) 1./(1+exp(-((V+35)/10)));                      % Rest p-value [-]
taup = @(V) taumax./(3.3*exp((V+35)/20)+exp(-(V+35)/20));  % Time-constant [s] for p 
f1rtV.('a_p') = @(V) pinf(V)./taup(V);
f1rtV.('apb_p') = @(V) 1./taup(V);
end
if MODEL == 3 || MODEL == 4 || MODEL == 8
sinf = @(V) 1./(1+exp(-(V+Vx+57)/6.2));					  % Rest s-value [-]
uinf = @(V) 1./(1+exp((V+Vx+81)/4));				          % Rest u-value [-] 
taus = @(V) 10^(-3)*((1/3.7)*(0.612+1./(exp(-((V+Vx+132)/16.7))+exp((V+Vx+16.8)/18.2))));  % Time-constant [s] for s
tauu = @(V) 10^(-3)*funtauu(V,Vx);                        % Time-constant [s] for u
f1rtV.('a_s') = @(V) sinf(V)./taus(V); f1rtV.('a_u') = @(V) uinf(V)./tauu(V);
f1rtV.('apb_s') = @(V) 1./taus(V); f1rtV.('apb_u') = @(V) 1./tauu(V);
end
if MODEL == 5
sinf = @(V) 1./(1+exp(-(V+52)/7.4));                       % Rest s-value [-]
uinf = @(V) 1./(1+exp((V+80)/5));                          % rest u-value [-]
taus = @(V) 10^(-3)*(1+0.33./(exp((V+27)/10)+exp(-(V+102)/15))); % Time constant [s] for s
tauu = @(V) 10^(-3)*(28.3+0.33./(exp((V+48)/4)+exp(-(V+407)/50))); % Time constant [s] for u
f1rtV.('a_s') = @(V) sinf(V)./taus(V); f1rtV.('a_u') = @(V) uinf(V)./tauu(V);
f1rtV.('apb_s') = @(V) 1./taus(V); f1rtV.('apb_u') = @(V) 1./tauu(V);
end
if MODEL == 4
winf = @(V) 1./(1+exp((V+75)/5.5));                        % Rest w-value [-]    
tauw = @(V) 10^(-3)./(exp(-14.59-0.086*V)+exp(-1.87+0.0701*V));  % Time-constant [s] for w
f1rtV.('a_w') = @(V) winf(V)./tauw(V);
f1rtV.('apb_w') = @(V) 1./tauw(V);
end
if MODEL == 9
minf = @(V) 1./(1+exp(-(V+40)./8));            % rest m-value [-]
hinf = @(V) 1./(1+exp((V+45.5)./6.4));          % rest h-value [-]
ninf = @(V) 1./(1+exp(-(V+41)./14));            % rest n-value [-]
pinf = @(V) 1./(1+exp(-(V+56)./6.7));          % rest p-value [-]
qinf = @(V) 1./(1+exp((V+85)./5.8));           % rest q-value [-]
ainf = @(V) 1./(1+exp(-(V+45)./14.7));          % rest a-value [-]
binf = @(V) 1./(1+exp((V+90)./7.5));            % rest b-value [-]
cinf = @(V) 1./(1+exp(-(V+30.6)./5));           % rest c-value [-]
d1inf = @(V) 1./(1+exp((V+60)./7.5));          % rest d1-value [-]

taum = @(V) (10^(-3))*(0.2+3./(1+exp((V+53)/0.7)));   % Time constant [s] for m
tauh = @(V) (10^(-3))*(24.5./(exp((V+50)/15)+exp(-(V+50)/16))); % Time-constant [s] for h
taun = @(V) (10^(-3))*(11./(exp((V+40)/40)+exp(-(V+40)/50))); % Time-constant [s] for n
taup = @(V) (10^(-3))*(5+0.33./(exp((V+27)/10)+exp(-(V+102)/15)));  % Time constant [s] for p
tauq = @(V) (10^(-3))*(400./(exp((V+50)/15)+exp(-(V+50)/16))); % Time constant [s] for q
taua = @(V) (10^(-3))*(1+1./(1+exp((V+40)/0.5))); % Time-constant [s] for a
taub = @(V) (10^(-3))*(200./(exp((V+60)/30)+exp(-(V+40)/10))); % Time constant [s] for b
tauc = @(V) (10^(-3))*(45+10./(exp((V+27)/20)+exp(-(V+50)/15))); % Time constant [s] for c
taud1 = @(V) (10^(-3))*(400+500./(exp((V+40)/15)+exp(-(V+20)/20))); % Time constant [s] for d1

f1rtV.('a_m') = @(V) minf(V)./taum(V); 
f1rtV.('a_h') = @(V) hinf(V)./tauh(V);
f1rtV.('a_n') = @(V) ninf(V)./taun(V);
f1rtV.('a_p') = @(V) pinf(V)./taup(V);
f1rtV.('a_q') = @(V) qinf(V)./tauq(V);
f1rtV.('a_a') = @(V) ainf(V)./taua(V);
f1rtV.('a_b') = @(V) binf(V)./taub(V);
f1rtV.('a_c') = @(V) cinf(V)./tauc(V);
f1rtV.('a_d1') = @(V) d1inf(V)./taud1(V);

f1rtV.('apb_m') = @(V) 1./taum(V); f1rtV.('apb_h') = @(V) 1./tauh(V);
f1rtV.('apb_n') = @(V) 1./taun(V); f1rtV.('apb_p') = @(V) 1./taup(V);
f1rtV.('apb_q') = @(V) 1./tauq(V); f1rtV.('apb_a') = @(V) 1./taua(V);
f1rtV.('apb_b') = @(V) 1./taub(V); f1rtV.('apb_c') = @(V) 1./tauc(V);
f1rtV.('apb_d1') = @(V) 1./taud1(V);
end
if MODEL == 10
hinf = @(V) 1./(1+exp((V+41)/4));            % rest m-value [-]
tauah = @(V) 10^(3)*0.128*exp(-(V+46)/18);      
taubh = @(V) 10^(3)*4./(1+exp(-(V+23)/5));
tauh = @(V) 1./(tauah(V)+taubh(V));             % Time constant [s] for h
rinf = @(V) 1./(1+exp((V+84)/4));            % rest r-value [-]
taur = @(V) (10^(-3))*0.15*(28+exp(-(V+25)/10.5));      % Time constant [s] for r

f1rtV.('a_h') = @(V) hinf(V)./tauh(V); f1rtV.('apb_h') = @(V) 1./tauh(V);
f1rtV.('a_r') = @(V) rinf(V)./taur(V); f1rtV.('apb_r') = @(V) 1./taur(V);
end
if MODEL == 11 || MODEL == 12
hinf = @(V) 1./(1+exp((V+58)./12));            % rest h-value [-]
ninf = @(V) 1./(1+exp(-(V+50)./14));           % rest n-value [-]
rinf = @(V) 1./(1+exp((V+70)./2));            % rest r-value [-]
tauh = @(V) (10^(-3)/0.05)*(0.05+0.27./(1+exp((V+40)/12))); % Time constant [s] for h
taun = @(V) (10^(-3)/0.1)*(0.05+0.27./(1+exp((V+40)/12)));  % Time constant [s] for n
taur = @(V) (10^(-3))*15;           % Time constant [s] for r

f1rtV.('a_h') = @(V) hinf(V)./tauh(V); f1rtV.('apb_h') = @(V) 1./tauh(V);
f1rtV.('a_n') = @(V) ninf(V)./taun(V); f1rtV.('apb_n') = @(V) 1./taun(V);
f1rtV.('a_r') = @(V) rinf(V)./taur(V); f1rtV.('apb_r') = @(V) 1./taur(V);
end
if MODEL == 13
am = @(V) 1000*0.32*(4*double((V+54)==0)+(54+V)./(1-double((V+54)~=0).*exp(-(V+54)/4))); % Rate constant alpha_m [1/s]
bm = @(V) 1000*0.28*(5*double((V+27)==0)+(27+V)./(-1+double((V+27)~=0).*exp((V+27)/5))); % Rate constant beta_m [1/s]
ah = @(V) 1000*0.128*exp(-(V+50)/18); % Rate constant alpha_h [1/s]
bh = @(V) 1000*4./(1+exp(-(V+27)/5)); % Rate constant beta_h [1/s]
an = @(V) 1000*0.032*(5*double((V+52)==0)+(52+V)./(1-double((V+52)~=0).*exp(-(V+52)/5))); % Rate constant alpha_n [1/s]
bn = @(V) 1000*0.5*exp(-(V+57)/40); % Rate constant beta_n [1/s]
ap = @(V) 1000*3.209*10^(-4)*(9*double((V+30)==0)+(30+V)./(1-double((V+30)~=0).*exp(-(V+30)/9))); % Rate constant alpha_p [1/s]
bp = @(V) -1000*3.209*10^(-4)*((-9)*double((V+30)==0)+(30+V)./(1-double((V+30)~=0).*exp((V+30)/9))); % Rate constant beta_p [1/s]

f1rtV.('a_m') = am; f1rtV.('a_h') = ah; f1rtV.('a_n') = an; f1rtV.('a_p') = ap;
f1rtV.('apb_m') = @(V) am(V)+bm(V); f1rtV.('apb_h') = @(V) ah(V)+bh(V); 
f1rtV.('apb_n') = @(V) an(V)+bn(V); f1rtV.('apb_p') = @(V) ap(V)+bp(V);
end
if MODEL == 14
am = @(V) 1000*Qam.^(0.1.*(Temp-Temp0)).*((2.5-0.1.*(V-Vm0))./(double((V-Vm0)~=25).*exp(2.5-0.1.*(V-Vm0))-1)+((V-Vm0)==25));
bm = @(V) 1000*Qbm.^(0.1.*(Temp-Temp0)).*4.*exp(-(V-Vm0)./18);
an = @(V) 1000*Qan.^(0.1.*(Temp-Temp0)).*(((0.1-0.01.*(V-Vm0))./(double((V-Vm0)~=10).*exp(1-0.1.*(V-Vm0))-1))+0.1*((V-Vm0)==10));
bn = @(V) 1000*Qbn.^(0.1.*(Temp-Temp0)).*0.125.*exp(-(V-Vm0)./80);
ah = @(V) 1000*Qah.^(0.1.*(Temp-Temp0)).*0.07.*exp(-(V-Vm0)./20);
bh = @(V) 1000*Qbh.^(0.1.*(Temp-Temp0)).*(1./(exp(3-0.1.*(V-Vm0))+1));

f1rtV.('a_m') = am; f1rtV.('a_h') = ah; f1rtV.('a_n') = an;
f1rtV.('apb_m') = @(V) am(V)+bm(V); f1rtV.('apb_h') = @(V) ah(V)+bh(V); 
f1rtV.('apb_n') = @(V) an(V)+bn(V); 
end
% ---------------------------ULTRASONIC SOURCE-----------------------------
if USdc == 1
USstep = @(t) double(t>=USpstart&t<=USpd+USpstart); 
else
USprp = (1/USprf);      % Pulse repetition period (s)
USstep = @(t) double(mod(t-USpstart,USprp)<=USdc*USprp).*double(t>=USpstart&t<=USpd+USpstart);
end
% --------------------------ELECTRONIC SOURCE------------------------------
if ESdc == 1
ESstep = @(t) double(t>=ESpstart&t<=ESpd+ESpstart); 
else
ESprp = (1/ESprf);      % Pulse repetition period (s)
ESstep = @(t) double(mod(t-ESpstart,ESprp)<=ESdc*ESprp).*double(t>=ESpstart&t<=ESpd+ESpstart);
end
ESi = @ (t) ESipa*ESstep(t);  % [A/m^2]
% -------------------------------------------------------------------------
% 2a. Multicompartmental important derived params or functions
RSI = (rhoi/(2*pi*deffV))*log((a+a/sqrt(fBLS))/a);           % Axial resistance (Ohm)

% 2b. General important functions
SONIC = load(['SONIC-' modelName '-QoscFourier' num2str(NFS) '-FourierIn' num2str(NFS) '.mat']); 
SONICtable = SONIC.SONICtable;
% 2.1 SONIC functions (rate, Veff, Zeff, Cmeff, ngend)
QmRange = SONICtable.QmRange; USPaRange = SONICtable.USPaRange; 
USfreqRange = SONICtable.USfreqRange; aBLSRange = SONICtable.aBLSRange;
DeltaQmRange = SONICtable.DeltaQmRange; psiQRange = SONICtable.psiQRange;

Veff6D = permute(SONICtable.Veff(:,:,:,:,1,:,:),[1 2 3 4 6 7 5]); Zeff6D = permute(SONICtable.Zeff(:,:,:,:,1,:,:),[1 2 3 4 6 7 5]); 
Cmeff6D = permute(SONICtable.Cmeff(:,:,:,:,1,:,:),[1 2 3 4 6 7 5]); ngend6D = permute(SONICtable.ngend(:,:,:,:,1,:,:),[1 2 3 4 6 7 5]);       % Note: Compartment 1 has full sonophore coverage -> don't use SONIC-xfs tables!
cfit6D = permute(SONICtable.cfit(:,:,:,:,1,:,:),[1 2 3 4 6 7 5]);

resh = [numel(QmRange),numel(USPaRange),numel(USfreqRange),numel(aBLSRange),repmat(numel(DeltaQmRange),[1,NFS]),repmat(numel(psiQRange),[1,NFS-1])];
reshC = num2cell(resh);
interpgrid = horzcat({QmRange},{USPaRange},{USfreqRange},{aBLSRange},repmat({DeltaQmRange},[1,NFS]),repmat({psiQRange},[1,NFS-1]));

Veff2NFSp3D = reshape(Veff6D,resh); Zeff2NFSp3D = reshape(Zeff6D,resh); Cmeff2NFSp3D = reshape(Cmeff6D,resh);
ngend2NFSp3D = reshape(ngend6D,resh); cfit2NFSp3D = reshape(cfit6D,resh);
cfit2NFSp3Dcell = permute(mat2cell(cell2mat(cellfun(@(X)permute(X',[2*NFS+4,(2:2*NFS+3),1]),cfit2NFSp3D,'UniformOutput',0)),reshC{:},ones(2*NFS+1,1)),[2*NFS+4,(1:2*NFS+3)]);

% rate 6D sonic tables
SONICfields = fieldnames(SONICtable);
SONICrates = sort(SONICfields(cellfun(@(X) contains(X,'a_')|contains(X,'apb_'),SONICfields)));
SONICgates = cellfun(@(X) X(3:end),SONICrates(cellfun(@(X) contains(X,'a_'),SONICrates)),'UniformOutput',0); 
rt = struct;
for i = 1:length(SONICrates)
rt.(SONICrates{i}) = min(SONICtable.(SONICrates{i}),maxRate);
rt.(SONICrates{i}) = reshape(permute(rt.(SONICrates{i})(:,:,:,:,1,:,:),[1 2 3 4 6 7 5]),resh);
end
tempf6Veff = @(queryC) interpn(interpgrid{:},Veff2NFSp3D,queryC{:},'linear');
tempf6Zeff = @(queryC) interpn(interpgrid{:},Zeff2NFSp3D,queryC{:},'linear');
tempf6Cmeff = @(queryC) interpn(interpgrid{:},Cmeff2NFSp3D,queryC{:},'linear');
tempf6ngend = @(queryC) interpn(interpgrid{:},ngend2NFSp3D,queryC{:},'linear');
tempf6cfit = @(queryC) cellfun(@(X) interpn(interpgrid{:},X,queryC{:},'linear'),cfit2NFSp3Dcell); 

f6Veff = @(Qm,USPa,USfreq,aBLS,DeltaQm,psiQ) tempf6Veff(num2cell(vertcat(Qm,USPa,USfreq,aBLS,DeltaQm,psiQ))); % DeltaQm/psiQ are NFSx1 and (NFS-1)x1 column vectors
f6Zeff = @(Qm,USPa,USfreq,aBLS,DeltaQm,psiQ) tempf6Zeff(num2cell(vertcat(Qm,USPa,USfreq,aBLS,DeltaQm,psiQ)));
f6Cmeff = @(Qm,USPa,USfreq,aBLS,DeltaQm,psiQ) tempf6Cmeff(num2cell(vertcat(Qm,USPa,USfreq,aBLS,DeltaQm,psiQ)));
f6ngend= @(Qm,USPa,USfreq,aBLS,DeltaQm,psiQ) tempf6ngend(num2cell(vertcat(Qm,USPa,USfreq,aBLS,DeltaQm,psiQ)));  
f6cfit = @(Qm,USPa,USfreq,aBLS,DeltaQm,psiQ) tempf6cfit(num2cell(vertcat(Qm,USPa,USfreq,aBLS,DeltaQm,psiQ)));

f5Veff = @(Qm,USPa,USfreq,DeltaQm,psiQ) f6Veff(Qm,USPa,USfreq,a,DeltaQm,psiQ); 
f5Zeff = @(Qm,USPa,USfreq,DeltaQm,psiQ) f6Zeff(Qm,USPa,USfreq,a,DeltaQm,psiQ); %#ok<*NASGU>
f5Cmeff = @(Qm,USPa,USfreq,DeltaQm,psiQ) f6Cmeff(Qm,USPa,USfreq,a,DeltaQm,psiQ);
f5ngend = @(Qm,USPa,USfreq,DeltaQm,psiQ) f6ngend(Qm,USPa,USfreq,a,DeltaQm,psiQ);
f5cfit = @(Qm,USPa,USfreq,DeltaQm,psiQ) f6cfit(Qm,USPa,USfreq,a,DeltaQm,psiQ);

f6rt = struct; f5rt = struct; tempf6rt = struct;
VecVeffPa = zeros(1,length(QmRange));
VecrtPa = struct; f1rt0 = struct; f1rtPa = struct;
for i = 1:length(SONICrates)
tempf6rt.(SONICrates{i}) = @(queryC) interpn(interpgrid{:},rt.(SONICrates{i}),queryC{:},'linear');
f6rt.(SONICrates{i}) =  @(Qm,USPa,USfreq,aBLS,DeltaQm,psiQ) tempf6rt.(SONICrates{i})(num2cell(vertcat(Qm,USPa,USfreq,aBLS,DeltaQm,psiQ)));   
f5rt.(SONICrates{i}) = @(Qm,USPa,USfreq,DeltaQm,psiQ) f6rt.(SONICrates{i})(Qm,USPa,USfreq,a,DeltaQm,psiQ); 
end
fVCa = @(cCai) 10^(3)*((Rg*Temp)/(2*Far))*log(cCae./cCai); % Nernst equation for Ca-potential [mV] (if not assumed constant)

% 3. Initial conditions and timespan
% (order:) Y = [Q,SONICgates [ alphabetically ] , other gates [r0,d20,w0], wLock, hProtein, cCai]  (other gates are cCai/hprotein dependent)
Qm0 = Cm0*(10^(-3)*Vm0);            % Approximation of rest charge
Y0 = zeros(round(length(SONICrates)/2)+1,1);
Y0(1) = Qm0;
for i = 1:length(SONICgates)
Y0(i+1) = f1rtV.(['a_' SONICgates{i}])(Vm0)./f1rtV.(['apb_' SONICgates{i}])(Vm0);      
end

if MODEL == 4 || MODEL == 5
cCai0 = 240*10^(-6);             % Rest concentration Ca (mol/m^3)
Y0 = horzcat(Y0,cCai0);
elseif MODEL == 9
cCai0 = 5*10^(-6);
d20 = d2inf(cCai0);
r0 = rinf(cCai0);
Y0 = vertcat(Y0,[r0,d20,cCai0]');
elseif MODEL == 11 || MODEL == 12
CA0 = 4.32e-07;           % Rest concentration Ca (uM = umol/l)
Y0 = horzcat(Y0,CA0);
end
if MODEL == 4
% dw/dt=alpha(V)*(1-w-wLock)-beta(V)*w = (winf*(1-wLock)-w)/tauw
% w0 = winf*(1-wLock0) = winf*(1-((k3*P0)/k4)w0)
% w0*(1+winf*(k3*P0)/k4) = winf -> w0 = winf/(1+winf*(k3*P0)/k4)
% dP/dt = k1*(1-P)*cCa^2-k2*P
% dwLock/dt = k3*w*P-k4*wLock
hProtein0 = (k1*cCai0^2)/(k2+k1*cCai0^2); % Intitial condition for w-gate protein (-)
winf_in_Vm0 = f1rtV.('a_w')(Vm0)./f1rtV.('apb_w')(Vm0);
w0 = winf_in_Vm0/((1+winf_in_Vm0)*((k3*hProtein0)/k4));
wLock0 = (k3/k4)*(w0*hProtein0);   % Initial condition for locked w-gate
Y0 = vertcat(Y0,[w0,wLock0,hProtein0]');
end

index = @(A,ind) A(ind);
fc2DeltaQm = @(X) sqrt(sum(reshape(X,2,[]).^2))';            % Fourier components to DeltaQm
fc2phiQ = @(X) index(-atan2(X(~~repmat([0;1],NFS,1)),X(~~repmat([1;0],NFS,1))),(2:NFS));   % Fourier components to psiQ
FC0 = zeros(2*NFS,1);           % Initial fourier components [A1;B1;A2;B2;...] of the bilayer sonophore
flipFC = reshape(flip(reshape((1:2*NFS),2,[])),[],1);  

fQosc = @(Qm,USPa,X) X-10^(-3)*repmat([1;-1],NFS,1).*(index(index(f5cfit(Qm,USPa,USfreq,fc2DeltaQm(X),fc2phiQ(X)),(2:2*NFS+1)),flipFC)+1000*(fBLS/(1-fBLS))*X(flipFC)./Cm0)./...
    (cumsum(repmat([1;0],NFS,1))*(2*pi*USfreq)*RSI*pi*aBLS^2);              % Charge oscillation equation       (Note: second term with + sign, because FCs_{compartment2} ~ -FCs_{compartment1}
tPeriod = (0:0.025/USfreq:1/USfreq)'; 
for j = 1:length(SONICrates)
f1rtVQosc.(SONICrates{j}) = @(Veff,FC)  (1/(tPeriod(end)-tPeriod(1)))*trapz(tPeriod,f1rtV.(SONICrates{j})(1000*((Cm0*1e-3*Veff).*ones(size(tPeriod))+cos(2*pi*USfreq*bsxfun(@times,(1:1:numel(fc2DeltaQm(FC))),tPeriod)+[0,fc2phiQ(FC)'])*fc2DeltaQm(FC))./Cm0));
end

% 4. Solver
global reverseStr;  %#ok<*TLEV>
if DISPLAY || DISPLAYbase, fprintf('\nSolver started\n'); end
SearchRange = [USibegin USiend]; PrecisionCheck = 0;
if SearchMode % SearchMode = 1
IIpa = (SearchRange(1)+SearchRange(2))/2; % W/m^2
else % SearchMode = 0
if ~any(SearchRange)
IIpa = 1000;          
else
IIpa = SearchRange(~~SearchRange)*10^(3-2*find(SearchRange));
end
end
while ~SearchMode || ~PrecisionCheck 
    reverseStr = '';
if MODE == 1
USPa = sqrt(2*rhol*c*IIpa);
USPaT = @(t) USPa*USstep(t);
elseif MODE == 2
USPa = sqrt(2*rhol*c*USipa);
USPaT = @ (t) USPa*USstep(t); 
end

Y0f = [Y0;Y0]; FCm = [];
TvaluesY = [0]; Y = Y0f'; %#ok<NBRAK>
tCURRENT=0;
DispPart = 0;
DispNrPart=double(USpstart>=Tsim||(USpstart+USpd)<=0||(USpstart==0&&USpd==Tsim))+...
   2*double(USpstart<=0&&((USpstart+USpd)>0)&&((USpstart+USpd)<Tsim))+...
   2*double(((USpstart+USpd)>=Tsim)&&(USpstart>0&&USpstart<Tsim))+...
   3*double(USpstart>0&&((USpstart+USpd)<Tsim));

while tCURRENT < Tsim
    DispPart = DispPart+1;
    if DISPLAYbase
    disp(' ');
    disp(['Part ' num2str(DispPart) ' of ' num2str(DispNrPart)]);
    disp(' ');
    end
    if tCURRENT < USpstart 
        ZONE = 1;
        tSONICc = [tCURRENT, min(USpstart,Tsim)]; % Complete (c) tSONIC interval
    elseif tCURRENT >= USpstart && tCURRENT < USpstart+USpd
        ZONE = 2;
        tSONICc = [tCURRENT, min(USpstart+USpd,Tsim)];
    elseif tCURRENT >= USpstart+USpd
        ZONE = 1;
        tSONICc = [tCURRENT,Tsim];
    end
    
if ZONE == 2
nUP = ceil((tSONICc(2)-tSONICc(1))/Tupdate); 
reverseStr2 = '';
else
nUP = 1; 
reverstreStr2 = ''; 
end


for iUP = 1:nUP
if ZONE == 2 && DISPLAYbase == 1
Progress = 100*iUP/nUP;  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr2, msg]);
reverseStr2 = repmat(sprintf('\b'), 1, length(msg));  
end
if ZONE == 2    
tspan = tSONICc(1)+[(iUP-1)*Tupdate iUP*Tupdate];
if iUP==nUP
    tspan(2) = tSONICc(2);
end
FC = fminsearch(@(X) sum(fQosc(Qm0,USPa,X).^2),FC0);
for k=1:length(SONICrates)
f1rtVp.(SONICrates{k}) = @(V) f1rtVQosc.(SONICrates{k})(V,-(fBLS/(1-fBLS))*FC);         % f1rtV at the proteins with oscillations
end
else 
FC=FC0;
tspan = tSONICc;
f1rtVp = f1rtV;
end

if pretabulSONIC
for i = 1:length(QmRange)
VecVeffPa(i) = f5Veff(QmRange(i),USPa,USfreq,fc2DeltaQm(FC),fc2phiQ(FC));
end
for j = 1:length(SONICrates)
VecrtPa.(SONICrates{j}) = zeros(1,length(QmRange));
for i = 1:length(QmRange)
VecrtPa.(SONICrates{j})(i) = f5rt.(SONICrates{j})(QmRange(i),USPa,USfreq,fc2DeltaQm(FC),fc2phiQ(FC));
end
f1rt0.(SONICrates{j}) = @(Q) f1rtV.(SONICrates{j})(1000*Q/Cm0);
f1rtPa.(SONICrates{j}) = @(Q) nakeinterp1(QmRange',VecrtPa.(SONICrates{j}),Q);
end
f1Veff0 = @(Q) 1000*Q/Cm0;
f1VeffPa = @(Q) nakeinterp1(QmRange',VecVeffPa,Q);
else          % No pretabulation
f1Veff0 = @(Q) 1000*Q/Cm0;
f1VeffPa = @(Q) f5Veff(Q,USPa,USfreq,fc2DeltaQm(FC),fc2phiQ(FC));
for j = 1:length(SONICrates)
f1rt0.(SONICrates{j}) = @(Q) f1rtV.(SONICrates{j})(1000*Q/Cm0);
f1rtPa.(SONICrates{j}) = @(Q) f5rt.(SONICrates{j})(Q,USPa,USfreq,fc2DeltaQm(FC),fc2phiQ(FC));
end
end

OdeOpts=odeset('MaxStep',dt,'AbsTol',atol,'RelTol',rtol); 
% The SONIC_MODELNAME_nanoMC files are written to be simultaneously compatible with
% SONICrun_nanoMC and SONICrun_nanoMC_Qosc 
%--------------------------------------------------------------------------
% ---------------REGULAR OR FAST SPIKING NEURONS---------------------------
%--------------------------------------------------------------------------
if MODEL == 1 || MODEL == 2 || MODEL == 6 || MODEL == 7
[t,U] = ode15s(@(t,U) SONIC_RSFS_nanoMC(ESi,USPaT,DISPLAY,tspan,t,U(1),U(2),U(3),U(4),U(5),...
    U(6),U(7),U(8),U(9),U(10),Gna,Vna,Gk,Vk,Gm,Gl,Vl,f1Veff0,f1VeffPa,f1rt0,f1rtPa,f1rtVp,SONICgates,...
    Cm0,a,fBLS,RSI,proteinMode,gateMultip),tspan,Y0f,OdeOpts);
%--------------------------------------------------------------------------
%-------------------LOW THRESHOLD SPIKING NEURONS--------------------------
%--------------------------------------------------------------------------
elseif MODEL == 3 || MODEL == 8
[t,U] = ode23s(@(t,U) SONIC_LTS_nanoMC(ESi,USPaT,DISPLAY,tspan,t,U(1),U(2),U(3),U(4),U(5),U(6),U(7),...
    U(8),U(9),U(10),U(11),U(12),U(13),U(14),Gna,Vna,Gk,Vk,Gm,GT,VCa,Gl,Vl,f1Veff0,f1VeffPa,f1rt0,f1rtPa,f1rtVp,SONICgates,...
    Cm0,a,fBLS,RSI,proteinMode,gateMultip),tspan,Y0f,OdeOpts);
%--------------------------------------------------------------------------
%-------------------THALAMOCORTICAL NEURONS--------------------------------
%--------------------------------------------------------------------------
elseif MODEL == 4
[t,U] = ode23s(@(t,U) SONIC_TC_nanoMC(ESi,USPaT,DISPLAY,tspan,t,U(1),U(2),U(3),U(4),U(5),U(6),U(7),U(8),U(9),U(10),...
    U(11),U(12),U(13),U(14),U(15),U(16),U(17),U(18),U(19),U(20),Gna,Vna,Gk,Vk,GT,fVCa,Gl,Vl,GKL,Gh,ginc,...
    Vh,k1,k2,k3,k4,Far,deffCa,tauCa,f1Veff0,f1VeffPa,f1rt0,f1rtPa,f1rtVp,SONICgates,...
    Cm0,a,fBLS,RSI,proteinMode,gateMultip),tspan,Y0f,OdeOpts);
%--------------------------------------------------------------------------
%------------------NUCLEUS RETICULARIS NEURONS-----------------------------
%--------------------------------------------------------------------------
elseif MODEL == 5
[t,U] = ode23s(@(t,U) SONIC_RE_nanoMC(ESi,USPaT,DISPLAY,tspan,t,U(1),U(2),U(3),U(4),U(5),U(6),U(7),...
    U(8),U(9),U(10),U(11),U(12),U(13),U(14),Gna,Vna,Gk,Vk,GT,fVCa,Gl,Vl,Far,deffCa,tauCa,f1Veff0,f1VeffPa,f1rt0,f1rtPa,f1rtVp,SONICgates,...
    Cm0,a,fBLS,RSI,proteinMode,gateMultip),tspan,Y0f,OdeOpts);
% -------------------------------------------------------------------------
% -----------------------SUBTHALAMIC NUCLEUS MODEL-------------------------
% -------------------------------------------------------------------------
elseif MODEL == 9
[t,U] = ode23s(@(t,U) SONIC_STN_nanoMC(ESi,USPaT,DISPLAY,tspan,t,...
    U(1),U(2),U(3),U(4),U(5),U(6),U(7),U(8),U(9),U(10),U(11),U(12),U(13),...
    U(14),U(15),U(16),U(17),U(18),U(19),U(20),U(21),U(22),U(23),U(24),U(25),U(26),...
    Gna,Vna,Gk,Vk,Gl,Vl,GT,fVCa,GCa,GA,GL,Far,tauCa,f1Veff0,f1VeffPa,f1rt0,f1rtPa,...
    rinf,d2inf,taur,taud2,f1rtVp,SONICgates,Cm0,a,fBLS,RSI,proteinMode,gateMultip),tspan,Y0f,OdeOpts);
% -------------------------------------------------------------------------
% -----------------------THALAMUS RUBIN-TERMAN MODEL-----------------------
% -------------------------------------------------------------------------
elseif MODEL == 10
[t,U] = ode23s(@(t,U) SONIC_ThRT_nanoMC(ESi,USPaT,DISPLAY,tspan,t,U(1),U(2),U(3),...
    U(4),U(5),U(6),Gna,Vna,Gk,Vk,Gl,Vl,GT,VT,minf,pinf,f1Veff0,f1VeffPa,f1rt0,f1rtPa,f1rtVp,SONICgates,...
    Cm0,a,fBLS,RSI,proteinMode,gateMultip),tspan,Y0f,OdeOpts);
% -------------------------------------------------------------------------
% ----------------------GLOBUS PALLIDUS INTERNUS NUCLEUS ------------------
% -------------------------------------------------------------------------
elseif MODEL == 11
[t,U] = ode23s(@(t,U) SONIC_GPi_nanoMC(ESi,USPaT,DISPLAY,tspan,t,U(1),U(2),U(3),U(4),U(5),...
    U(6),U(7),U(8),U(9),U(10),Gna,Vna,Gk,Vk,Gl,Vl,GT,VT,GCa,VCa,Gahp,minf,ainf,sinf,f1Veff0,f1VeffPa,f1rt0,f1rtPa,f1rtVp,SONICgates,...
    Cm0,a,fBLS,RSI,proteinMode,gateMultip),tspan,Y0f,OdeOpts);
% -------------------------------------------------------------------------
% ----------------------GLOBUS PALLIDUS EXTERNUS NUCLEUS-------------------
% -------------------------------------------------------------------------
elseif MODEL == 12
[t,U] = ode23s(@(t,U) SONIC_GPe_nanoMC(ESi,USPaT,DISPLAY,tspan,t,U(1),U(2),U(3),U(4),U(5),...
    U(6),U(7),U(8),U(9),U(10),Gna,Vna,Gk,Vk,Gl,Vl,GT,VT,GCa,VCa,Gahp,minf,ainf,sinf,f1Veff0,f1VeffPa,f1rt0,f1rtPa,f1rtVp,SONICgates,...
    Cm0,a,fBLS,RSI,proteinMode,gateMultip),tspan,Y0f,OdeOpts);
% -------------------------------------------------------------------------
% -----------------------MEDIUM SPINY STRIATUM NEURONS---------------------
% -------------------------------------------------------------------------
elseif MODEL == 13
[t,U] = ode23s(@(t,U) SONIC_MSN_nanoMC(ESi,USPaT,DISPLAY,tspan,t,U(1),U(2),U(3),U(4),U(5),...
    U(6),U(7),U(8),U(9),U(10),Gna,Vna,Gk,Vk,Gl,Vl,Vm,f1Veff0,f1VeffPa,f1rt0,f1rtPa,f1rtVp,SONICgates,...
    Cm0,a,fBLS,RSI,proteinMode,gateMultip),tspan,Y0f,OdeOpts);
% -------------------------------------------------------------------------
% --------------------------HODGKIN-HUXLEY NEURONS-------------------------
% -------------------------------------------------------------------------
elseif MODEL == 14
[t,U] = ode23s(@(t,U) SONIC_HH_nanoMC(ESi,USPaT,DISPLAY,tspan,t,U(1),U(2),U(3),U(4),...
    U(5),U(6),U(7),U(8),Gna,Vna,Gk,Vk,Gl,Vl,f1Veff0,f1VeffPa,f1rt0,f1rtPa,f1rtVp,SONICgates,...
    Cm0,a,fBLS,RSI,proteinMode,gateMultip),tspan,Y0f,OdeOpts);
end
TvaluesY = vertcat(TvaluesY,t); Y = vertcat(Y,U); FC0 = FC; Qm0 = Y(end,1); Y0f = Y(end,:); FCm = horzcat(FCm,FC); %#ok<AGROW>
tCURRENT=t(end);
end
end
Y1 = Y(:,1:end/2); Y2 = Y(:,end/2+1:end);

APindex1 = (10^5*Y1(:,1)>Qthresh)&(circshift(10^5*Y1(:,1),1,1)<Qthresh);
APindex1(1) = 0; % Remove circshift artifact
APtimes1 = TvaluesY(APindex1);  % AP times [s]
clear APindex1; % Save all memory that can be saved...
if (threshMode == 0)
NeuronActivated1 = ~isempty(APtimes1); % Bool: 1 if neuron is activated
elseif (threshMode == 1)
NeuronActivated1 = ~isempty(APtimes1(APtimes1>=USpstart&APtimes1<=(USpstart+USpd)));
end
APindex2 = (10^5*Y2(:,1)>Qthresh)&(circshift(10^5*Y2(:,1),1,1)<Qthresh);
APindex2(1) = 0; % Remove circshift artifact
APtimes2 = TvaluesY(APindex2);  % AP times [s]
clear APindex2; % Save all memory that can be saved...
if (threshMode == 0)
NeuronActivated2 = ~isempty(APtimes2); % Bool: 1 if neuron is activated
elseif (threshMode == 1)
NeuronActivated2 = ~isempty(APtimes2(APtimes2>=USpstart&APtimes2<=(USpstart+USpd)));
end
NeuronActivated = NeuronActivated1||NeuronActivated2;
APtimes = {APtimes1';APtimes2'};

if MODE == 2
SaveStr=['nanoMC-APtimes(' modelName ')-Tsim=' num2str(Tsim) '-US(' num2str(USpstart) ',' num2str(USpd) ',' ...
        num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' USisppa ...
        ')-ES(' num2str(ESpstart) ',' num2str(ESpd) ',' num2str(ESdc) ',' ...
        num2str(ESprf) ',' ESisppa ')-aBLS=(' num2str(aBLS) ')-fBLS=(' num2str(fBLS) ')' modeStr '.mat'];
save(SaveStr,'APtimes');

SaveStr2=['nanoMC-Chargevt(' modelName ')-Tsim=' num2str(Tsim) '-US(' num2str(USpstart) ',' num2str(USpd) ',' ...
        num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' USisppa ...
        ')-ES(' num2str(ESpstart) ',' num2str(ESpd) ',' num2str(ESdc) ',' ...
        num2str(ESprf) ',' ESisppa ')-aBLS=(' num2str(aBLS) ')-fBLS=(' num2str(fBLS) ')' modeStr '.mat'];
saveChargeSample = [TvaluesY, 10^5*Y1(:,1),10^5*Y2(:,1)];
save(SaveStr2,'saveChargeSample');
break;
elseif MODE == 1
SearchRange(NeuronActivated+1) = IIpa;
if SearchMode % SearchMode = 1
IIpa = (SearchRange(1)+SearchRange(2))/2;
PrecisionCheck=((SearchRange(2)-SearchRange(1))<...
    10^(floor(log10(SearchRange(2)))-SearchPrecision));
else % SearchMode = 0
if ~any(~SearchRange)
SearchMode = 1;
IIpa = (SearchRange(1)+SearchRange(2))/2;
PrecisionCheck=((SearchRange(2)-SearchRange(1))<...
    10^(floor(log10(SearchRange(2)))-SearchPrecision));
else
IIpa = IIpa*10^(3-2*find(SearchRange));
end
end
if DISPLAY
disp(' '); %#ok<*UNRCH>
Checkpoint=['nanoMC-CP(' modelName ')-Tsim=' num2str(Tsim) '-US(' num2str(USpstart) ',' num2str(USpd) ',' ...
        num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' USisppa ...
        ')-ES(' num2str(ESpstart) ',' num2str(ESpd) ',' num2str(ESdc) ',' ...
        num2str(ESprf) ',' ESisppa '):' num2str(SearchRange) '---aBLS=(' num2str(aBLS) ')-fBLS=(' num2str(fBLS) ')' modeStr];
disp(Checkpoint);   
fprintf('\n'); 
end
end
end
if MODE == 1
SaveStr=['nanoMC-Thresh(' modelName ')-Tsim=' num2str(Tsim) '-US(' num2str(USpstart) ','...
    num2str(USpd) ',' num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' ...
    USisppa ')-ES(' num2str(ESpstart) ',' num2str(ESpd) ',' num2str(ESdc) ',' ...
        num2str(ESprf) ',' ESisppa ')-aBLS=(' num2str(aBLS) ')-fBLS=(' num2str(fBLS) ')' modeStr '.mat'];    
save(SaveStr,'IIpa');
end
TTime = toc;
if DISPLAY
disp(' ');
disp(['Program finished in ' num2str(round(TTime,1)) 's']);
fprintf('Post-processing...\n'); 
end

% 5. Results
if PLOT
Charge1 = 10^5*Y1(:,1); Charge2 = 10^5*Y2(:,1); % Charge [nC/cm^2]
Charge = {Charge1;Charge2};
Chargetot = (Charge1*fBLS+Charge2*(1-fBLS));
tline = (0:0.025/USfreq:Tsim)'; QoscBLS = zeros(size(tline));
for i = 1:numel(tline)
if ~(tline(i) <=USpstart || tline(i) >=USpstart+USpd)  
j = ceil((tline(i)-USpstart)/Tupdate);
QoscBLS(i) = cos(2*pi*USfreq*tline(i)+[0,fc2phiQ(FCm(:,j))'])*fc2DeltaQm(FCm(:,j));
end
end
QoscP = -(fBLS/(1-fBLS))*QoscBLS; Qosc = {QoscBLS;QoscP};

if MODEL == 4
SONICgatesN = vertcat(SONICgates,{'w0','wLock','hProtein'}');
elseif MODEL == 9
SONICgatesN = vertcat(SONICgates,{'r0','d20'}');    
end
if any(MODEL==[4,5,9,11,12])
SONICgatesN = vertcat(SONICgatesN,'cCai');
end
if ~exist('SONICgatesN','var'), SONICgatesN = SONICgates; end
if PLOT == 2
SaveDataStr=['Qosc_nanoMC-Data(' modelName ')-Tsim=' num2str(Tsim) '-US(' num2str(USpstart) ','...
num2str(USpd) ',' num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' ...
USisppa ')-ES(' num2str(ESpstart) ',' num2str(ESpd) ',' num2str(ESdc) ',' ...
    num2str(ESprf) ',' ESisppa ')--(aBLS,fBLS)=(' num2str(aBLS) ',' num2str(fBLS) ')' modeStr '.mat'];
TvaluesYms = 10^(3)*TvaluesY'; % [ms]
saveData.TvaluesYms = TvaluesYms; saveData.Charge1 = Charge1; saveData.Charge2=Charge2; saveData.Chargetot = Chargetot;
saveData.Y1 = Y1; saveData.Y2 = Y2;
saveData.tline = tline; saveData.Qosc = Qosc; saveData.FCm = FCm;
saveData.f5Veff = f5Veff; saveData.f5Zeff = f5Zeff;
saveData.f5Cmeff = f5Cmeff; saveData.f5ngend = f5ngend;
saveData.f5cfit = f5cfit;
save(SaveDataStr,'saveData','-v7.3');
end
if PLOT == 1
figure;
nrGates = size(Y1,2)-1; 
TvaluesYms = 10^(3)*TvaluesY'; % [ms]

% VeffSample = {zeros(size(TvaluesY)),10^(-2)*Charge2./Cm0};   % (mV)
% for i = 1:length(TvaluesY)
% if USstep(TvaluesY(i)) == 0
% VeffSample{1}(i) = f1Veff0(10^(-5)*Charge{1}(i));
% else
% VeffSample{1}(i) = f1VeffPa(10^(-5)*Charge{1}(i));
% end
% end

subplot(1+ceil(nrGates/2),1,1);
hold on;
yyaxis left;
plot(TvaluesYms,Charge1,'color','b','linestyle','-','linewidth',2);
plot(TvaluesYms,Charge2,'color','b','linestyle','--','linewidth',2);
plot(TvaluesYms,Chargetot,'color','k','linestyle','-','linewidth',2);
ylim([-100,100]);
ylabel('Charge [nC/cm^2]');
% yyaxis right;
% plot(TvaluesYms,VeffSample{1},'color','r','linestyle','-');
% plot(TvaluesYms,VeffSample{2},'color','r','linestyle','--');
% ylim([-100,100]);
% ylabel('V_{eff} [mV]');
hold off;
legend({'Bilayer sonophore','Protein islands'});
for i=1:floor(nrGates/2)
subplot(1+ceil(nrGates/2),1,i+1);
hold on;
yyaxis left;
if proteinMode == 0
plot(TvaluesYms,Y1(:,2*i),'b','linestyle','-');
end
plot(TvaluesYms,Y2(:,2*i),'b','linestyle','--'); 
ylim([0,1]);
ylabel(SONICgatesN{2*i-1});
yyaxis right;
if proteinMode == 0
plot(TvaluesYms,Y1(:,2*i+1),'r','linestyle','-');
end
plot(TvaluesYms,Y2(:,2*i+1),'r','linestyle','--');
ylim([0,1]);
ylabel(SONICgatesN{2*i});
hold off;
set(gcf,'color','w');
end
if mod(nrGates,2)==1
subplot(1+ceil(nrGates/2),1,1+ceil(nrGates/2));
hold on;
if proteinMode == 0
plot(TvaluesYms,Y1(:,end),'b','linestyle','-');
end
plot(TvaluesYms,Y2(:,end),'b','linestyle','--');
ylim([0,1]);
ylabel(SONICgatesN{nrGates});
hold off;
end
xlabel('Time [ms]');
set(findall(gcf,'-property','FontSize'),'FontSize',14)
end
end