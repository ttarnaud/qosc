function SONICrun(Tsim,MODE,USpstart,USpd,USfreq,USdc,USprf,USisppa,ESpstart,ESpd,...
    ESdc,ESprf,ESisppa,PLOT,model,USibegin,USiend,SearchMode)
coder.extrinsic('nakeinterp1');
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
% 3. Units: SI-units are used, except for the membrane voltage [mV].
% -> All input variables are strings, because this works more fluently
% -PBS files for HPC simulations:
Tsim = str2double(Tsim); MODE = str2double(MODE); USpstart = str2double(USpstart);
USpd = str2double(USpd); USfreq = str2double(USfreq); USdc= str2double(USdc);
USprf = str2double(USprf); USipa = str2double(USisppa); ESpstart = str2double(ESpstart);
ESpd = str2double(ESpd); ESdc = str2double(ESdc); ESprf = str2double(ESprf);
ESipa = str2double(ESisppa); PLOT = str2double(PLOT); 
USibegin = str2double(USibegin); USiend = str2double(USiend);
SearchMode = str2double(SearchMode);

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

% Piezoelectric solver (see Plaksin et al.,2014; Plaksin et al.,2016;
% Krasovitski et al., 2011 and PhD thesis Plaksin, 2016).

% Remark:
% The constant ks (areal modulus) is taken as the ratio of the total
% tension to the areal strain. The tension T of one leaflet is taken as
% half the total tension.

% Remark 2: Memory considerations
% The piezoelectric solver implemented here applies a VSVO-solver (ode113)
% As a consequence the solution can not be preallocated.
% Memory can however be used more efficiently, by disposing unnecessary
% solutions.
MemorySaveMode = 3;
if MemorySaveMode == 0
    disp('Note: memory-save-mode is off');
    disp(' ');
elseif MemorySaveMode == 1
    disp('Note: memory-save-mode is on (1)');
    disp(' ');
elseif MemorySaveMode == 2
    disp('Note: memory-save-mode is on (2)');
    disp(' ');
elseif MemorySaveMode == 3
    disp('Note: memory-save-mode is on (3)');
    disp(' ');
end
% If 1, the X matrix holding Z,dZ,na,Pin and the concentration of air will
% be cleared each iteration. This saves memory but as a result no plots
% will be available of X (including no capacitance and voltage plots).
% If 2, everything in 1 applies + the Z values are interpolated to the
% TvaluesX axis and stored to allow potential plots
% If 3, only data necessary to plot Q values will be stored (low memory option)

% Remark 3: Speed considerations
% SpeedUp -> NICE-model, UpSpeed -> BLS-model
SpeedUp = 5;
UpSpeed = 0;
CORRTres = 0.99;     % Threshold for normalized unbiased autocorrelation in
% periodicity analysis 
% For long simulation times Tsim, pulse durations PD and short minimum step times
% dt, the default program (SpeedUp=0) is rather slow (f.i.
% Tsim=50ms,PD=40ms,dt=0.025/freq takes about 15 hours for one simulation
% configuration). 
% SpeedUp = 1 speeds the solver up, by skipping periodic calculations of X
% SpeedUp = 2 speeds the solver up, by using lookup tables with linear interpolation 
% SpeedUp = 3 speeds the solver up, by using lookup tables with the interp1qr accelerated function
% SpeedUp = 4 speeds the solver up, by using C++ accelerated (.MEX) nakeinterp1 interpolation
% SpeedUp = 5 speeds the solver up, by additionally exploiting periodicity before interpolation
% SpeedUp = 6 -> Same as SpeedUp = 5 but with interp1qr (slower than SpeedUp = 5)
% Note: in SpeedUp == 3 || 4, interp1qr||C++ interpolations will also be applied on all bottlenecks, identified by the profiler.
% Note2: interp1qr-linear is faster than interp1-NearestNeighbour interpolation
% UpSpeed = 1 speeds the solver up by lookup tables for the intramolecular integral
% UpSpeed = 2 speeds the solver up by lookup tables for the intramembrane volume
% UpSpeed  = 3 speeds the solver up by lookup tables for the membrane capacitance

DISPLAY = 1;
% Display level. Note: higher display level will give more runtime information but will slow the program 
% DISPLAY = 0 -> No information displayed (use this option for HPC simulations)
% DISPLAY = 1 -> Display progress based on update nr. alone
% DISPLAY = 2 -> Display progress inside update nr.

ERROR = 1; % If 1, show error for unphysical situation Z<=-delta/2, if 0 extrapolate value for Pa-r
% (Use 0 for HPC simulations) 
if ERROR == 0
    disp('Note 2: solver will not check for unphysical molecular pressure');
    disp(' ');
elseif ERROR == 1
    disp('Note 2: solver will check for unphysical molecular pressure');
    disp(' ');
end
Charges = 0; % This mode is meant for investigation on the mechanism of US-STN stimulation:
% All charges are calculated, but the solver is slowed down
if Charges == 1
    disp('Note 3: Charge-mode is on for STN-mechanism investigation!');
    disp(' ');
end

tic;
DiffusionApproximation = 1; % Bool: if 1, apply approximation (*) on diffusion of air.
%*Krasovitski 2011 doesn't use this approximation, Plaksin 2014 does use it. 

% 1. Parameters
Qthresh = 0;            % Threshold for Q for AP-discrimination [nC/cm^2]
xfs = 1;                % Active area fraction (-)
a=32*10^(-9);			% radius leaflet boundary (m)
ksi = 0.5*10^(-9);		% Boundary layer length leaflet-medium (m)
deltaR = ksi;           % discretisation step outside the BLS (m)
Rsim = 100*ksi;         % Total length simulation domain outside BLS (m)
eps0 = 8.854187817*10^(-12);% absolute vacuumpermittivity (F/m) 
kappa = 1;              % Polytropic constant (-)
epsr = 1;				% Relative permittivity intramembrane cavity (-)
c = 1515;				% Speed of sound surrounding medium (m/s)
rhol = 1028;			% Density surrounding medium (kg/m^3)
Po = 10^5;				% Static pressure surrounding medium (Pa)
USPa = sqrt(2*rhol*c*USipa);
omega = 2*pi*USfreq;  	% Radial frequency (Hz)
ks = 0.24; 				% Areal modulus of bilayer membrane (N/m)
delta0=2*10^(-9);		% Thickness leaflet (m)
mus = 0.035;			% Dynamic viscosity leaflets (Pa*s)
mul = 0.7*10^(-3);		% Dynamic viscosity surrounding medium (Pa*s)	
Cm0 = 0.01;				% Rest capacitance (F/m^2)	 
Rg = 8.314;             % Universal gas constant (J/(K*mol))
Far = 96485.3329;       % Faraday constant (C/mol) 
%Temp = 309.15; 		    % Surrounding medium temperature (K)	
%Temp = 306.15;
Temp = 36+273.15;
Da = 3*10^(-9);			% Diffusion coefficient of air in surrounding medium (m^2/s)
Ca = 0.62;				% Molar air concentration external medium (mol/m^3)
ka = 1.613*10^5;			% Henry constant dissolved air external medium (Pa*m^3/mol) 
Ar = 10^5;				% Attraction/repulsion pressure coefficient (Pa) 
deltax = 1.4*10^(-9); 	% Initial gap between leaflets (no charge) (m) 
x = 5;					% Repulsion exponent
y = 3.3;				% Attraction exponent
Nout = floor(Rsim/deltaR);% Number of steps outside of BLS
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
VT = -56.2;             % Spike treshold adjustment parameter (mV)
taumax = 608*10^(-3);   % Decay time constant slow non-inactivating K+ (s)
Vm0 = -71.9;            % rest membrane potential (mV)
delta = 1.26*10^(-9);   % Initial gap between leaflets (charge) (m)
elseif MODEL == 2       % Fast spiking neuron
Gna = 580;				% Maximal conductance of the Na-channel (S/m^2)
Vna = 50;               % Na nernst potential (mV)
Gk = 39;				% Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Gm = 0.787;				% Maximal conductance of the slow non-inactivating K-channel (S/m^2)
Vk = -90;				% Potassium nernst potential (mV)
Gl = 0.38;				% Maximal conductance of the non-voltage-dependent non-specific ion channel (S/m^2)
Vl = -70.4;				% Leak nernst potential (mV)
VT = -57.9;             % Spike treshold adjustment parameter (mV)
taumax = 502*10^(-3);   % Decay time constant slow non-inactivating K+ (s)
Vm0 = -71.4;            % Rest membrane potential (mV)
delta = 1.26*10^(-9);   % Initial gap between leaflet (charge) (m)  
elseif MODEL == 3       % Low threshold spiking neuron
Gna = 500;				% Maximal conductance of the Na-channel (S/m^2)
Vna = 50;               % Na nernst potential (mV)
Gk = 40;				% Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Gm = 0.28;				% Maximal conductance of the slow non-inactivating K-channel (S/m^2)
Vk = -90;				% Potassium nernst potential (mV)
Gl = 0.19;				% Maximal conductance of the non-voltage-dependent non-specific ion channel (S/m^2)
Vl = -50;				% Leak nernst potential (mV)
VT = -50;               % Spike treshold adjustment parameter (mV)
taumax = 4;             % Decay time constant slow non-inactivating K+ (s)
GT = 4;                 % Maximal conductance of low-threshold Ca2+ channels (S/m^2)
VCa = 120;              % Nernst potential of Ca2+ (mV)
Vx = -7;				% Shift Ca2+ voltage (mV)
Vm0 = -54;              % Rest membrane potential (mV)
delta = 1.3*10^(-9);    % Initial gap between leaflets (charge) (m)
elseif MODEL == 4       % Thalamocortical neuron
Gna = 900;				% Maximal conductance of the Na-channel (S/m^2)
Vna = 50;               % Na nernst potential (mV)
Gk = 100;				% Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Vk = -90;				% Potassium nernst potential (mV)
Gl = 0.1;				% Maximal conductance of the non-voltage-dependent non-specific ion channel (S/m^2)
Vl = -70;				% Leak nernst potential (mV)
VT = -52;               % Spike treshold adjustment parameter (mV)
GT = 20;                % Maximal conductance of low-threshold Ca2+ channels
GKL = 0.138;            % Maximal conductance of leak potassium currents
Gh = 0.175;             % Maximal conductance of hyperpolarization-activated mixed cationic current
ginc = 2;               % Locked gate relative conductance (-)
Vh = -40;               % Reversal potential of a hyperpolarization-activated mixed cationic urrent
Vx = 0;                 % Shift Ca2+ voltage (mV)
Vm0 = -63.4;            % Rest membrane potential (mV)
delta = 1.28*10^(-9);   % Initial gap between leaflet (charge) (m)  
tauCa = 5*10^(-3);      % Calcium decay time constant (s)
elseif MODEL == 5       % Reticular thalamus neuron
Gna = 2000;				% Maximal conductance of the Na-channel (S/m^2)
Vna = 50;               % Na nernst potential (mV)
Gk = 200;				% Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Vk = -90;				% Potassium nernst potential (mV)
Gl = 0.5;				% Maximal conductance of the non-voltage-dependent non-specific ion channel (S/m^2)
Vl = -90;				% Leak nernst potential (mV)
VT = -67;               % Spike treshold adjustment parameter (mV)
GT = 30;                % Maximal conductance of low-threshold Ca2+ channels
Vm0 = -89.5;            % Rest membrane potential (mV)
delta = 1.21*10^(-9);   % Initial gap between leaflet (charge) (m)
tauCa = 5*10^(-3);      % Calcium decay time constant (s)
elseif MODEL == 6       % Regular spiking ferret visual cortex neuron
Gna = 500;				% Maximal conductance of the Na-channel (S/m^2)
Vna = 50;               % Na nernst potential (mV)
Gk = 50;				% Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Vk = -90;				% Potassium nernst potential (mV)
Gm = 0.7;				% Maximal conductance of the slow non-inactivating K-channel (S/m^2
Gl = 1;				    % Maximal conductance of the non-voltage-dependent non-specific ion channel (S/m^2)
Vl = -70;				% Leak nernst potential (mV)
taumax = 1;             % Decay time constant slow non-inactivating K+ (s)
VT = -55;               % Spike treshold adjustment parameter (mV)
Vm0 = -70.4;            % Rest membrane potential (mV)
delta = 1.26*10^(-9);   % Initial gap between leaflet (charge) (m)
elseif MODEL == 7       % Fast spiking ferret visual cortex neuron
Gna = 500;				% Maximal conductance of the Na-channel (S/m^2)
Vna = 50;               % Na nernst potential (mV)
Gk = 100;				% Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Vk = -90;				% Potassium nernst potential (mV)
Gm = 0;			    	% Maximal conductance of the slow non-inactivating K-channel (S/m^2
Gl = 1.5;				% Maximal conductance of the non-voltage-dependent non-specific ion channel (S/m^2)
Vl = -70;				% Leak nernst potential (mV)
taumax = 1;             % Decay time constant slow non-inactivating K+ (s)
VT = -55;               % Spike treshold adjustment parameter (mV)
Vm0 = -70;              % Rest membrane potential (mV)
delta = 1.26*10^(-9);   % Initial gap between leaflet (charge) (m)    
elseif MODEL == 8       % Low threshold spiking cat association cortex neuron
Gna = 500;				% Maximal conductance of the Na-channel (S/m^2)
Vna = 50;               % Na nernst potential (mV)
Gk = 50;				% Maximal conductance of the delayed-rectifier K-channel (S/m^2)
Vk = -90;				% Potassium nernst potential (mV)
Gm = 0.3;				% Maximal conductance of the slow non-inactivating K-channel (S/m^2
Gl = 0.1;				% Maximal conductance of the non-voltage-dependent non-specific ion channel (S/m^2)
Vl = -85;				% Leak nernst potential (mV)
taumax = 1;             % Decay time constant slow non-inactivating K+ (s)
VT = -55;               % Spike treshold adjustment parameter (mV)
GT = 4;                 % Maximal conductance of low-threshold Ca2+ channels
VCa = 120;              % Nernst potential of Ca2+ (mV)
Vx = -2;				% Shift Ca2+ voltage (mV)
Vm0 = -84.6;            % Rest membrane potential (mV)
delta = 1.3*10^(-9);    % Initial gap between leaflet (charge) (m)
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
delta=1.2923e-09;       % Initial gap between leaflets (charge) (m) (calculated from Vm0 defined supra)
tauCa = 0.5*10^(-3);    % Calcium decay time constant (s) 
elseif MODEL == 10      % Thalamic Rubin-Terman based neuron
Gl=0.5;                 % Maximal conductance of non-specific non-voltage dependent ion channel (S/m^2)  
Vl=-70;                 % Leak nernst potential (mV)
Gna=30;                 % Maximal conductance of the Na-channel (S/m^2)
Vna=50;                 % Na nernst potential (mV)
Gk=50;                  % Maximal conductance of the K-channel (S/m^2)
Vk=-75;                 % K nernst potential (mV)
GT=50;                  % Maximal conductance of T-type low threshold Ca2+ channels (S/m^2)
VT=0;                   % Nernst potential of T-type low threshold Ca2+ channel (mV)
Vm0=-65;                % Resting potential (mV)
delta=1.2744e-09;       % Initial gap betweel leaflets (m)
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
delta=1.2713e-09;       % Initial gap between leaflets (charge) (m)
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
delta=1.2713e-09;       % Initial gap between leaflets (charge) (m)
elseif MODEL == 13      % Striatum medium spiny neuron
Gl=1;                   % Maximal conductance of the non-voltage dependent non-specific ion channel (S/m^2)
Vl=-67;                 % Leak nernst potential (mV)
Gna=1000;               % Maximal conductance of the Na-channel (S/m^2)
Vna=50;                 % Na nernst potential (mV)
Gk=800;                 % Maximal conductance of the delayed rectifier K-channel
Vk=-100;                % K delayed rectifier nernst potential (mV)
Vm=-100;                % K non-inactivating nernst potential (mV)
Vm0=-87;                % Resting potential (mV)
delta=1.2177e-09;       % Initial gap between leaflets (charge) (m)
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
delta = 1.26e-9;        % Initial gap between leaflets (charge) (m)
end
global Hmin;
Hmin = 1*10^(-12);     % Linear approximation Rayleigh-Plesset
Hmax = 1*10^(-9);        % Hmax for fzero
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

Ap1 = (Da./(((0:Nout-1)*deltaR+a)*deltaR)).*((1:Nout)+a/deltaR);    % 1 diagonal of A
A0 = -2*(Da./(((1:Nout)*deltaR+a)*deltaR)).*((1:Nout)+a/deltaR); % 0 diagonal of A
Am1 = (Da./(((2:Nout+1)*deltaR+a)*deltaR)).*((1:Nout)+a/deltaR);  % -1 diagonal of A
A = spdiags([Am1' A0' Ap1'],[-1 0 1],Nout,Nout); % A-discretisation matrix for outside air diffusion
B1 = @(Pin) (Da/((a+deltaR)*deltaR))*(a/deltaR)*(Pin/ka);% First B element
Bend = (Da/((Nout*deltaR+a)*deltaR))*(Nout+1+a/deltaR)*Ca;% Final B element
B = @(Pin) [B1(Pin); zeros(Nout-2,1); Bend]; % B-discretisation matrix for outside air diffusion

% 2. Important functions
% 2.1 Rate functions 
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 4 || MODEL == 5 || MODEL == 6 ||...
        MODEL == 7 || MODEL == 8
am = @(V) -1000*0.32*((-4)*double((V-VT-13)==0)+...
    ((V-VT-13)/(double((V-VT-13)~=0)*exp(-((V-VT-13)/4))-1)));% Rate constant alpha_m [1/s]
bm = @(V) 1000*0.28*(5*double((V-VT-40)==0)+...
    (V-VT-40)/(double((V-VT-40)~=0)*exp(((V-VT-40)/5))-1)); % Rate constant beta_m [1/s]
an = @(V) -1000*0.032*((-5)*double((V-VT-15)==0)+...
    (V-VT-15)/(double((V-VT-15)~=0)*exp(-((V-VT-15)/5))-1));% Rate constant alpha_n [1/s]
bn = @(V) 1000*0.5*exp(-(V-VT-10)/40);                    % Rate constant beta_n [1/s]
ah = @(V) 1000*0.128*exp(-((V-VT-17)/18));                % Rate constant alpha_h [1/s]
bh = @(V) (1000*4)/(1+exp(-((V-VT-40)/5)));               % Rate constant beta_h [1/s]
end
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 6 || MODEL == 7 || MODEL == 8 
pinf = @(V) 1/(1+exp(-((V+35)/10)));                      % Rest p-value [-]
taup = @(V) taumax/(3.3*exp((V+35)/20)+exp(-(V+35)/20));  % Time-constant [s] for p 
end
if MODEL == 3 || MODEL == 4 || MODEL == 8
sinf = @(V) 1/(1+exp(-(V+Vx+57)/6.2));					  % Rest s-value [-]
uinf = @(V) 1/(1+exp((V+Vx+81)/4));				          % Rest u-value [-] 
taus = @(V) 10^(-3)*((1/3.7)*(0.612+1/(exp(-((V+Vx+132)/16.7))+exp((V+Vx+16.8)/18.2))));  % Time-constant [s] for s
tauu = @(V) 10^(-3)*funtauu(V,Vx);                        % Time-constant [s] for u
end
if MODEL == 5
sinf = @(V) 1/(1+exp(-(V+52)/7.4));                       % Rest s-value [-]
uinf = @(V) 1/(1+exp((V+80)/5));                          % rest u-value [-]
taus = @(V) 10^(-3)*(1+0.33/(exp((V+27)/10)+exp(-(V+102)/15))); % Time constant [s] for s
tauu = @(V) 10^(-3)*(28.3+0.33/(exp((V+48)/4)+exp(-(V+407)/50))); % Time constant [s] for u
end
if MODEL == 4
winf = @(V) 1/(1+exp((V+75)/5.5));                        % Rest w-value [-]    
tauw = @(V) 10^(-3)/(exp(-14.59-0.086*V)+exp(-1.87+0.0701*V));  % Time-constant [s] for w
end
if MODEL == 9
minf = @(V) 1./(1+exp(-(V+40)./8));            % rest m-value [-]
hinf = @(V) 1./(1+exp((V+45.5)./6.4));          % rest h-value [-]
ninf = @(V) 1./(1+exp(-(V+41)./14));            % rest n-value [-]
pinf = @(V) 1./(1+exp(-(V+56)./6.7));          % rest p-value [-]
qinf = @(V) 1./(1+exp((V+85)./5.8));           % rest q-value [-]
rinf = @(kCA) 1./(1+exp(-(kCA-0.17)./0.08));  % rest r-value [-]
ainf = @(V) 1./(1+exp(-(V+45)./14.7));          % rest a-value [-]
binf = @(V) 1./(1+exp((V+90)./7.5));            % rest b-value [-]
cinf = @(V) 1./(1+exp(-(V+30.6)./5));           % rest c-value [-]
d1inf = @(V) 1./(1+exp((V+60)./7.5));          % rest d1-value [-]
d2inf = @(kCA) 1./(1+exp((kCA-0.1)./0.02));        % rest d2-value [-]

taum = @(V) (10^(-3))*(0.2+3./(1+exp((V+53)/0.7)));   % Time constant [s] for m
tauh = @(V) (10^(-3))*(24.5./(exp((V+50)/15)+exp(-(V+50)/16))); % Time-constant [s] for h
taun = @(V) (10^(-3))*(11./(exp((V+40)/40)+exp(-(V+40)/50))); % Time-constant [s] for n
taup = @(V) (10^(-3))*(5+0.33./(exp((V+27)/10)+exp(-(V+102)/15)));  % Time constant [s] for p
tauq = @(V) (10^(-3))*(400./(exp((V+50)/15)+exp(-(V+50)/16))); % Time constant [s] for q
taur = @(kCA) (10^(-3))*2;  % Time-constant [s] for r
taua = @(V) (10^(-3))*(1+1./(1+exp((V+40)/0.5))); % Time-constant [s] for a
taub = @(V) (10^(-3))*(200./(exp((V+60)/30)+exp(-(V+40)/10))); % Time constant [s] for b
tauc = @(V) (10^(-3))*(45+10./(exp((V+27)/20)+exp(-(V+50)/15))); % Time constant [s] for c
taud1 = @(V) (10^(-3))*(400+500./(exp((V+40)/15)+exp(-(V+20)/20))); % Time constant [s] for d1
taud2 = @(kCA) (10^(-3))*130; % Time constant [s] for d2
end
if MODEL == 10
minf = @(V) 1./(1+exp(-(V+37)/7));            % rest m-value [-]
hinf = @(V) 1./(1+exp((V+41)/4));            % rest m-value [-]
tauah = @(V) 10^(3)*0.128*exp(-(V+46)/18);      
taubh = @(V) 10^(3)*4./(1+exp(-(V+23)/5));
tauh = @(V) 1./(tauah(V)+taubh(V));             % Time constant [s] for h
pinf = @(V) 1./(1+exp(-(V+60)/6.2));            % rest p-value [-]
rinf = @(V) 1./(1+exp((V+84)/4));            % rest r-value [-]
taur = @(V) (10^(-3))*0.15*(28+exp(-(V+25)/10.5));      % Time constant [s] for r
end
if MODEL == 11 || MODEL == 12
minf = @(V) 1./(1+exp(-(V+37)./10));            % rest m-value [-]
hinf = @(V) 1./(1+exp((V+58)./12));            % rest h-value [-]
ninf = @(V) 1./(1+exp(-(V+50)./14));           % rest n-value [-]
ainf = @(V) 1./(1+exp(-(V+57)./2));            % rest a-value [-]
rinf = @(V) 1./(1+exp((V+70)./2));            % rest r-value [-]
sinf = @(V) 1./(1+exp(-(V+35)./2));             % rest s-value [-]
tauh = @(V) (10^(-3)/0.05)*(0.05+0.27./(1+exp((V+40)/12))); % Time constant [s] for h
taun = @(V) (10^(-3)/0.1)*(0.05+0.27./(1+exp((V+40)/12)));  % Time constant [s] for n
taur = @(V) (10^(-3))*15;           % Time constant [s] for r
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
end
if MODEL == 14
am = @(V) 1000*Qam.^(0.1.*(Temp-Temp0)).*((2.5-0.1.*(V-Vm0))./(double((V-Vm0)~=25).*exp(2.5-0.1.*(V-Vm0))-1)+((V-Vm0)==25));
bm = @(V) 1000*Qbm.^(0.1.*(Temp-Temp0)).*4.*exp(-(V-Vm0)./18);
an = @(V) 1000*Qan.^(0.1.*(Temp-Temp0)).*(((0.1-0.01.*(V-Vm0))./(double((V-Vm0)~=10).*exp(1-0.1.*(V-Vm0))-1))+0.1*((V-Vm0)==10));
bn = @(V) 1000*Qbn.^(0.1.*(Temp-Temp0)).*0.125.*exp(-(V-Vm0)./80);
ah = @(V) 1000*Qah.^(0.1.*(Temp-Temp0)).*0.07.*exp(-(V-Vm0)./20);
bh = @(V) 1000*Qbh.^(0.1.*(Temp-Temp0)).*(1./(exp(3-0.1.*(V-Vm0))+1));
end

if SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6 % Create look-up (LU) tables to increase speed
Vresol = 1;                            % Voltage resolution (mV) 
kcCairesol = 1e-3;                      % 10^3*Cai concentration resolution (mol/m^3)
VLIMs = [-750; 500];            % Upper and lower table limit (mV)
kcCaiLIMs = [0; 1];             % Upper and lower table limit (10
% Note: The idea of look-up tables comes from Hines (1983), a voltage
% resolution of 1 mV results in 4 digit accuracy for the rates.
VrangeLU = (VLIMs(1):Vresol:VLIMs(2));          % Look-up range for V (mV)
kcCaiLU = (kcCaiLIMs(1):kcCairesol:kcCaiLIMs(2));            % Look-up range for 10^3*cCai (mol/m^3)
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 4 || MODEL == 5 || MODEL == 6 ||...
        MODEL == 7 || MODEL == 8 || MODEL == 13 || MODEL == 14
amLU = zeros(length(VrangeLU),1); bmLU = zeros(length(VrangeLU),1);
anLU = zeros(length(VrangeLU),1); bnLU = zeros(length(VrangeLU),1);
ahLU = zeros(length(VrangeLU),1); bhLU = zeros(length(VrangeLU),1);
end
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 6 || MODEL == 7 || MODEL == 8
pinfLU = zeros(length(VrangeLU),1); taupLU = zeros(length(VrangeLU),1);
end
if MODEL == 3 || MODEL == 4 || MODEL == 5 || MODEL == 8
sinfLU = zeros(length(VrangeLU),1); tausLU = zeros(length(VrangeLU),1);
uinfLU = zeros(length(VrangeLU),1); tauuLU = zeros(length(VrangeLU),1);
end
if MODEL == 4
winfLU = zeros(length(VrangeLU),1); tauwLU = zeros(length(VrangeLU),1);
end
if MODEL == 9
minfLU = zeros(length(VrangeLU),1); hinfLU = zeros(length(VrangeLU),1);
ninfLU = zeros(length(VrangeLU),1); ainfLU = zeros(length(VrangeLU),1);
binfLU = zeros(length(VrangeLU),1); rinfLU = zeros(length(kcCaiLU),1);
cinfLU = zeros(length(VrangeLU),1); pinfLU = zeros(length(VrangeLU),1);
qinfLU = zeros(length(VrangeLU),1); d1infLU = zeros(length(VrangeLU),1);
d2infLU = zeros(length(kcCaiLU),1);

taumLU = zeros(length(VrangeLU),1); tauhLU = zeros(length(VrangeLU),1); 
taunLU = zeros(length(VrangeLU),1); taupLU = zeros(length(VrangeLU),1);
tauqLU = zeros(length(VrangeLU),1); taurLU = zeros(length(kcCaiLU),1);
taucLU = zeros(length(VrangeLU),1); tauaLU = zeros(length(VrangeLU),1);
taubLU = zeros(length(VrangeLU),1); taud1LU = zeros(length(VrangeLU),1);
taud2LU = zeros(length(kcCaiLU),1);
end
if MODEL == 10
minfLU = zeros(length(VrangeLU),1); hinfLU=zeros(length(VrangeLU),1);
tauhLU = zeros(length(VrangeLU),1); pinfLU = zeros(length(VrangeLU),1);
rinfLU = zeros(length(VrangeLU),1); taurLU = zeros(length(VrangeLU),1);
end
if MODEL == 11 || MODEL == 12
minfLU = zeros(length(VrangeLU),1); hinfLU = zeros(length(VrangeLU),1); 
ninfLU = zeros(length(VrangeLU),1); ainfLU = zeros(length(VrangeLU),1); 
rinfLU = zeros(length(VrangeLU),1); sinfLU = zeros(length(VrangeLU),1); 
tauhLU = zeros(length(VrangeLU),1); taunLU = zeros(length(VrangeLU),1); 
taurLU = zeros(length(VrangeLU),1);  
end
if MODEL == 13
apLU = zeros(length(VrangeLU),1); bpLU = zeros(length(VrangeLU),1);
end
for iLU = 1:length(VrangeLU)
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 4 || MODEL == 5 || MODEL == 6 ||...
        MODEL == 7 || MODEL == 8 || MODEL == 13 || MODEL == 14
amLU(iLU) = am(VrangeLU(iLU)); bmLU(iLU) = bm(VrangeLU(iLU)); 
anLU(iLU) = an(VrangeLU(iLU)); bnLU(iLU) = bn(VrangeLU(iLU));
ahLU(iLU) = ah(VrangeLU(iLU)); bhLU(iLU) = bh(VrangeLU(iLU));
end
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 6 || MODEL == 7 || MODEL == 8
pinfLU(iLU) = pinf(VrangeLU(iLU)); taupLU(iLU) = taup(VrangeLU(iLU));
end
if MODEL == 3 || MODEL == 4  || MODEL == 5 || MODEL == 8
sinfLU(iLU) = sinf(VrangeLU(iLU)); tausLU(iLU) = taus(VrangeLU(iLU));
uinfLU(iLU) = uinf(VrangeLU(iLU)); tauuLU(iLU) = tauu(VrangeLU(iLU));
end
if MODEL == 4
winfLU(iLU) = winf(VrangeLU(iLU)); tauwLU(iLU) = tauw(VrangeLU(iLU));
end
if MODEL == 9
minfLU(iLU) = minf(VrangeLU(iLU)); hinfLU(iLU) = hinf(VrangeLU(iLU));
ninfLU(iLU) = ninf(VrangeLU(iLU)); ainfLU(iLU) = ainf(VrangeLU(iLU));
binfLU(iLU) = binf(VrangeLU(iLU)); 
cinfLU(iLU) = cinf(VrangeLU(iLU)); pinfLU(iLU) = pinf(VrangeLU(iLU));
qinfLU(iLU) = qinf(VrangeLU(iLU)); d1infLU(iLU) = d1inf(VrangeLU(iLU));

taumLU(iLU) = taum(VrangeLU(iLU)); tauhLU(iLU) = tauh(VrangeLU(iLU)); 
taunLU(iLU) = taun(VrangeLU(iLU)); taupLU(iLU) = taup(VrangeLU(iLU));
tauqLU(iLU) = tauq(VrangeLU(iLU)); 
taucLU(iLU) = tauc(VrangeLU(iLU)); tauaLU(iLU) = taua(VrangeLU(iLU));
taubLU(iLU) = taub(VrangeLU(iLU)); taud1LU(iLU) = taud1(VrangeLU(iLU));
end
if MODEL == 10
minfLU(iLU) = minf(VrangeLU(iLU)); hinfLU(iLU) = hinf(VrangeLU(iLU));
tauhLU(iLU) = tauh(VrangeLU(iLU)); pinfLU(iLU) = pinf(VrangeLU(iLU));
rinfLU(iLU) = rinf(VrangeLU(iLU)); taurLU(iLU) = taur(VrangeLU(iLU));
end
if MODEL == 11 || MODEL == 12
minfLU(iLU) = minf(VrangeLU(iLU)); hinfLU(iLU) = hinf(VrangeLU(iLU)); 
ninfLU(iLU) = ninf(VrangeLU(iLU)); ainfLU(iLU) = ainf(VrangeLU(iLU)); 
rinfLU(iLU) = rinf(VrangeLU(iLU)); sinfLU(iLU) = sinf(VrangeLU(iLU)); 
tauhLU(iLU) = tauh(VrangeLU(iLU)); taunLU(iLU) = taun(VrangeLU(iLU)); 
taurLU(iLU) = taur(VrangeLU(iLU));  
end
if MODEL == 13
apLU(iLU) = ap(VrangeLU(iLU)); bpLU(iLU) = bp(VrangeLU(iLU));
end
end
for iLU = 1:length(kcCaiLU)
if MODEL == 9
d2infLU(iLU) =  d2inf(kcCaiLU(iLU)); rinfLU(iLU) = rinf(kcCaiLU(iLU));
taud2LU(iLU) = taud2(kcCaiLU(iLU)); taurLU(iLU) = taur(kcCaiLU(iLU));
end
end
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 4 || MODEL == 5 || MODEL == 6 ||...
        MODEL == 7 || MODEL == 8 || MODEL == 13 || MODEL == 14
ampbmLU = amLU+bmLU; anpbnLU = anLU+bnLU; ahpbhLU = ahLU+bhLU;
end
if MODEL == 13
appbpLU = apLU+bpLU;
end
if SpeedUp == 2
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 4 || MODEL == 5 || MODEL == 6 ||...
        MODEL == 7 || MODEL == 8 || MODEL == 13 || MODEL == 14
am = @(V) interp1(VrangeLU',amLU,V); ampbm = @(V) interp1(VrangeLU',ampbmLU,V);
an = @(V) interp1(VrangeLU',anLU,V); anpbn = @(V) interp1(VrangeLU',anpbnLU,V);
ah = @(V) interp1(VrangeLU',ahLU,V); ahpbh = @(V) interp1(VrangeLU',ahpbhLU,V);
end
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 6 || MODEL == 7 || MODEL == 8
pinf = @(V) interp1(VrangeLU',pinfLU,V); taup = @(V) interp1(VrangeLU',taupLU,V);
end
if MODEL == 3 || MODEL == 4 || MODEL == 5 || MODEL == 8
sinf = @(V) interp1(VrangeLU',sinfLU,V); taus = @(V) interp1(VrangeLU',tausLU,V);
uinf = @(V) interp1(VrangeLU',uinfLU,V); tauu = @(V) interp1(VrangeLU',tauuLU,V);
end
if MODEL == 4
winf = @(V) interp1(VrangeLU',winfLU,V); tauw = @(V) interp1(VrangeLU',tauwLU,V);
end
if MODEL == 9
minf = @(V) interp1(VrangeLU',minfLU,V); hinf = @(V) interp1(VrangeLU',hinfLU,V);
ninf = @(V) interp1(VrangeLU',ninfLU,V); ainf = @(V) interp1(VrangeLU',ainfLU,V);
binf = @(V) interp1(VrangeLU',binfLU,V); rinf = @(kcCai) interp1(kcCaiLU',rinfLU,kcCai);
cinf = @(V) interp1(VrangeLU',cinfLU,V); pinf = @(V) interp1(VrangeLU',pinfLU,V);
qinf = @(V) interp1(VrangeLU',qinfLU,V); d1inf = @(V) interp1(VrangeLU',d1infLU,V);
d2inf = @(kcCai) interp1(kcCaiLU',d2infLU,kcCai);

taum = @(V) interp1(VrangeLU',taumLU,V); tauh = @(V) interp1(VrangeLU',tauhLU,V); 
taun = @(V) interp1(VrangeLU',taunLU,V); taup = @(V) interp1(VrangeLU',taupLU,V);
tauq = @(V) interp1(VrangeLU',tauqLU,V); taur = @(kcCai) interp1(kcCaiLU',taurLU,kcCai);
tauc = @(V) interp1(VrangeLU',taucLU,V); taua = @(V) interp1(VrangeLU',tauaLU,V);
taub = @(V) interp1(VrangeLU',taubLU,V); taud1 = @(V) interp1(VrangeLU',taud1LU,V);
taud2 = @(kcCai) interp1(kcCaiLU',taud2LU,kcCai);
end
if MODEL == 10
minf = @(V) interp1(VrangeLU',minfLU,V);  hinf = @(V) interp1(VrangeLU',hinfLU,V);
tauh = @(V) interp1(VrangeLU',tauhLU,V); pinf = @(V) interp1(VrangeLU',pinfLU,V);
rinf = @(V) interp1(VrangeLU',rinfLU,V); taur = @(V) interp1(VrangeLU',taurLU,V);
end
if MODEL == 11 || MODEL == 12
minf = @(V) interp1(VrangeLU',minfLU,V); hinf = @(V) interp1(VrangeLU',hinfLU,V); 
ninf = @(V) interp1(VrangeLU',ninfLU,V); ainf = @(V) interp1(VrangeLU',ainfLU,V); 
rinf = @(V) interp1(VrangeLU',rinfLU,V); sinf = @(V) interp1(VrangeLU',sinfLU,V); 
tauh = @(V) interp1(VrangeLU',tauhLU,V); taun = @(V) interp1(VrangeLU',taunLU,V); 
taur = @(V) interp1(VrangeLU',taurLU,V);  
end
if MODEL == 13
ap = @(V) interp1(VrangeLU',apLU,V); appbp = @(V) interp1(VrangeLU',appbpLU,V);
end
elseif SpeedUp == 3 || SpeedUp == 6
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 4 || MODEL == 5 || MODEL == 6 ||...
        MODEL == 7 || MODEL == 8 || MODEL == 13 || MODEL == 14
am = @(V) interp1qr(VrangeLU',amLU,V); ampbm = @(V) interp1qr(VrangeLU',ampbmLU,V);
an = @(V) interp1qr(VrangeLU',anLU,V); anpbn = @(V) interp1qr(VrangeLU',anpbnLU,V);
ah = @(V) interp1qr(VrangeLU',ahLU,V); ahpbh = @(V) interp1qr(VrangeLU',ahpbhLU,V);
end
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 6 || MODEL == 7 || MODEL == 8
pinf = @(V) interp1qr(VrangeLU',pinfLU,V); taup = @(V) interp1qr(VrangeLU',taupLU,V);
end
if MODEL == 3 || MODEL == 4 || MODEL == 5 || MODEL == 8
sinf = @(V) interp1qr(VrangeLU',sinfLU,V); taus = @(V) interp1qr(VrangeLU',tausLU,V);
uinf = @(V) interp1qr(VrangeLU',uinfLU,V); tauu = @(V) interp1qr(VrangeLU',tauuLU,V);
end
if MODEL == 4
winf = @(V) interp1qr(VrangeLU',winfLU,V); tauw = @(V) interp1qr(VrangeLU',tauwLU,V);
end
if MODEL == 9
minf = @(V) interp1qr(VrangeLU',minfLU,V); hinf = @(V) interp1qr(VrangeLU',hinfLU,V);
ninf = @(V) interp1qr(VrangeLU',ninfLU,V); ainf = @(V) interp1qr(VrangeLU',ainfLU,V);
binf = @(V) interp1qr(VrangeLU',binfLU,V);rinf = @(kcCai) interp1qr(kcCaiLU',rinfLU,kcCai);
cinf = @(V) interp1qr(VrangeLU',cinfLU,V);pinf = @(V) interp1qr(VrangeLU',pinfLU,V);
qinf = @(V) interp1qr(VrangeLU',qinfLU,V);d1inf = @(V) interp1qr(VrangeLU',d1infLU,V);
d2inf = @(kcCai) interp1qr(kcCaiLU',d2infLU,kcCai);

taum = @(V) interp1qr(VrangeLU',taumLU,V); tauh = @(V) interp1qr(VrangeLU',tauhLU,V); 
taun = @(V) interp1qr(VrangeLU',taunLU,V); taup = @(V) interp1qr(VrangeLU',taupLU,V);
tauq = @(V) interp1qr(VrangeLU',tauqLU,V); taur = @(kcCai) interp1qr(kcCaiLU',taurLU,kcCai);
tauc = @(V) interp1qr(VrangeLU',taucLU,V); taua = @(V) interp1qr(VrangeLU',tauaLU,V);
taub = @(V) interp1qr(VrangeLU',taubLU,V); taud1 = @(V) interp1qr(VrangeLU',taud1LU,V);
taud2 = @(kcCai) interp1qr(kcCaiLU',taud2LU,kcCai);
end
if MODEL == 10
minf = @(V) interp1qr(VrangeLU',minfLU,V);  hinf = @(V) interp1qr(VrangeLU',hinfLU,V);
tauh = @(V) interp1qr(VrangeLU',tauhLU,V); pinf = @(V) interp1qr(VrangeLU',pinfLU,V);
rinf = @(V) interp1qr(VrangeLU',rinfLU,V); taur = @(V) interp1qr(VrangeLU',taurLU,V);
end
if MODEL == 11 || MODEL == 12
minf = @(V) interp1qr(VrangeLU',minfLU,V); hinf = @(V) interp1qr(VrangeLU',hinfLU,V); 
ninf = @(V) interp1qr(VrangeLU',ninfLU,V); ainf = @(V) interp1qr(VrangeLU',ainfLU,V); 
rinf = @(V) interp1qr(VrangeLU',rinfLU,V); sinf = @(V) interp1qr(VrangeLU',sinfLU,V); 
tauh = @(V) interp1qr(VrangeLU',tauhLU,V); taun = @(V) interp1qr(VrangeLU',taunLU,V); 
taur = @(V) interp1qr(VrangeLU',taurLU,V);  
end
if MODEL == 13
ap = @(V) interp1qr(VrangeLU',apLU,V); appbp = @(V) interp1qr(VrangeLU',appbpLU,V);
end
elseif SpeedUp == 4 || SpeedUp == 5
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 4 || MODEL == 5 || MODEL == 6 ||...
        MODEL == 7 || MODEL == 8 || MODEL == 13 || MODEL == 14
am = @(V) nakeinterp1(VrangeLU',amLU,V); ampbm = @(V) nakeinterp1(VrangeLU',ampbmLU,V);
an = @(V) nakeinterp1(VrangeLU',anLU,V); anpbn = @(V) nakeinterp1(VrangeLU',anpbnLU,V);
ah = @(V) nakeinterp1(VrangeLU',ahLU,V); ahpbh = @(V) nakeinterp1(VrangeLU',ahpbhLU,V);
end
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 6 || MODEL == 7 || MODEL == 8
pinf = @(V) nakeinterp1(VrangeLU',pinfLU,V); taup = @(V) nakeinterp1(VrangeLU',taupLU,V);
end
if MODEL == 3 || MODEL == 4  || MODEL == 5 || MODEL == 8
sinf = @(V) nakeinterp1(VrangeLU',sinfLU,V); taus = @(V) nakeinterp1(VrangeLU',tausLU,V);
uinf = @(V) nakeinterp1(VrangeLU',uinfLU,V); tauu = @(V) nakeinterp1(VrangeLU',tauuLU,V);
end
if MODEL == 4
winf = @(V) nakeinterp1(VrangeLU',winfLU,V); tauw = @(V) nakeinterp1(VrangeLU',tauwLU,V);
end
if MODEL == 9
minf = @(V) nakeinterp1(VrangeLU',minfLU,V); hinf = @(V) nakeinterp1(VrangeLU',hinfLU,V);
ninf = @(V) nakeinterp1(VrangeLU',ninfLU,V); ainf = @(V) nakeinterp1(VrangeLU',ainfLU,V);
binf = @(V) nakeinterp1(VrangeLU',binfLU,V);rinf = @(kcCai) nakeinterp1(kcCaiLU',rinfLU,kcCai);
cinf = @(V) nakeinterp1(VrangeLU',cinfLU,V);pinf = @(V) nakeinterp1(VrangeLU',pinfLU,V);
qinf = @(V) nakeinterp1(VrangeLU',qinfLU,V);d1inf = @(V) nakeinterp1(VrangeLU',d1infLU,V);
d2inf = @(kcCai) nakeinterp1(kcCaiLU',d2infLU,kcCai);

taum = @(V) nakeinterp1(VrangeLU',taumLU,V); tauh = @(V) nakeinterp1(VrangeLU',tauhLU,V); 
taun = @(V) nakeinterp1(VrangeLU',taunLU,V); taup = @(V) nakeinterp1(VrangeLU',taupLU,V);
tauq = @(V) nakeinterp1(VrangeLU',tauqLU,V); taur = @(kcCai) nakeinterp1(kcCaiLU',taurLU,kcCai);
tauc = @(V) nakeinterp1(VrangeLU',taucLU,V); taua = @(V) nakeinterp1(VrangeLU',tauaLU,V);
taub = @(V) nakeinterp1(VrangeLU',taubLU,V); taud1 = @(V) nakeinterp1(VrangeLU',taud1LU,V);
taud2 = @(kcCai) nakeinterp1(kcCaiLU',taud2LU,kcCai);
end
if MODEL == 10
minf = @(V) nakeinterp1(VrangeLU',minfLU,V);  hinf = @(V) nakeinterp1(VrangeLU',hinfLU,V);
tauh = @(V) nakeinterp1(VrangeLU',tauhLU,V); pinf = @(V) nakeinterp1(VrangeLU',pinfLU,V);
rinf = @(V) nakeinterp1(VrangeLU',rinfLU,V); taur = @(V) nakeinterp1(VrangeLU',taurLU,V);
end
if MODEL == 11 || MODEL == 12
minf = @(V) nakeinterp1(VrangeLU',minfLU,V); hinf = @(V) nakeinterp1(VrangeLU',hinfLU,V); 
ninf = @(V) nakeinterp1(VrangeLU',ninfLU,V); ainf = @(V) nakeinterp1(VrangeLU',ainfLU,V); 
rinf = @(V) nakeinterp1(VrangeLU',rinfLU,V); sinf = @(V) nakeinterp1(VrangeLU',sinfLU,V); 
tauh = @(V) nakeinterp1(VrangeLU',tauhLU,V); taun = @(V) nakeinterp1(VrangeLU',taunLU,V); 
taur = @(V) nakeinterp1(VrangeLU',taurLU,V);  
end
if MODEL == 13
ap = @(V) nakeinterp1(VrangeLU',apLU,V); appbp = @(V) nakeinterp1(VrangeLU',appbpLU,V);
end
end
end

% 2.2. BLS functions
Cm = @(Z) FunCm(Z,Cm0,a,delta); % Membrane capacitance (F/m^2)
if xfs ~= 1
Cm = @(Z) xfs*Cm(Z)+(1-xfs)*Cm0; % Membrane capacitance with partial coverage (F/m^2)   
end
Pec = @(V,Z) (-a^2/(Z.^2+a^2))*((Cm(Z)*0.001*V)^2/(2*eps0*epsr)); % Electrostatic pressure (Pa)
PecQ = @(Q,Z) (-a^2/(Z.^2+a^2))*(Q^2/(2*eps0*epsr));
Va = @(Z) pi*a^2*delta*(1+(Z./(3*delta)).*(Z.^2/a^2+3)); % Intraleaflet-volume (m^3)
PS = @(Z) (2*ks*Z.^3)./(a^2*(a^2+Z.^2)); % Tension-pressure [Pa]
R = @(Z) (a^2+Z.^2)/(2.*Z); % Radius of curvature [m]
z = @(r,Z) Funz(r,Z,R(Z)); % Curvature [m]
f = @(r,Z) Ar*((deltax./(2*z(r,Z)+delta)).^x-(deltax./(2*z(r,Z)+delta)).^y); % [Pa]
Pm = @(Z) (2./(Z.^2+a^2)).*integral(@(r)(r.*f(r,Z)),0,a); % Attraction/repulsion moulecular pressure [Pa]
% Pm(Z) is not properly defined for Z < -delta/2. Z should be > -delta/2!
if UpSpeed == 1 || UpSpeed == 2 || UpSpeed == 3
if DISPLAY, disp('Initialising molecular force vector...'); end
ll = -delta/2; infinites = 1*10^(-16); ss = 10^(-12); rl = 50*10^(-9);
PmRange = (ll+infinites:ss:rl); 
PmI = zeros(length(PmRange),1);
for iR = 1:length(PmRange)
PmI(iR) = Pm(PmRange(iR));
end
if SpeedUp == 6 || SpeedUp == 3
Pm = @(Z) interp1qr(PmRange',PmI,Z);
else
Pm = @(Z) nakeinterp1(PmRange',PmI,Z);
end
if DISPLAY 
disp('...Initialisation complete');
end
end
if UpSpeed == 2 || UpSpeed == 3
if DISPLAY, disp('Initialising intraleaflet volume...'); end
ll2 = -delta/2; infinites2 = 1*10^(-16); ss2 = 10^(-12); rl2 = 50*10^(-9);
VaRange = (ll2+infinites2:ss2:rl2); 
VaI = zeros(length(VaRange),1);
for iVa = 1:length(VaRange)
VaI(iVa) = Va(VaRange(iVa));
end
if SpeedUp == 6 || SpeedUp == 3
Va = @(Z) interp1qr(VaRange',VaI,Z);
else
Va = @(Z) nakeinterp1(VaRange',VaI,Z);
end
if DISPLAY 
disp('...Initialisation complete');
end    
end
if UpSpeed == 3
if DISPLAY, disp('Initialising capacitance...'); end
ll = -delta/2; infinites = 1*10^(-16); ss = 10^(-12); rl = 50*10^(-9);
CmRange = (ll+infinites:ss:rl); 
CmI = zeros(length(CmRange),1);
for iR = 1:length(CmRange)
CmI(iR) = Cm(CmRange(iR));
end
if SpeedUp == 6 || SpeedUp == 3
Cm = @(Z) interp1qr(CmRange',CmI,Z);
else
Cm = @(Z) nakeinterp1(CmRange',CmI,Z);
end
if DISPLAY 
disp('...Initialisation complete');
end
end    
if ERROR
Pm = @(Z) FunPm(Z,delta,Pm);
end
dCmdZ = @(Z) FundCmdZ(Z,Cm0,a,delta); % derivative dCm/dZ [F/m]
PinLIN = @(Z) Po*(1+(Z./(3*delta)).*(Z.^2/a^2+3)).^(-kappa); %Linear approximation of Pin 
PmdPin = @(Z) Pm(Z)*(1+(Z./(3*delta))*(Z.^2/a^2+3))^(kappa)/Po; % Pm divided by Pinlin
Pin = @(na,Z) (na*Rg*Temp)./(Va(Z));% Internal pressure [Pa]
S = @(Z) pi*(a^2+Z.^2);% Surface of a single leaflet [m^2]
fVCa = @(cCai) 10^(3)*((Rg*Temp)/(2*Far))*log(cCae./cCai); % Nernst equation for Ca-potential [mV] (if not assumed constant) 

% 3. Initial conditions and timespan
% Y = [V,m,n,p,h,s,u,Z,dZ,Pin,na,Ca]
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 4 || MODEL == 5 || MODEL == 6 ||...
        MODEL == 7 || MODEL == 8 || MODEL == 14
m0 = am(Vm0)/(am(Vm0)+bm(Vm0));
n0 = an(Vm0)/(an(Vm0)+bn(Vm0));
h0 = ah(Vm0)/(ah(Vm0)+bh(Vm0));
end
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 6 || MODEL == 7 || MODEL == 8
p0 = pinf(Vm0);
end
if MODEL == 3 || MODEL == 4 || MODEL == 5 || MODEL == 8
s0 = sinf(Vm0);
u0 = uinf(Vm0);
end
if MODEL == 9
%     Vm0fun = @(V) Gl*(V-Vl)+Gna*minf(V).^3.*hinf(V).*(V-Vna)+Gk*ninf(V).^4.*(V-Vk)+...
%      GT*pinf(V).^2.*qinf(V).*(V-fVCa(5e-6))+GCa.*rinf(V).^2.*(V-Vk)+GA*ainf(V).^2.*binf(V).*(V-Vk)+...
%      +GL*cinf(V).^2.*d1inf(V).*d2inf(5e-6).*(V-fVCa(5e-6));
%  plot(VrangeLU',Vm0fun(VrangeLU'));
%  Vm0 = fzero(Vm0fun,[-200,100]); -> Vm0 = -42 mV

% % ---------VALIDATION --- PLOT OF STEAD STATE CURRENT ----------------------
% Vrange = (-90:1:-20); Icontrol = zeros(length(Vrange),1); Ittx = zeros(length(Vrange),1);
% Iconmttx = zeros(length(Vrange),1);
% funcCai = @(V,cCai) -0.001*tauCa*(GT*pinf(V).^2.*qinf(V).*(V-fVCa(cCai))+...
%     GL*cinf(V).^2.*d1inf(V).*d2inf(cCai).*(V-fVCa(cCai)))/(2*Far*10236*10^(-9))-cCai;
% 
% for i = 1:length(Vrange)
%     Pot = Vrange(i);
% cCaiSol = fzero(@(cCai) funcCai(Pot,cCai),[1e-12 1e-3]);
% Icontrol(i) = 0.001*(Gl*(Pot-Vl)+Gna*minf(Pot).^3.*hinf(Pot).*(Pot-Vna)+Gk*ninf(Pot).^4.*(Pot-Vk)+...
%   GT*pinf(Pot).^2.*qinf(Pot).*(Pot-fVCa(cCaiSol))+GCa.*rinf(cCaiSol).^2.*(Pot-Vk)+GA*ainf(Pot).^2.*binf(Pot).*(Pot-Vk)+...
%  +GL*cinf(Pot).^2.*d1inf(Pot).*d2inf(cCaiSol).*(Pot-fVCa(cCaiSol)));
% Ittx(i) = Icontrol(i)-0.001*Gna*minf(Pot).^3.*hinf(Pot).*(Pot-Vna);
% Iconmttx(i) = 0.001*Gna*minf(Pot).^3.*hinf(Pot).*(Pot-Vna);
% end
% figure;
% hold on;
% plot(Vrange',Icontrol);
% plot(Vrange',Ittx);
% plot(Vrange',Iconmttx);
% xlabel('Potential [mV]');
% ylabel('Steady-state current [A/m^2]');
% xlim([-90,-20]);
% %ylim([-210,300]);
% hold off;
% 
m0 = minf(Vm0); h0 = hinf(Vm0); n0 = ninf(Vm0); a0 = ainf(Vm0);
b0 = binf(Vm0); c0 = cinf(Vm0); p0 = pinf(Vm0); 
q0 = qinf(Vm0); d10 = d1inf(Vm0);
end
if MODEL == 10
%Vm0fun = @(V) Gl*(V-Vl)+Gna*minf(V).^3.*hinf(V).*(V-Vna)+...
%    Gk*(0.75*(1-hinf(V))).^4.*(V-Vk)+GT*pinf(V).^2.*rinf(V).*(V-VT);
%  plot(VrangeLU',Vm0fun(VrangeLU'));
%  Vm0 = fzero(Vm0fun,[-200,100]);  -> Vm0 = -65

m0 = minf(Vm0); h0 = hinf(Vm0); p0 = pinf(Vm0); r0 = rinf(Vm0);
end
if MODEL == 11 || MODEL == 12
%Ca = @(V) (-1/15)*(0.1*GT*ainf(V).^3.*rinf(V).*(V-VT)+0.1*GCa.*sinf(V).^2.*(V-VCa));
%Vm0fun = @(V) Gl*(V-Vl)+Gna*minf(V).^3.*hinf(V).*(V-Vna)+...
%    Gk*ninf(V).^4.*(V-Vk)+GT*ainf(V).^3.*rinf(V).*(V-VT)+...
%    GCa.*sinf(V).^2.*(V-VCa)+Gahp*(V-Vk).*(Ca(V)./(Ca(V)+10));
%  plot(VrangeLU',Vm0fun(VrangeLU'));
%  Vm0 = fzero(Vm0fun,[-200,100]);    -> vm0 =-66 
%  CA0 = Ca(Vm0);   -> CA0 = 4.32e-7

m0 = minf(Vm0); h0 = hinf(Vm0); n0 = ninf(Vm0); a0 = ainf(Vm0); r0 = rinf(Vm0);
s0 = sinf(Vm0); 
end
if MODEL == 13
% minf = @(V) am(V)./(am(V)+bm(V)); hinf = @(V) ah(V)./(ah(V)+bh(V)); 
% ninf = @(V) an(V)./(an(V)+bn(V)); pinf = @(V) ap(V)./(ap(V)+bp(V)); 
% Vm0fun = @(V) Gl*(V-Vl)+Gna*minf(V).^3.*hinf(V).*(V-Vna)+...
%    Gk*ninf(V).^4.*(V-Vk)+Gk*pinf(V).*(V-Vm);
%  plot(VrangeLU',Vm0fun(VrangeLU'));
%  Vm0 = fzero(Vm0fun,[-200,100]);  -> Vm0 = -87
m0 = am(Vm0)/(am(Vm0)+bm(Vm0)); n0 = an(Vm0)/(an(Vm0)+bn(Vm0));
h0 = ah(Vm0)/(ah(Vm0)+bh(Vm0)); p0 = ap(Vm0)/(ap(Vm0)+bp(Vm0));
end
if MODEL == 4 || MODEL == 5
cCai0 = 240*10^(-6);             % Rest concentration Ca (mol/m^3)
elseif MODEL == 9
cCai0 = 5*10^(-6);
d20 = d2inf(cCai0);
r0 = rinf(cCai0);
elseif MODEL == 11 || MODEL == 12
CA0 = 4.32e-07;           % Rest concentration Ca (uM = umol/l)
end
if MODEL == 4
% dw/dt=alpha(V)*(1-w-wLock)-beta(V)*w = (winf*(1-wLock)-w)/tauw
% w0 = winf*(1-wLock0) = winf*(1-((k3*P0)/k4)w0)
% w0*(1+winf*(k3*P0)/k4) = winf -> w0 = winf/(1+winf*(k3*P0)/k4)
% dP/dt = k1*(1-P)*cCa^2-k2*P
% dwLock/dt = k3*w*P-k4*wLock
hProtein0 = (k1*cCai0^2)/(k2+k1*cCai0^2); % Intitial condition for w-gate protein (-)
w0 = winf(Vm0)/(1+winf(Vm0)*((k3*hProtein0)/k4));
wLock0 = (k3/k4)*(w0*hProtein0);   % Initial condition for locked w-gate
end
Z0 = 0;
dZ0 = 0;
Pin0 = Po;
na0 = (Pin0*Va(Z0))/(Rg*Temp);  % Mole air in BLS
Ca0 = Ca.*ones(1,Nout);         % External molar air concentration (mole/m^3)
Tupdate = 500*10^(-6);          % Update period:  every Tupdate (s) the charge Cm*Vm is updated (Plaksin et al., 2014, appendix)
if Tupdate < 10/USfreq
    Wnum = floor(Tupdate*USfreq);
    disp(' ');
    disp(['Warning! Only ' num2str(Wnum) ' periods per update cycle.'])
    disp(' ');
end
if USdc ~= 1
maxErrorDC = 100*((USdc*USprp+Tupdate)/(USprp)-USdc)/USdc;      % Max error on duty cycle [%]
if maxErrorDC >= 2.5
    disp(' '); 
    disp(['Warning! Maximum discretisation error on duty cycle is ' num2str(maxErrorDC) '%']);
    disp(' ');
end
end
dtUS = 0.025/USfreq;            % Discretisation time (s)
%dtUS = 0.1/USfreq;
%dtUS = 0.005/USfreq;
dtES = 10*10^(-6);              % ES discr. time (s)
dtlin = 0.05*dtUS;

Q0 = Cm(Z0)*(10^(-3)*Vm0);      % Update charge [C]

if MODEL == 1 || MODEL == 2 || MODEL == 6 || MODEL == 7 || MODEL == 13 
Y0 = [Q0 m0 n0 p0 h0]';
elseif MODEL == 3 || MODEL == 8
Y0 = [Q0 m0 n0 p0 h0 s0 u0]';
elseif MODEL == 4
Y0 = [Q0 m0 n0 h0 s0 u0 w0 wLock0 hProtein0 cCai0]';
elseif MODEL == 5
Y0 = [Q0 m0 n0 h0 s0 u0 cCai0]';    
elseif MODEL == 9
    if Charges
Y0 = [Q0 m0 n0 p0 h0 q0 r0 a0 b0 c0 d10 d20 cCai0 zeros(1,7)]';
    else
Y0 = [Q0 m0 n0 p0 h0 q0 r0 a0 b0 c0 d10 d20 cCai0]';
    end
elseif MODEL == 10
Y0 = [Q0 h0 r0]';
elseif MODEL == 11 || MODEL == 12
Y0 = [Q0 n0 h0 r0 CA0]';
elseif MODEL == 14
Y0 = [Q0 m0 n0 h0]';
end
X0 = [Z0 dZ0 Pin0 na0 Ca0]';
%Tvalueslin = (max(0,USpstart):dtlin:ceil(min(Tsim,USpstart+USpd)/dtlin)*dtlin);
Tvalueslin = @(it) max(0,USpstart)+(it-1)*dtlin;

% 4. Solver
global reverseStr; 
global temp; 
if DISPLAY, disp('Solver started'); end
SearchRange = [USibegin USiend]; PrecisionCheck = 0;
if SearchMode % SearchMode = 1
IIpa = (SearchRange(1)+SearchRange(2))/2; % W/m^2
else % SearchMode = 0
if ~any(SearchRange)
IIpa = 1000;            % Don't choose too small, because it will slow down calculation -> Linear BLS is relatively slow
else
IIpa = SearchRange(~~SearchRange)*10^(3-2*find(SearchRange));
end
end
while ~SearchMode || ~PrecisionCheck 
if MODE == 1
USPaT = @(t) sqrt(2*rhol*c*IIpa)*USstep(t);
elseif MODE == 2
USPaT = @ (t) USPa*USstep(t); 
end
BLSlin =  @(t,Z,V) PmdPin(Z)+1-(1+(Z./(3*delta))*(Z.^2/a^2+3))^(kappa)+...
    (USPaT(t)/Po)*sin(omega*t)*(1+(Z./(3*delta))*(Z.^2/a^2+3))^(kappa)+...
    (Pec(V,Z)/Po)*(1+(Z./(3*delta))*(Z.^2/a^2+3))^(kappa);  % Differential-algebraic equation
BLSlinQ =  @(t,Z,Q) PmdPin(Z)+1-(1+(Z./(3*delta))*(Z.^2/a^2+3))^(kappa)+...
    (USPaT(t)/Po)*sin(omega*t)*(1+(Z./(3*delta))*(Z.^2/a^2+3))^(kappa)+...
    (PecQ(Q,Z)/Po)*(1+(Z./(3*delta))*(Z.^2/a^2+3))^(kappa);
iV0 = Vm0; iZ0 = Z0; X = X0'; TvaluesY = [0]; %#ok<*NBRAK>
Q0 = Cm(iZ0)*(10^(-3)*iV0);      % Update charge [C]
U0 = Y0; Y=Y0';
if MemorySaveMode == 3
    Y = Y0(1);
end
tCURRENT = 0; % Current time (s)
DispPart = 0;
DispNrPart=double(USpstart>=Tsim||(USpstart+USpd)<=0||(USpstart==0&&USpd==Tsim))+...
   2*double(USpstart<=0&&((USpstart+USpd)>0)&&((USpstart+USpd)<Tsim))+...
   2*double(((USpstart+USpd)>=Tsim)&&(USpstart>0&&USpstart<Tsim))+...
   3*double(USpstart>0&&((USpstart+USpd)<Tsim));

while tCURRENT < Tsim
    DispPart = DispPart+1;
    if DISPLAY
    disp(' ');
    disp(['Part ' num2str(DispPart) ' of ' num2str(DispNrPart)]);
    disp(' ');
    end
    if tCURRENT < USpstart 
        ZONE = 1;
        tNICEc = [tCURRENT, min(USpstart,Tsim)]; % Complete (c) tNICE interval
    elseif tCURRENT >= USpstart && tCURRENT < USpstart+USpd
        ZONE = 2;
        tNICEc = [tCURRENT, min(USpstart+USpd,Tsim)];
    elseif tCURRENT >= USpstart+USpd
        ZONE = 1;
        tNICEc = [tCURRENT,Tsim];
    end
if SpeedUp == 1 || SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
    temp = [];
end
if ZONE == 2
nUP = ceil((tNICEc(2)-tNICEc(1))/Tupdate); it=2; epsT = 2/USfreq;
if (nUP-1)*Tupdate+tNICEc(1)>=(tNICEc(2)-epsT), nUP = nUP-1; end % Line added to solve numerical division bug 
itB=it; iCa0 = Ca0; itT = 0; 
Event = abs(iZ0)<abs(Hmin); Checkpoint = 0; Checkpoint2 = 0; 
acorrpass = 0; acorrpass2 = 0; Pass=0; PassPoint = 1;
TvaluesX = Tvalueslin(1); 
MemX = [Z0]; 
reverseStr3 = '';
else
nUP = 1; 
reverstreStr3 = ''; 
end

for iUP=1:nUP
if ZONE == 2
    dt = dtUS;
    if MemorySaveMode 
    TvaluesX = TvaluesX(end);
    X = X(end,:); it=2; itB=it;
    end
    if DISPLAY == 2
UpdateStr = ['Solver updated: update nr. ' num2str(iUP)];
disp(' ');
disp(UpdateStr);
    elseif DISPLAY == 1
Progress3 = 100*iUP/nUP;  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress3); 
fprintf([reverseStr3, msg]);
reverseStr3 = repmat(sprintf('\b'), 1, length(msg));  
    end
    % A. Calculate Z,dZ,na,Pin
    % A.1. Linear initialization
tspan = tNICEc(1)+[(iUP-1)*Tupdate iUP*Tupdate];
if iUP==nUP
    tspan(2) = tNICEc(2);
end
infinit = 1*10^(-12);           % Small number for fzero range
if Event
Zold=iZ0;
while Tvalueslin(it+itT)<=tspan(2) && Event
Z = fzero(@(Z) BLSlinQ(Tvalueslin(it+itT),Z,Q0),[-delta/2+infinit,Hmax]);
if Z <= -delta/2
    error('Unphysical solution: Z<=-delta/2');
end
dZ = (Z-Zold)/dtlin;       % Backward Euler (O(dt))
X(it,1)=Z; X(it,2)=dZ;
it=it+1; Zold=Z; Event=abs(Z)<abs(Hmin);
end
itE = it-1; 
TvaluesX = horzcat(TvaluesX,Tvalueslin(itB+itT:itE+itT));
X(itB:itE,3) = PinLIN(X(itB:itE,1));
X(itB:itE,4) = (X(itB:itE,3).*Va(X(itB:itE,1)))./(Rg*Temp);
% Solve diffusion equation 
if ~DiffusionApproximation
diffLin = @(t,Ca) A*Ca+B(interp1(Tvalueslin(itT+itB-1:itT+itE)',X(itB-1:itE,3)',t));
[~, Ca] = ode113(diffLin,Tvalueslin(itT+itB-1:itT+itE),iCa0');
if itB==itE
X(itB:itE,5:4+Nout) = Ca(end,:);
else
X(itB:itE,5:4+Nout) = Ca(2:end,:);
end
end
if ~Event
fitT = toc;
if DISPLAY == 2, disp(['Fsolve iteration finished in ' num2str(round(fitT,1)) 's']); end
end
end 

if ~Event&&~(Tvalueslin(itT+it)<=tspan(2))
    Pass = 1; PassPoint = 0;
end
if ~Event&&~Pass
reverseStr = '';    % String to tracks progress.
% A.2 Non-linear BLS equation
if DISPLAY==2, disp('Solve Rayleigh-Plesset equation'); end
if Checkpoint
if acorrpass 
tBLS = [acorrSPAN1 tspan(2)];
acorrpass = 0;
else
tBLS = [tspan(1) tspan(2)];
end
else
if PassPoint 
tBLS=[Tvalueslin(itE+itT) tspan(2)];
else
tBLS = [Tvalueslin(itB-1+itT) tspan(2)];
end
Checkpoint = 1;
end
W0 = [X(end,1:2), X(end,4:4+Nout)]'; % W = [Z dZ na Ca]
if DiffusionApproximation
%BLS = @(t,Z,dZ,na) BLS1Q(DISPLAY,tBLS,t,Z,dZ,na,R,rhol,PecQ,Q0,Pin,Pm,Po,PaT,omega,PS,delta0,mus,mul,S,Da,Ca,ka,ksi); %#ok<*UNRCH>
if SpeedUp == 1 || SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
EventFcn = @(t,W) EventFcn1(t,W(1),W(2),W(3),USfreq,dt,CORRTres);
OdeOpts = odeset('MaxStep',dt,'Events',EventFcn,'AbsTol',[1e-13;1e-9;1e-24],'RelTol',1e-4);
else
OdeOpts = odeset('MaxStep',dt);
end
[t,W] = ode113(@(t,W) BLS1Q(DISPLAY,tBLS,t,W(1),W(2),W(3),R,rhol,PecQ,Q0,Pin,Pm,Po,USPaT,omega,PS,delta0,mus,mul,S,Da,Ca,ka,ksi),tBLS,W0(1:3),OdeOpts);
else
%BLS = @(t,Z,dZ,na,Ca) BLS2Q(DISPLAY,tBLS,t,Z,dZ,na,Ca,R,rhol,PecQ,Q0,Pin,Pm,Po,PaT,omega,PS,delta0,mus,mul,S,Da,ka,A,B,deltaR);
if SpeedUp == 1 || SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
EventFcn = @(t,W) EventFcn2(t,W(1),W(2),W(3),W(4:3+Nout),USfreq,dt,CORRTres);
OdeOpts = odeset('MaxStep',dt,'Events',EventFcn);
else
OdeOpts = odeset('MaxStep',dt);    
end
[t,W] = ode113(@(t,W) BLS2Q(DISPLAY,tBLS,t,W(1),W(2),W(3),W(4:3+Nout),R,rhol,PecQ,Q0,Pin,Pm,Po,USPaT,omega,PS,delta0,mus,mul,S,Da,ka,A,B,deltaR),tBLS,W0,OdeOpts);
end
if t(end) == tBLS(end)
TvaluesX = horzcat(TvaluesX,t(2:end)'); %#ok<*AGROW>
InTb = size(X,1)+1; InTe = InTb+length(t)-2;
X(InTb:InTe,1:4)=[W(2:end,1:2), Pin(W(2:end,3),W(2:end,1)), W(2:end,3)]; 
if ~DiffusionApproximation
X(InTb:InTe,5:4+Nout)=W(2:end,4:3+Nout);
end
else
if DISPLAY==2
disp(' ');
disp(['Periodicity found in solution with autocorrelation >' num2str(CORRTres)]);
disp('Ending Rayleigh-Plesset solver...');
end
% We now fill the matrix X and TvaluesX by replicating the periodic
% function
% 1. W is (approximately with CORRTres) periodic with T=1/freq
% -> It is important to replicate EXACTLY the period, as small errors will
%    accumulate quickly
% -> Very important: the endconditions on X (or W) have to be conserved
%    over an integer times the period (otherwise discontinuities will
%    crash the program!)
% 2. length Dt = t(end)-t(firstIn) (approx<)= 1/freq
% 3. length to fill: tBLS(end)-t(end) 
tline = t-tBLS(1); Wline = W(:,1);
firstIn = find((t(end)-t)<(1/USfreq),1);
dtUNPER = (t(end)-t(1))-(1/USfreq);       % Timeperiod where BLS is unperiodic 
tPeriod = [(t(end)-1/USfreq), t(firstIn:end)'];
Ntperiod = tPeriod-t(1);
Wperiod = vertcat(W(end,:),W(firstIn:end,:));
Nfill = ceil((tBLS(end)-t(end))*USfreq);
TvaluesX = horzcat(TvaluesX,t(2:end)'); %#ok<*AGROW>
TvaluesX = horzcat(TvaluesX,...
    t(end)+cumsum(repmat(diff(tPeriod),1,Nfill)));
W = vertcat(W,repmat(W(firstIn:end,:),Nfill,1));
InTb = size(X,1)+1; InTe = InTb+size(W,1)-2;
X(InTb:InTe,1:4)=[W(2:end,1:2), Pin(W(2:end,3),W(2:end,1)), W(2:end,3)]; 
if ~DiffusionApproximation
X(InTb:InTe,5:4+Nout)=W(2:end,4:3+Nout);
end
acorrpass = 1;
acorrSPAN1 = TvaluesX(end);
acorrNSPAN2 = acorrSPAN1;
end
end
temp = [];
% B. Calculate V,m,n,p,h,s,u
% I. Event = 1: Tvalueslin(itB-1)->Tvalueslin(itE) (Y0 = Y(end,:))
% II. Event = 0, first time: Tvalueslin(itB-1)->tspan(2) (Y0=Y(end,:))
% III. Event = 0, tspan(1)->tspan(2) (Y0=Y(end,:))
if Event||Pass        % I.
tNICE = [Tvalueslin(itT+itB-1) Tvalueslin(itT+itE)];
if SpeedUp == 3 || SpeedUp == 6
Zinterp = @(t) interp1qr(Tvalueslin(itT+itB-1:itT+itE),X(itB-1:itE,1),t);
elseif SpeedUp == 4 || SpeedUp == 5
Zinterp = @(t) nakeinterp1(Tvalueslin(itT+itB-1:itT+itE),X(itB-1:itE,1),t);
elseif SpeedUp ~= 3 && SpeedUp ~= 4 && SpeedUp ~= 5 && SpeedUp ~= 6
Zinterp = @(t) interp1(Tvalueslin(itT+itB-1:itT+itE),X(itB-1:itE,1),t);
end
CmR = @(t) Cm(Zinterp(t));
else  % II. 
if ~Checkpoint2
tNICE = [Tvalueslin(itB+itT-1) tspan(2)];
if acorrpass
tNICE(2) = acorrNSPAN2;
acorrpass2 = 1;
acorrNSPAN1 = acorrNSPAN2;
end
if SpeedUp == 3 || (SpeedUp == 6 && acorrpass == 0)
Zinterp = @(t) interp1qr(TvaluesX(itB-1:end)',X(itB-1:end,1),t);
elseif SpeedUp == 4 || (SpeedUp == 5 && acorrpass == 0) 
Zinterp = @(t) nakeinterp1(TvaluesX(itB-1:end)',X(itB-1:end,1),t); 
elseif (SpeedUp == 5 && acorrpass == 1) || (SpeedUp == 6 && acorrpass == 1)
NNeg = sum(double(TvaluesX(itB-1:end)'-tBLS(1)<tline(1))); 
if NNeg >= 1
Wline = vertcat(X(itB-1:itB-2+NNeg,1),Wline);
tline = vertcat(TvaluesX(itB-1:itB-2+NNeg)'-tBLS(1),tline);
end
tMOD = @(t) (t-tBLS(1))*double((t-tBLS(1))<=dtUNPER)+...
    double((t-tBLS(1))>dtUNPER)*(mod(((t-tBLS(1))-dtUNPER),1/USfreq)+dtUNPER);
if SpeedUp == 5
Zinterp = @(t) nakeinterp1(tline,Wline,tMOD(t));
elseif SpeedUp == 6
Zinterp = @(t) interp1qr(tline,Wline,tMOD(t));
end
elseif SpeedUp ~= 3 && SpeedUp ~= 4 && SpeedUp ~= 5 && SpeedUp ~= 6
Zinterp = @(t) interp1(TvaluesX(itB-1:end),X(itB-1:end,1),t);
end
CmR = @(t) Cm(Zinterp(t));
Checkpoint2 = 1;
else % III.
tNICE = tspan;
if acorrpass2 
tNICE(1) = acorrNSPAN1;
acorrpass2 = 0;
end
if acorrpass
tNICE(2) = acorrNSPAN2;
acorrNSPAN1 = acorrNSPAN2;
acorrpass2 = 1;
end
INB = find(tNICE(1)==TvaluesX,1);
INE = find(tNICE(2)==TvaluesX,1);
if isempty(INB)||isempty(INE)
    error('Index empty - debug please');
end
if SpeedUp == 3 || (SpeedUp == 6 && acorrpass == 0) 
Zinterp =  @(t) interp1qr(TvaluesX(INB:INE)',X(end-(INE-INB):end,1),t);
elseif SpeedUp == 4 || (SpeedUp == 5 && acorrpass == 0)
Zinterp =  @(t) nakeinterp1(TvaluesX(INB:INE)',X(end-(INE-INB):end,1),t);  
elseif SpeedUp ~= 3 && SpeedUp ~= 4 && SpeedUp ~= 5 && SpeedUp ~= 6
Zinterp = @(t) interp1(TvaluesX(INB:INE),X(end-(INE-INB):end,1),t);
elseif SpeedUp == 5 && acorrpass == 1
tMOD = @(t) (t-tBLS(1))*double((t-tBLS(1))<=dtUNPER)+...
    double((t-tBLS(1))>dtUNPER)*(mod(((t-tBLS(1))-dtUNPER),1/USfreq)+dtUNPER);
Zinterp = @(t) nakeinterp1(tline,Wline,tMOD(t));
elseif SpeedUp == 6 && acorrpass == 1
tMOD = @(t) (t-tBLS(1))*double((t-tBLS(1))<=dtUNPER)+...
    double((t-tBLS(1))>dtUNPER)*(mod(((t-tBLS(1))-dtUNPER),1/USfreq)+dtUNPER);
Zinterp = @(t) interp1qr(tline,Wline,tMOD(t));
end
CmR = @(t) Cm(Zinterp(t));
end
end
if DISPLAY == 2
disp(' ');
disp('Solve generalized HH-type equations');
end
end
if ZONE == 1
CmR = @(t) Cm0;
dt = dtES; tNICE=tNICEc;
end
reverseStr = '';  
%--------------------------------------------------------------------------
% ---------------REGULAR OR FAST SPIKING NEURONS---------------------------
%--------------------------------------------------------------------------
if MODEL == 1 || MODEL == 2 || MODEL == 6 || MODEL == 7
    if SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
%NICE = @(t,Q,m,n,p,h)...
%    SimplNICERSLU(DISPLAY,tNICE,t,Q,m,n,p,h,CmR,Gna,Vna,Gk,Vk,Gm,Gl,Vl,am,ampbm,an,anpbn,pinf,taup,ah,ahpbh);
OdeOpts=odeset('MaxStep',dt,'AbsTol',1e-3,'RelTol',1e-3);
[t,U] = ode113(@(t,U) SimplNICERSFSLU(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),CmR,Gna,Vna,Gk,Vk,Gm,Gl,Vl,am,ampbm,an,anpbn,pinf,taup,ah,ahpbh,VLIMs),tNICE,U0,OdeOpts);
    else
%NICE = @(t,Q,m,n,p,h)...
%    SimplNICERS(DISPLAY,tNICE,t,Q,m,n,p,h,CmR,Gna,Vna,Gk,Vk,Gm,Gl,Vl,am,bm,an,bn,pinf,taup,ah,bh);
OdeOpts=odeset('MaxStep',dt,'AbsTol',1e-3,'RelTol',1e-3);
[t,U] = ode113(@(t,U) SimplNICERSFS(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),CmR,Gna,Vna,Gk,Vk,Gm,Gl,Vl,am,bm,an,bn,pinf,taup,ah,bh),tNICE,U0,OdeOpts);
    end
%--------------------------------------------------------------------------
%-------------------LOW THRESHOLD SPIKING NEURONS--------------------------
%--------------------------------------------------------------------------
elseif MODEL == 3 || MODEL == 8
    if SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
%NICE = @(t,Q,m,n,p,h,s,u)...
%    SimplNICELTSLU(DISPLAY,tNICE,t,Q,m,n,p,h,s,u,CmR,Gna,Vna,Gk,Vk,Gm,GT,VCa,Gl,Vl,am,ampbm,an,anpbn,pinf,taup,ah,ahpbh,sinf,taus,uinf,tauu);
OdeOpts=odeset('MaxStep',dt,'AbsTol',1e-3,'RelTol',1e-3);
[t,U] = ode113(@(t,U) SimplNICELTSLU(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),U(6),U(7),CmR,Gna,Vna,Gk,Vk,Gm,GT,VCa,Gl,Vl,am,ampbm,an,anpbn,pinf,taup,ah,ahpbh,sinf,taus,uinf,tauu,VLIMs),tNICE,U0,OdeOpts);
    else
%NICE = @(t,Q,m,n,p,h,s,u)...
%    SimplNICELTS(DISPLAY,tNICE,t,Q,m,n,p,h,s,u,CmR,Gna,Vna,Gk,Vk,Gm,GT,VCa,Gl,Vl,am,bm,an,bn,pinf,taup,ah,bh,sinf,taus,uinf,tauu);
OdeOpts=odeset('MaxStep',dt,'AbsTol',1e-3,'RelTol',1e-3);
[t,U] = ode113(@(t,U) SimplNICELTS(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),U(6),U(7),CmR,Gna,Vna,Gk,Vk,Gm,GT,VCa,Gl,Vl,am,bm,an,bn,pinf,taup,ah,bh,sinf,taus,uinf,tauu),tNICE,U0,OdeOpts);
    end
%--------------------------------------------------------------------------
%-------------------THALAMOCORTICAL NEURONS--------------------------------
%--------------------------------------------------------------------------
elseif MODEL == 4
    if SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
%NICE = @(t,Q,m,n,h,s,u,w,wLock,hProtein,cCai)...
%    SimplNICETCLU(DISPLAY,tNICE,t,Q,m,n,h,s,u,w,wLock,hProtein,cCai,...
%    CmR,Gna,Vna,Gk,Vk,GT,fVCa,Gl,Vl,GKL,Gh,ginc,Vh,am,ampbm,an,anpbn,ah,...
%    ahpbh,sinf,taus,uinf,tauu,winf,tauw,k1,k2,k3,k4,Far,deffCa,tauCa);
OdeOpts=odeset('MaxStep',dt);
[t,U] = ode15s(@(t,U) SimplNICETCLU(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),U(6),U(7),U(8),...
    U(9),U(10),CmR,Gna,Vna,Gk,Vk,GT,fVCa,Gl,Vl,GKL,Gh,ginc,Vh,am,ampbm,an,anpbn,ah,...
     ahpbh,sinf,taus,uinf,tauu,winf,tauw,k1,k2,k3,k4,Far,deffCa,tauCa,VLIMs),tNICE,U0,OdeOpts);
    else
%NICE = @(t,Q,m,n,p,h,s,u,w,wLock,hProtein,cCai)...
%    SimplNICELTC(DISPLAY,tNICE,t,Q,m,n,p,h,s,u,wLock,hProtein,cCai,...
%    CmR,Gna,Vna,Gk,Vk,GT,fVCa,Gl,Vl,GKL,Gh,ginc,Vh,am,bm,an,bn,ah,...
%    bh,sinf,taus,uinf,tauu,winf,tauw,k1,k2,k3,k4,Far,deffCa,tauCa)
OdeOpts=odeset('MaxStep',dt);
[t,U] = ode15s(@(t,U) SimplNICETC(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),U(6),U(7),U(8),...
    U(9),U(10),CmR,Gna,Vna,Gk,Vk,GT,fVCa,Gl,Vl,GKL,Gh,ginc,Vh,am,bm,an,bn,ah,...
     bh,sinf,taus,uinf,tauu,winf,tauw,k1,k2,k3,k4,Far,deffCa,tauCa),tNICE,U0,OdeOpts);
    end
%--------------------------------------------------------------------------
%------------------NUCLEUS RETICULARIS NEURONS-----------------------------
%--------------------------------------------------------------------------
elseif MODEL == 5
SparsityRE = sparse([1 1 1 1 1 1 1;
            1 1 0 0 0 0 0;
            1 0 1 0 0 0 0;
            1 0 0 1 0 0 0;
            1 0 0 0 1 0 0;
            1 0 0 0 0 1 0;
            1 0 0 0 1 1 1]);     % Sparsity matrix of Jacobian matrix (0 iff J=0)
    if SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
%NICE = @(t,Q,m,n,h,s,u,cCai)...
%    SimplNICERELU(DISPLAY,tNICE,t,Q,m,n,h,s,u,cCai,CmR,Gna,Vna,Gk,Vk,...
%   GT,fVCa,Gl,Vl,am,ampbm,an,anpbn,ah,...
%    ahpbh,sinf,taus,uinf,tauu,Far,deffCa,tauCa);
%OdeOpts=odeset('MaxStep',dt,'JPattern',SparsityRE);
OdeOpts = odeset('MaxStep',dt,'JPattern',SparsityRE);
[t,U] = ode15s(@(t,U) SimplNICERELU(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),U(6),U(7),...
    CmR,Gna,Vna,Gk,Vk,GT,fVCa,Gl,Vl,am,ampbm,an,anpbn,ah,ahpbh,sinf,taus,uinf,...
    tauu,Far,deffCa,tauCa,VLIMs),tNICE,U0,OdeOpts);
    else
%NICE = @(t,Q,m,n,h,s,u,cCai)...
%    SimplNICERE(DISPLAY,tNICE,t,Q,m,n,h,s,u,cCai,CmR,Gna,Vna,Gk,Vk,...
%   GT,fVCa,Gl,Vl,am,bm,an,bn,ah,...
%    bh,sinf,taus,uinf,tauu,Far,deffCa,tauCa);
OdeOpts=odeset('MaxStep',dt,'JPattern',SparsityRE);
[t,U] = ode15s(@(t,U) SimplNICERE(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),U(6),U(7),...
    CmR,Gna,Vna,Gk,Vk,GT,fVCa,Gl,Vl,am,bm,an,bn,ah,bh,sinf,taus,uinf,...
    tauu,Far,deffCa,tauCa),tNICE,U0,OdeOpts);
    end
% -------------------------------------------------------------------------
% -----------------------SUBTHALAMIC NUCLEUS MODEL-------------------------
% -------------------------------------------------------------------------
elseif MODEL == 9
    if SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
%NICE = @(t,Q,m,n,p,h)...
%    SimplNICESTNLU(DISPLAY,tNICE,t,Q,m,n,p,h,q,r,a,b,c,d1,d2,cCai,CmR,...
%    Gna,Vna,Gk,Vk,Gl,Vl,GT,fVCa,GCa,GA,GL,minf,ninf,pinf,hinf,qinf,rinf,ainf,...
%    binf,cinf,d1inf,d2inf,taum,taun,taup,tauh,tauq,taur,taua,taub,tauc,taud1,taud2,...
%    Far,tauCa)
if Charges
OdeOpts=odeset('MaxStep',dt);
[t,U] = ode113(@(t,U) SimplNICESTNLUcharges(ESi,DISPLAY,tNICE,t,...
    U(1),U(2),U(3),U(4),U(5),U(6),U(7),U(8),U(9),U(10),U(11),U(12),U(13),...
    U(14),U(15),U(16),U(17),U(18),U(19),U(20),...
    CmR,Gna,Vna,Gk,Vk,Gl,Vl,GT,fVCa,GCa,GA,GL,minf,ninf,...
    pinf,hinf,qinf,rinf,ainf,binf,cinf,d1inf,d2inf,taum,taun,taup,tauh,tauq,...
    taur,taua,taub,tauc,taud1,taud2,Far,tauCa,VLIMs,kcCaiLIMs),tNICE,U0,OdeOpts);
else
OdeOpts=odeset('MaxStep',dt,'AbsTol',1e-7,'RelTol',1e-4);
[t,U] = ode113(@(t,U) SimplNICESTNLU(ESi,DISPLAY,tNICE,t,...
    U(1),U(2),U(3),U(4),U(5),U(6),U(7),U(8),U(9),U(10),U(11),U(12),U(13),...
    CmR,Gna,Vna,Gk,Vk,Gl,Vl,GT,fVCa,GCa,GA,GL,minf,ninf,...
    pinf,hinf,qinf,rinf,ainf,binf,cinf,d1inf,d2inf,taum,taun,taup,tauh,tauq,...
    taur,taua,taub,tauc,taud1,taud2,Far,tauCa,VLIMs,kcCaiLIMs),tNICE,U0,OdeOpts);
end
    else
%NICE = @(t,Q,m,n,p,h)...
%    SimplNICESTN(DISPLAY,tNICE,t,Q,m,n,p,h,q,r,a,b,c,d1,d2,cCai,CmR,...§
%    Gna,Vna,Gk,Vk,Gl,Vl,GT,fVCa,GCa,GA,GL,minf,ninf,pinf,hinf,qinf,rinf,ainf,...
%    binf,cinf,d1inf,d2inf,taum,taun,taup,tauh,tauq,taur,taua,taub,tauc,taud1,taud2,...
%    Far,tauCa);
OdeOpts=odeset('MaxStep',dt);
[t,U] = ode113(@(t,U) SimplNICESTN(ESi,DISPLAY,tNICE,t,...
    U(1),U(2),U(3),U(4),U(5),U(6),U(7),U(8),U(9),U(10),U(11),U(12),U(13),...
    CmR,Gna,Vna,Gk,Vk,Gl,Vl,GT,fVCa,GCa,GA,GL,minf,ninf,pinf,hinf,qinf,rinf,ainf,...
    binf,cinf,d1inf,d2inf,taum,taun,taup,tauh,tauq,taur,taua,taub,tauc,taud1,taud2,...
    Far,tauCa),...
    tNICE,U0,OdeOpts);
    end
% -------------------------------------------------------------------------
% -----------------------THALAMUS RUBIN-TERMAN MODEL-----------------------
% -------------------------------------------------------------------------
elseif MODEL == 10
    if SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
%NICE = @(t,Q,m,n,p,h)...
%    SimplNICEThRTLU(DISPLAY,tNICE,t,Q,h,r,CmR,Gna,Vna,Gk,Vk,Gl,Vl,...
%    GT,VT,pinf,taup,rinf,taur);
OdeOpts=odeset('MaxStep',dt);
[t,U] = ode113(@(t,U) SimplNICEThRTLU(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),CmR,Gna,Vna,...
    Gk,Vk,Gl,Vl,GT,VT,hinf,tauh,rinf,taur,minf,pinf,VLIMs),tNICE,U0,OdeOpts);
    else
%NICE = @(t,Q,m,n,p,h)...
%    SimplNICEThRT(DISPLAY,tNICE,t,Q,h,r,CmR,Gna,Vna,Gk,Vk,Gl,Vl,...
%    GT,VT,pinf,taup,rinf,taur);
OdeOpts=odeset('MaxStep',dt);
[t,U] = ode113(@(t,U) SimplNICEThRT(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),CmR,Gna,Vna,Gk,Vk,Gl,Vl,...
    GT,VT,pinf,tauh,hinf,taur,minf,pinf),tNICE,U0,OdeOpts);
    end
% -------------------------------------------------------------------------
% ----------------------GLOBUS PALLIDUS INTERNUS NUCLEUS ------------------
% -------------------------------------------------------------------------
elseif MODEL == 11
    if SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
%NICE = @(t,Q,m,n,p,h)...
%    SimplNICEGPiLU(DISPLAY,tNICE,t,Q,n,h,r,CA,CmR,Gna,Vna,Gk,Vk,Gl,Vl,...
%    GT,VT,GCa,VCa,Gahp,minf,ainf,sinf,ninf,hinf,rinf);
OdeOpts=odeset('MaxStep',dt);
[t,U] = ode113(@(t,U) SimplNICEGPiLU(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),...
    CmR,Gna,Vna,Gk,Vk,Gl,Vl,GT,VT,GCa,VCa,Gahp,minf,ainf,sinf,ninf,hinf,rinf,...
    taun,tauh,taur,VLIMs),...
    tNICE,U0,OdeOpts);
    else
%NICE = @(t,Q,m,n,p,h)...
%    SimplNICEGPi(DISPLAY,tNICE,t,Q,n,h,r,CA,CmR,Gna,Vna,Gk,Vk,Gl,Vl,...
%    GT,VT,GCa,VCa,Gahp,minf,ainf,sinf,ninf,hinf,rinf);
OdeOpts=odeset('MaxStep',dt);
[t,U] = ode113(@(t,U) SimplNICEGPi(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),...
    CmR,Gna,Vna,Gk,Vk,Gl,Vl,GT,VT,GCa,VCa,Gahp,minf,ainf,sinf,ninf,hinf,rinf,...
    taun,tauh,taur),...
    tNICE,U0,OdeOpts);
    end
% -------------------------------------------------------------------------
% ----------------------GLOBUS PALLIDUS EXTERNUS NUCLEUS-------------------
% -------------------------------------------------------------------------
elseif MODEL == 12
    if SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
%NICE = @(t,Q,m,n,p,h)...
%    SimplNICEGPeLU(DISPLAY,tNICE,t,Q,n,h,r,CA,CmR,Gna,Vna,Gk,Vk,Gl,Vl,...
%    GT,VT,GCa,VCa,Gahp,minf,ainf,sinf,ninf,hinf,rinf);
OdeOpts=odeset('MaxStep',dt);
[t,U] = ode113(@(t,U) SimplNICEGPeLU(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),...
    CmR,Gna,Vna,Gk,Vk,Gl,Vl,GT,VT,GCa,VCa,Gahp,minf,ainf,sinf,ninf,hinf,rinf,...
    taun,tauh,taur,VLIMs),...
    tNICE,U0,OdeOpts);
    else
%NICE = @(t,Q,m,n,p,h)...
%    SimplNICEGPe(DISPLAY,tNICE,t,Q,n,h,r,CA,CmR,Gna,Vna,Gk,Vk,Gl,Vl,...
%    GT,VT,GCa,VCa,Gahp,minf,ainf,sinf,ninf,hinf,rinf);
OdeOpts=odeset('MaxStep',dt);
[t,U] = ode113(@(t,U) SimplNICEGPe(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),...
    CmR,Gna,Vna,Gk,Vk,Gl,Vl,GT,VT,GCa,VCa,Gahp,minf,ainf,sinf,ninf,hinf,rinf,...
    taun,tauh,taur),...
    tNICE,U0,OdeOpts);
    end
% -------------------------------------------------------------------------
% -----------------------MEDIUM SPINY STRIATUM NEURONS---------------------
% -------------------------------------------------------------------------
elseif MODEL == 13
    if SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
%NICE = @(t,Q,m,n,p,h)...
%    SimplNICEMSNLU(DISPLAY,tNICE,t,Q,m,n,p,h,CmR,Gna,Vna,Gk,Vk,Gl,Vl,Vm,...
%    am,ampbm,an,anpbn,ap,appbp,ah,ahpbh);
OdeOpts=odeset('MaxStep',dt);
[t,U] = ode113(@(t,U) SimplNICEMSNLU(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),...
    CmR,Gna,Vna,Gk,Vk,Gl,Vl,Vm,am,ampbm,an,anpbn,ap,appbp,ah,ahpbh,VLIMs),tNICE,U0,OdeOpts);
    else
%NICE = @(t,Q,m,n,p,h)...
%    SimplNICEMSN(DISPLAY,tNICE,t,Q,m,n,p,h,CmR,Gna,Vna,Gk,Vk,Gl,Vl,Vm,...
%    am,bm,an,bn,ap,bp,ah,bh);
OdeOpts=odeset('MaxStep',dt);
[t,U] = ode113(@(t,U) SimplNICEMSN(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),CmR,...
    Gna,Vna,Gk,Vk,Gl,Vl,Vm,am,bm,an,bn,ap,bp,ah,bh),tNICE,U0,OdeOpts);
    end
% -------------------------------------------------------------------------
% --------------------------HODGKIN-HUXLEY NEURONS-------------------------
% -------------------------------------------------------------------------
elseif MODEL == 14
    if SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
%NICE = @(t,Q,m,n,p,h)...
%    SimplNICEMSNLU(DISPLAY,tNICE,t,Q,m,n,p,h,CmR,Gna,Vna,Gk,Vk,Gl,Vl,Vm,...
%    am,ampbm,an,anpbn,ap,appbp,ah,ahpbh);
OdeOpts=odeset('MaxStep',dt,'AbsTol',1e-3,'RelTol',1e-3);
[t,U] = ode113(@(t,U) SimplNICEHHLU(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),...
    CmR,Gna,Vna,Gk,Vk,Gl,Vl,am,ampbm,an,anpbn,ah,ahpbh,VLIMs),tNICE,U0,OdeOpts);
    else
%NICE = @(t,Q,m,n,p,h)...
%    SimplNICEMSN(DISPLAY,tNICE,t,Q,m,n,p,h,CmR,Gna,Vna,Gk,Vk,Gl,Vl,Vm,...
%    am,bm,an,bn,ap,bp,ah,bh);
OdeOpts=odeset('MaxStep',dt,'AbsTol',1e-3,'RelTol',1e-3);
[t,U] = ode113(@(t,U) SimplNICEHH(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),CmR,...
    Gna,Vna,Gk,Vk,Gl,Vl,am,bm,an,bn,ah,bh),tNICE,U0,OdeOpts);
    end
end
TvaluesY = horzcat(TvaluesY,t(2:end)');
InTb = size(Y,1)+1; InTe = InTb+length(t)-2;
if MemorySaveMode == 3
Y(InTb:InTe,1) = U(2:end,1);

Y = interp1(TvaluesY',Y(:,1),(0:dtES:t(end))');
TvaluesY = (0:dtES:t(end));
else
Y(InTb:InTe,:) = U(2:end,:);
end
U0 = U(end,:)'; Q0 = U(end,1);
if size(Y,1) ~= size(TvaluesY,2)
    error('Error: size of Y-matrix and TvaluesY not equal');
end
if ZONE == 2
if Event||Pass
itB = it; iCa0 = Ca(end,:);
if MemorySaveMode 
itT = itB-2;
end
end
iZ0 = X(end,1);
Pass = 0;
if MemorySaveMode == 2
% Beware for possibility of making MemX unphysical when using non-linear
% interpolation!
[TvaluesXinterp,interpInd] = unique(TvaluesX(end-size(X,1)+1:end));
MemX=vertcat(MemX,interp1(TvaluesXinterp,X(interpInd,1),TvaluesY(InTb:InTe),'linear','extrap')');    
if iUP == 1, InTbb = InTb-1; end
if iUP == nUP, InTee = InTe; end    
if sum(X(:,1)<=-delta/2)
error('Unphysical solution in X');
end
if sum(MemX<=-delta/2)
error('Unphysical solution in MemX');
end
end
end
end
tCURRENT = t(end);
end
APindex = (10^5*Y(:,1)>Qthresh)&(circshift(10^5*Y(:,1),1,1)<Qthresh);
APindex(1) = 0; % Remove circshift artifact
APtimes = TvaluesY(APindex);  % AP times [s]
clear APindex; % Save all memory that can be saved...
NeuronActivated = ~isempty(APtimes); % Bool: 1 if neuron is activated
   switch MODEL
       case 1, MODELstr='RS';
       case 2, MODELstr='FS';
       case 3, MODELstr='LTS';
       case 4, MODELstr='TC';
       case 5, MODELstr='RE';
       case 6, MODELstr='RS_FVX';
       case 7, MODELstr='FS_FVX';
       case 8, MODELstr='LTS_CAX';
       case 9, MODELstr='STN';
       case 10, MODELstr='Th-RT';
       case 11, MODELstr='GPi';
       case 12, MODELstr='GPe';
       case 13, MODELstr='Str-MSN';
       case 14, MODELstr='HH';
   end
if MODE == 2
SaveStr=['APtimes(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USpstart) ',' num2str(USpd) ',' ...
        num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' USisppa ...
        ')-ES(' num2str(ESpstart) ',' num2str(ESpd) ',' num2str(ESdc) ',' ...
        num2str(ESprf) ',' ESisppa ').mat'];
save(SaveStr,'APtimes');

SaveStr2=['Chargevt(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USpstart) ',' num2str(USpd) ',' ...
        num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' USisppa ...
        ')-ES(' num2str(ESpstart) ',' num2str(ESpd) ',' num2str(ESdc) ',' ...
        num2str(ESprf) ',' ESisppa ').mat'];
saveChargeSample = [TvaluesY', 10^5*Y(:,1)];
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
Checkpoint=['CP(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USpstart) ',' num2str(USpd) ',' ...
        num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' USisppa ...
        ')-ES(' num2str(ESpstart) ',' num2str(ESpd) ',' num2str(ESdc) ',' ...
        num2str(ESprf) ',' ESisppa '):' num2str(SearchRange)];
disp(Checkpoint);   
end
end
if MODE == 1
SaveStr=['Thresh(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USpstart) ','...
    num2str(USpd) ',' num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' ...
    USisppa ')-ES(' num2str(ESpstart) ',' num2str(ESpd) ',' num2str(ESdc) ',' ...
        num2str(ESprf) ',' ESisppa ').mat'];    
save(SaveStr,'IIpa');
end
TTime = toc;
if DISPLAY
disp(' ');
disp(['Program finished in ' num2str(round(TTime,1)) 's']);
disp('Post-processing...');
end

% 5. Results
if PLOT
Charge = 10^5*Y(:,1); % Charge [nC/cm^2]
end
if PLOT == 2
SaveDataStr=['Data(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USpstart) ','...
num2str(USpd) ',' num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' ...
USisppa ')-ES(' num2str(ESpstart) ',' num2str(ESpd) ',' num2str(ESdc) ',' ...
    num2str(ESprf) ',' ESisppa ').mat'];
TvaluesYms = 10^(3)*TvaluesY'; % [ms]
saveData.TvaluesYms = TvaluesYms; saveData.Charge = Charge;
saveData.Y = Y; 
    if MemorySaveMode == 2 
    if exist('MemX') %#ok<*EXIST>
    CapacitanceUS = zeros(length(MemX),1);
    for i=1:length(MemX), CapacitanceUS(i) = Cm(MemX(i)); end
    Capacitance = [Cm0*ones(InTbb-1,1); CapacitanceUS; Cm0*ones(length(TvaluesYms)-InTee,1)];
    else
    Capacitance = Cm0;
    end
    Potential = 0.01*Charge./Capacitance; % [mV]
    saveData.Potential = Potential;
    saveData.Capacitance = Capacitance; %#ok<STRNU>
    end
save(SaveDataStr,'saveData','-v7.3');
end
if PLOT == 1
if MemorySaveMode == 3
TvaluesYms = 10^(3)*TvaluesY'; % [ms]
figure;
hold on;
plot(TvaluesYms,Charge);
ylabel('Charge [nC/cm^2]');
xlabel('Time [ms]');
set(gca,'fontsize',17);
hold off;
else
figure;
GateNames = {'m','n','p','h','q','r','a','b','c','d1','d2','Ca'};
nrGates = size(Y,2)-1; 
if Charges
    nrGates = nrGates-7;
end
TvaluesYms = 10^(3)*TvaluesY'; % [ms]

subplot(1+ceil(nrGates/2),1,1);
hold on;
plot(TvaluesYms,Charge);
ylabel('Charge [nC/cm^2]');
hold off;
for i=1:floor(nrGates/2)
subplot(1+ceil(nrGates/2),1,i+1);
hold on;
yyaxis left;
plot(TvaluesYms,Y(:,2*i));
ylabel(GateNames{2*i-1});
yyaxis right;
plot(TvaluesYms,Y(:,2*i+1));
ylabel(GateNames{2*i});
hold off;
set(gcf,'color','w');
end
if mod(nrGates,2)==1
subplot(1+ceil(nrGates/2),1,1+ceil(nrGates/2));
hold on;
plot(TvaluesYms,Y(:,end));
ylabel(GateNames{nrGates});
hold off;
end
xlabel('Time [ms]');
set(findall(gcf,'-property','FontSize'),'FontSize',14)
if MemorySaveMode == 2 && 1
figure;
if exist('MemX') %#ok<*EXIST>
CapacitanceUS = zeros(length(MemX),1);
for i=1:length(MemX), CapacitanceUS(i) = Cm(MemX(i)); end
Capacitance = [Cm0*ones(InTbb-1,1); CapacitanceUS; Cm0*ones(length(TvaluesYms)-InTee,1)];
else
Capacitance = Cm0;
end
Potential = 0.01*Charge./Capacitance; % [mV]
plot(TvaluesYms,Potential);
ylabel('Potential [mV]');
xlabel('Time [ms]');
set(gca,'fontsize',17);
set(gcf,'color','w'); box off;
if exist('MemX')
figure;
plot(TvaluesYms(InTbb:InTee),10^(9)*MemX);
xlabel('Time [ms]');
ylabel('Displacement [nm]')
set(gca,'fontsize',17);
set(gcf,'color','w'); box off;
end
figure;
plot(TvaluesYms,100*Capacitance);
xlabel('Time [ms]');
ylabel('Capacitance [\muF/cm^2]');
set(gca,'fontsize',17);
set(gcf,'color','w'); box off;
% ------------------ PLOT CURRENTS FOR SUBTHALAMIC NUCLEUS----------------
if MODEL == 9
mGate = Y(:,2); nGate = Y(:,3); pGate = Y(:,4); hGate = Y(:,5);
qGate = Y(:,6); rGate = Y(:,7); aGate = Y(:,8); bGate = Y(:,9); cGate = Y(:,10);
d1Gate = Y(:,11); d2Gate = Y(:,12); CaVal = Y(:,13);

iNaVal = 0.1*Gna*mGate.^3.*hGate.*(Potential-Vna); % Na-current [uA/cm^2]
iKVal = 0.1*Gk*nGate.^4.*(Potential-Vk);    % K-current [uA/cm^2]
iTVal = 0.1*GT*pGate.^2.*qGate.*(Potential-fVCa(CaVal)); % T-current [uA/cm^2]
ilVal = 0.1*Gl.*(Potential-Vl);         % Leak-current [uA/cm^2]
iKCaVal = 0.1*GCa.*rGate.^2.*(Potential-Vk); % Ca-dependent K-current [uA/cm^2]
iAVal = 0.1*GA*aGate.^2.*bGate.*(Potential-Vk); % A-type K-current [uA/cm^2]
iLVal = 0.1*GL*cGate.^2.*d1Gate.*d2Gate.*(Potential-fVCa(CaVal)); % L-current [uA/cm^2]
itotVal = 0.1*(Gl*(Potential-Vl)+Gna*mGate.^3.*hGate.*(Potential-Vna)+Gk*nGate.^4.*(Potential-Vk)+...
   GT*pGate.^2.*qGate.*(Potential-fVCa(CaVal))+GCa.*rGate.^2.*(Potential-Vk)+...
   GA*aGate.^2.*bGate.*(Potential-Vk)+GL*cGate.^2.*d1Gate.*d2Gate.*(Potential-fVCa(CaVal)));
% Total current [uA/cm^2]
if USpstart < Tsim && USpstart+USpd > 0 && Charges
    QNa = 10^5*Y(:,14);  % [nC/cm^2]
    QK = 10^5*Y(:,15);
    QT = 10^5*Y(:,16);
    QL = 10^5*Y(:,17);
    QKCa = 10^5*Y(:,18);
    QA = 10^5*Y(:,19);
    Ql = 10^5*Y(:,20);
BoolRange = TvaluesYms>=10^3*max(USpstart,0)&TvaluesYms<=10^3*min((USpstart+USpd),Tsim);
TvaluesYmsR = TvaluesYms(BoolRange);

QNa = QNa(BoolRange); QK = QK(BoolRange);
QT = QT(BoolRange); Ql = Ql(BoolRange);
QKCa = QKCa(BoolRange); QA = QA(BoolRange);
QL = QL(BoolRange); 
QNa = QNa-QNa(1); QK = QK-QK(1); QT = QT-QT(1);
Ql = Ql-Ql(1); QKCa = QKCa-QKCa(1); QA = QA-QA(1);
QL = QL-QL(1);
Qtot = QNa+QK+QT+Ql+QKCa+QA+QL;

% % Displaying all end charges to research mechanisms (uncomment)
% [QNa(end), QK(end), QT(end), Ql(end), QKCa(end), QA(end), ...
% QL(end), Qtot(end), QNa(end)+QK(end), ...
% QT(end)+Ql(end)+QL(end)+QA(end)+QKCa(end)]

figure;
subplot(4,2,1);
hold on;
plot(TvaluesYmsR,QNa);
ylabel('Q_{Na} [nC/cm^2]');
hold off;
subplot(4,2,2);
hold on;
plot(TvaluesYmsR,QK);
ylabel('Q_{K} [nC/cm^2]');
hold off;
subplot(4,2,3);
hold on;
plot(TvaluesYmsR,QT);
ylabel('Q_T [nC/cm^2]');
hold off;
subplot(4,2,4);
hold on;
plot(TvaluesYmsR,QL);
ylabel('Q_L [nC/cm^2]');
hold off;
subplot(4,2,5);
hold on;
plot(TvaluesYmsR,QKCa);
ylabel('Q_{K-Ca} [nC/cm^2]');
hold off;
subplot(4,2,6);
hold on;
plot(TvaluesYmsR,QA);
ylabel('Q_A [nC/cm^2]');
hold off;
subplot(4,2,7);
hold on;
plot(TvaluesYmsR,Ql);
ylabel('Q_l [nC/cm^2]');
xlabel('Time [ms]');
hold off;
subplot(4,2,8);
hold on;
plot(TvaluesYmsR,Qtot);
ylabel('Q_{tot} [nC/cm^2]');
xlabel('Time [ms]');
hold off;
set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',14);

figure;
subplot(1,2,1);
hold on;
plot(TvaluesYmsR,QNa+QK);
ylabel('Q_{Na}+Q_K [nC/cm^2]');
xlabel('Time [ms]');
hold off;
subplot(1,2,2);
hold on;
plot(TvaluesYmsR,QT+Ql+QKCa+QA+QL);
ylabel('Q_T+Q_l+Q_{KCa}+Q_A+Q_L [nC/cm^2]');
xlabel('Time [ms]');
hold off;
set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',17);
end

figure;
hold on;
subplot(2,2,1);
yyaxis left;
plot(TvaluesYms,iNaVal);
ylabel('i_{Na} [\muA/cm^2]');
yyaxis right;
plot(TvaluesYms,iKVal);
ylabel('i_{K} [\muA/cm^2]');
subplot(2,2,2);
yyaxis left;
plot(TvaluesYms,iTVal);
ylabel('i_T [\muA/cm^2]');
yyaxis right;
plot(TvaluesYms,iLVal);
ylabel('i_L [\muA/cm^2]');
subplot(2,2,3);
yyaxis left;
plot(TvaluesYms,iKCaVal);
ylabel('i_{K-Ca} [\muA/cm^2]');
yyaxis right;
plot(TvaluesYms,iAVal);
ylabel('i_A [\muA/cm^2]');
xlabel('Time [ms]');
subplot(2,2,4);
yyaxis left;
plot(TvaluesYms,ilVal);
ylabel('i_l [\muA/cm^2]');
yyaxis right;
plot(TvaluesYms,itotVal);
ylabel('i_{tot} [\muA/cm^2]');
xlabel('Time [ms]');
set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',14)

figure;
subplot(4,2,1);
hold on;
plot(TvaluesYms,iNaVal);
ylabel('i_{Na} [\muA/cm^2]');
hold off;
subplot(4,2,2);
hold on;
plot(TvaluesYms,iKVal);
ylabel('i_{K} [\muA/cm^2]');
hold off;
subplot(4,2,3);
hold on;
plot(TvaluesYms,iTVal);
ylabel('i_T [\muA/cm^2]');
hold off;
subplot(4,2,4);
hold on;
plot(TvaluesYms,iLVal);
ylabel('i_L [\muA/cm^2]');
hold off;
subplot(4,2,5);
hold on;
plot(TvaluesYms,iKCaVal);
ylabel('i_{K-Ca} [\muA/cm^2]');
hold off;
subplot(4,2,6);
hold on;
plot(TvaluesYms,iAVal);
ylabel('i_A [\muA/cm^2]');
hold off;
subplot(4,2,7);
hold on;
plot(TvaluesYms,ilVal);
ylabel('i_l [\muA/cm^2]');
xlabel('Time [ms]');
hold off;
subplot(4,2,8);
hold on;
plot(TvaluesYms,itotVal);
ylabel('i_{tot} [\muA/cm^2]');
xlabel('Time [ms]');
hold off;
set(gcf,'color','w');
set(findall(gcf,'-property','FontSize'),'FontSize',14)
end
end

% % - To plot charge from already plotted currents --
% axes = get(gcf,'children');
% for i = 1:length(axes)
% datastruct(i) = get(axes(length(axes)-i+1),'children');
% end
% for i = 1:length(datastruct)
% xdata{i} = get(datastruct(i),'XData');
% ydata{i} = get(datastruct(i),'YData');
% end
% ydata2 = cellfun(@(ydata,xdata) ydata(xdata>=1000&xdata<=2000),ydata,xdata,'UniformOutput',false);
% xdata2 = cellfun(@(xdata) xdata(xdata>=1000&xdata<=2000),xdata,'UniformOutput',false);
% Dxdata = cellfun(@(xdata) circshift(xdata,-1,2)-xdata,xdata2,'UniformOutput',false);
% Dxdata = cellfun(@(Dxdata) Dxdata(1:end-1),Dxdata,'UniformOutput',false);
% charge = cellfun(@(ydata2,Dxdata) [0 cumsum(ydata2(1:end-1).*Dxdata)],ydata2,Dxdata,'UniformOutput',false);
% figure;
% for i = 1:length(xdata2)
% subplot(4,2,i);
% plot(xdata2{i},charge{i});
% end

% ------------------ PLOT CURRENTS FOR HH MODEL --------------------------
if MODEL == 14
mGate = Y(:,2); nGate = Y(:,3); hGate = Y(:,4); 
iNaVal = 0.1*Gna*mGate.^3.*hGate.*(Potential-Vna); % Na-current [uA/cm^2]
iKVal = 0.1*Gk*nGate.^4.*(Potential-Vk);    % K-current [uA/cm^2]
ilVal = 0.1*Gl.*(Potential-Vl);         % Leak-current [uA/cm^2]
itotVal = iNaVal+iKVal+ilVal;           % Total current [uA/cm^2]
figure; hold on;
subplot(2,2,1);
plot(TvaluesYms,iNaVal);
ylabel('i_{Na} [\muA/cm^2]'); xlabel('time [ms]');
subplot(2,2,2);
plot(TvaluesYms,iKVal);
ylabel('i_{K} [\muA/cm^2]'); xlabel('time [ms]');
subplot(2,2,3);
plot(TvaluesYms,ilVal);
ylabel('i_l [\muA/cm^2]'); xlabel('time [ms]');
subplot(2,2,4);
plot(TvaluesYms,itotVal); xlabel('time [ms]');
ylabel('i_{tot} [\muA/cm^2]'); set(gcf,'color','w');
title(['Total membrane current. Injected current (' num2str(100*ESipa) ' \muA/cm^2) excluded']);
set(findall(gcf,'-property','FontSize'),'FontSize',14);
hold off;
end
end
end
end
