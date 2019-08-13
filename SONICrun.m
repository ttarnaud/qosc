function SONICrun(Tsim,MODE,USpstart,USpd,USfreq,USdc,USprf,USisppa,ESpstart,ESpd,...
    ESdc,ESprf,ESisppa,PLOT,model,USibegin,USiend,SearchMode)
coder.extrinsic('nakeinterp1');
% SONICrun is a speed-up version of funPES based on multi-scale
% optimization with effective parameters. Based on:
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

DISPLAY = 1;
% Display level. Note: higher display level will give more runtime information but will slow the program 
% DISPLAY = 0 -> No information displayed (use this option for HPC simulations)
% DISPLAY = 1 -> Display progress based on update nr. alone

Charges = 0; % This mode is meant for investigation on the mechanism of US-STN stimulation:
% All charges are calculated, but the solver is slowed down
if Charges == 1
    disp('Note 3: Charge-mode is on for STN-mechanism investigation!');
    disp(' ');
end

% Discretisation parameters
dt = 50e-6;            % Discretisation time (s)
atol = 1e-6; rtol = 1e-3; % absolute and relative VSVO-tolerances

tic;
% 1. Parameters
Qthresh = 0;            % Threshold for Q for AP-discrimination [nC/cm^2]
Cm0 = 0.01;				% Rest capacitance (F/m^2)	 
a=32*10^(-9);			% radius leaflet boundary (m)
c = 1515;				% Speed of sound surrounding medium (m/s)
rhol = 1028;			% Density surrounding medium (kg/m^3)
USPa = sqrt(2*rhol*c*USipa);
Rg = 8.314;             % Universal gas constant (J/(K*mol))
Far = 96485.3329;       % Faraday constant (C/mol) 
%Temp = 309.15; 		    % Surrounding medium temperature (K)	
%Temp = 306.15;
Temp = 36+273.15;
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
modelName = 'HH';
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

% 2. Important functions
SONIC = load(['SONIC-' modelName '.mat']); 
SONICtable = SONIC.SONICtable;
% 2.1 SONIC functions (rate, Veff, Zeff, Cmeff, ngend)
QmRange = SONICtable.QmRange; USPaRange = SONICtable.USPaRange; 
USfreqRange = SONICtable.USfreqRange; aBLSRange = SONICtable.aBLSRange;

Veff4D = SONICtable.Veff; Zeff4D = SONICtable.Zeff; Cmeff4D = SONICtable.Cmeff; ngend4D = SONICtable.ngend;

% rate 4D sonic tables
SONICfields = fieldnames(SONICtable);
SONICrates = SONICfields(cellfun(@(X) contains(X,'a_')|contains(X,'apb_'),SONICfields));
SONICgates = cellfun(@(X) X(3:end),SONICrates(cellfun(@(X) contains(X,'a_'),SONICrates)),'UniformOutput',0); 
rt = struct;
for i = 1:length(SONICrates)
rt.(SONICrates{i}) = SONICtable.(SONICrates{i});
end
f4Veff = @(Qm,USPa,USfreq,aBLS) interpn(QmRange,USPaRange,USfreqRange,aBLSRange,Veff4D,Qm,USPa,USfreq,aBLS,'linear');
f4Zeff = @(Qm,USPa,USfreq,aBLS) interpn(QmRange,USPaRange,USfreqRange,aBLSRange,Zeff4D,Qm,USPa,USfreq,aBLS,'linear');
f4Cmeff = @(Qm,USPa,USfreq,aBLS) interpn(QmRange,USPaRange,USfreqRange,aBLSRange,Cmeff4D,Qm,USPa,USfreq,aBLS,'linear');
f4ngend = @(Qm,USPa,USfreq,aBLS) interpn(QmRange,USPaRange,USfreqRange,aBLSRange,ngend4D,Qm,USPa,USfreq,aBLS,'linear');

f3Veff = @(Qm,USPa,USfreq) f4Veff(Qm,USPa,USfreq,a); 
f3Zeff = @(Qm,USPa,USfreq) f4Zeff(Qm,USPa,USfreq,a); %#ok<*NASGU>
f3Cmeff = @(Qm,USPa,USfreq) f4Cmeff(Qm,USPa,USfreq,a);
f3ngend = @(Qm,USPa,USfreq) f4ngend(Qm,USPa,USfreq,a);

f4rt = struct; f3rt = struct;
for i = 1:length(SONICrates)
f4rt.(SONICrates{i}) =  @(Qm,USPa,USfreq,aBLS) interpn(QmRange,USPaRange,USfreqRange,aBLSRange,rt.(SONICrates{i}),Qm,USPa,USfreq,aBLS,'linear');   
f3rt.(SONICrates{i}) = @(Qm,USPa,USfreq) f4rt.(SONICrates{i})(Qm,USPa,USfreq,a); 
end


VecVeff0 = zeros(1,length(QmRange)); VecVeffPa = zeros(1,length(QmRange));
Vecrt0 = struct; VecrtPa = struct; f1rt0 = struct; f1rtPa = struct;
for i = 1:length(QmRange)
VecVeff0(i) = f3Veff(QmRange(i),0,USfreq);
VecVeffPa(i) = f3Veff(QmRange(i),USPa,USfreq);
end
for j = 1:length(SONICrates)
Vecrt0.(SONICrates{j}) = zeros(1,length(QmRange));
VecrtPa.(SONICrates{j}) = zeros(1,length(QmRange));
for i = 1:length(QmRange)
Vecrt0.(SONICrates{j})(i) = f3rt.(SONICrates{j})(QmRange(i),0,USfreq);
VecrtPa.(SONICrates{j})(i) = f3rt.(SONICrates{j})(QmRange(i),USPa,USfreq);
end
f1rt0.(SONICrates{j}) = @(Q) nakeinterp1(QmRange',Vecrt0.(SONICrates{j}),Q);
f1rtPa.(SONICrates{j}) = @(Q) nakeinterp1(QmRange',VecrtPa.(SONICrates{j}),Q);
end
f1Veff0 = @(Q) nakeinterp1(QmRange',VecVeff0,Q);
f1VeffPa = @(Q) nakeinterp1(QmRange',VecVeffPa,Q);


fVCa = @(cCai) 10^(3)*((Rg*Temp)/(2*Far))*log(cCae./cCai); % Nernst equation for Ca-potential [mV] (if not assumed constant)

% 3. Initial conditions and timespan
% (order:) Y = [Q,SONICgates, other gates [d20,r0,w0], wLock, hProtein, cCai]  (other gates are cCai/hprotein dependent)
Qm0 = Cm0*(10^(-3)*Vm0);            % Approximation of rest charge
Y0 = zeros(round(length(SONICrates)/2)+1,1);
Y0(1) = Qm0;
for i = 1:length(SONICgates)
Y0(i+1) = f3rt.(['a_' SONICgates{i}])(Qm0,0,1e6)./f3rt.(['apb_' SONICgates{i}])(Qm0,0,1e6);      
end

if MODEL == 4 || MODEL == 5
cCai0 = 240*10^(-6);             % Rest concentration Ca (mol/m^3)
Y0 = horzcat(Y0,cCai0);
elseif MODEL == 9
cCai0 = 5*10^(-6);
d20 = d2inf(cCai0);
r0 = rinf(cCai0);
Y0 = horzcat(Y0,[d20,r0,cCai0]);
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
w0 = winf(Vm0)/(1+winf(Vm0)*((k3*hProtein0)/k4));
wLock0 = (k3/k4)*(w0*hProtein0);   % Initial condition for locked w-gate
Y0 = horzcat(Y0,[w0,wLock0,hProtein0]);
end


% 4. Solver
global reverseStr; reverseStr = '';  %#ok<*TLEV>
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
OdeOpts=odeset('MaxStep',dt,'AbsTol',atol,'RelTol',rtol); tNICE = [0,Tsim];
%--------------------------------------------------------------------------
% ---------------REGULAR OR FAST SPIKING NEURONS---------------------------
%--------------------------------------------------------------------------
if MODEL == 1 || MODEL == 2 || MODEL == 6 || MODEL == 7
[t,U] = ode113(@(t,U) SONIC_RSFS(ESi,USPaT,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),...
    Gna,Vna,Gk,Vk,Gm,Gl,Vl,f1Veff0,f1VeffPa,f1rt0,f1rtPa,SONICgates),tNICE,Y0,OdeOpts);
%--------------------------------------------------------------------------
%-------------------LOW THRESHOLD SPIKING NEURONS--------------------------
%--------------------------------------------------------------------------
elseif MODEL == 3 || MODEL == 8
[t,U] = ode113(@(t,U) SONIC_LTS(ESi,USPaT,USfreq,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),U(6),U(7),...
    Gna,Vna,Gk,Vk,Gm,GT,VCa,Gl,Vl,f1Veff0,f1VeffPa,f1rt0,f1rtPa,SONICgates),tNICE,Y0,OdeOpts);
%--------------------------------------------------------------------------
%-------------------THALAMOCORTICAL NEURONS--------------------------------
%--------------------------------------------------------------------------
elseif MODEL == 4
[t,U] = ode15s(@(t,U) SONIC_TC(ESi,USPaT,USfreq,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),U(6),U(7),U(8),...
    U(9),U(10),Gna,Vna,Gk,Vk,GT,fVCa,Gl,Vl,GKL,Gh,ginc,Vh,k1,k2,k3,k4,Far,deffCa,tauCa,f3Veff,f3rt),tNICE,Y0,OdeOpts);
%--------------------------------------------------------------------------
%------------------NUCLEUS RETICULARIS NEURONS-----------------------------
%--------------------------------------------------------------------------
elseif MODEL == 5
[t,U] = ode15s(@(t,U) SONIC_RE(ESi,USPaT,USfreq,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),U(6),U(7),...
    Gna,Vna,Gk,Vk,GT,fVCa,Gl,Vl,Far,deffCa,tauCa,f3Veff,f3rt),tNICE,U0,OdeOpts);
% -------------------------------------------------------------------------
% -----------------------SUBTHALAMIC NUCLEUS MODEL-------------------------
% -------------------------------------------------------------------------
elseif MODEL == 9
if Charges
[t,U] = ode113(@(t,U) SimplNICESTNLUcharges(ESi,DISPLAY,tNICE,t,...
    U(1),U(2),U(3),U(4),U(5),U(6),U(7),U(8),U(9),U(10),U(11),U(12),U(13),...
    U(14),U(15),U(16),U(17),U(18),U(19),U(20),...
    CmR,Gna,Vna,Gk,Vk,Gl,Vl,GT,fVCa,GCa,GA,GL,minf,ninf,...
    pinf,hinf,qinf,rinf,ainf,binf,cinf,d1inf,d2inf,taum,taun,taup,tauh,tauq,...
    taur,taua,taub,tauc,taud1,taud2,Far,tauCa,VLIMs,kcCaiLIMs),tNICE,U0,OdeOpts); %#ok<*UNRCH>
else
OdeOpts=odeset('MaxStep',dt,'AbsTol',1e-7,'RelTol',1e-4);
[t,U] = ode113(@(t,U) SimplNICESTNLU(ESi,DISPLAY,tNICE,t,...
    U(1),U(2),U(3),U(4),U(5),U(6),U(7),U(8),U(9),U(10),U(11),U(12),U(13),...
    CmR,Gna,Vna,Gk,Vk,Gl,Vl,GT,fVCa,GCa,GA,GL,minf,ninf,...
    pinf,hinf,qinf,rinf,ainf,binf,cinf,d1inf,d2inf,taum,taun,taup,tauh,tauq,...
    taur,taua,taub,tauc,taud1,taud2,Far,tauCa,VLIMs,kcCaiLIMs),tNICE,U0,OdeOpts);
end
% -------------------------------------------------------------------------
% -----------------------THALAMUS RUBIN-TERMAN MODEL-----------------------
% -------------------------------------------------------------------------
elseif MODEL == 10
[t,U] = ode113(@(t,U) SONIC_ThRT(ESi,USPaT,USfreq,DISPLAY,tNICE,t,U(1),U(2),U(3),Gna,Vna,...
    Gk,Vk,Gl,Vl,GT,VT,f3Veff,f3rt),tNICE,U0,OdeOpts);
% -------------------------------------------------------------------------
% ----------------------GLOBUS PALLIDUS INTERNUS NUCLEUS ------------------
% -------------------------------------------------------------------------
elseif MODEL == 11
[t,U] = ode113(@(t,U) SimplNICEGPiLU(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),...
    CmR,Gna,Vna,Gk,Vk,Gl,Vl,GT,VT,GCa,VCa,Gahp,minf,ainf,sinf,ninf,hinf,rinf,...
    taun,tauh,taur,VLIMs),...
    tNICE,U0,OdeOpts);
% -------------------------------------------------------------------------
% ----------------------GLOBUS PALLIDUS EXTERNUS NUCLEUS-------------------
% -------------------------------------------------------------------------
elseif MODEL == 12
[t,U] = ode113(@(t,U) SimplNICEGPeLU(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),...
    CmR,Gna,Vna,Gk,Vk,Gl,Vl,GT,VT,GCa,VCa,Gahp,minf,ainf,sinf,ninf,hinf,rinf,...
    taun,tauh,taur,VLIMs),...
    tNICE,U0,OdeOpts);
% -------------------------------------------------------------------------
% -----------------------MEDIUM SPINY STRIATUM NEURONS---------------------
% -------------------------------------------------------------------------
elseif MODEL == 13
[t,U] = ode113(@(t,U) SimplNICEMSNLU(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),U(5),...
    CmR,Gna,Vna,Gk,Vk,Gl,Vl,Vm,am,ampbm,an,anpbn,ap,appbp,ah,ahpbh,VLIMs),tNICE,U0,OdeOpts);
% -------------------------------------------------------------------------
% --------------------------HODGKIN-HUXLEY NEURONS-------------------------
% -------------------------------------------------------------------------
elseif MODEL == 14
[t,U] = ode113(@(t,U) SimplNICEHHLU(ESi,DISPLAY,tNICE,t,U(1),U(2),U(3),U(4),...
    CmR,Gna,Vna,Gk,Vk,Gl,Vl,am,ampbm,an,anpbn,ah,ahpbh,VLIMs),tNICE,U0,OdeOpts);
end
TvaluesY = t;

APindex = (10^5*U(:,1)>Qthresh)&(circshift(10^5*U(:,1),1,1)<Qthresh);
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

