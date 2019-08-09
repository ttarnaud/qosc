function [Zeff,Veff,a_i_eff,apb_i_eff,ngend,Cmeff] = SONICcalc(MODEL,Qm,USPa,USfreq,aBLS)
coder.extrinsic('nakeinterp1');
% SONICcalc calculates effective variables for one set of stimulus parameters

%% Input:
% 1) MODEL 
% Choose model: 1=RS,2=FS,3=LTS,4=TC,5=RE,6=RS_FVX,7=FS_FVX,8=LTS_CAX,
% 9=STN, 10=Th-RT, 11=GPi, 12=GPe, 13=Str-MSN, 14 = HH

% RS,FS,LTS,RS_FVX,FS_FVX and LTS_CAX are Pospischil models
% RE and TC are Destexhe-based models
% Th-RT, GPi, GPe are based on Rubin-Terman
% STN is based on Otsuka et al.
% Str-MSN is based on McCarthy et al.

% 2) Qm : charge on membrane [C/m^2]
% 3) USPa : instantaneous local ultrasonic pressure [Pa]
% 4) USfreq : ultrasonic local frequency [Hz]
% 5) aBLS : radius of bilayer sonophore [m]

%% Output:
% varargout contains: [Zeff,Veff,alpha_i_eff,beta_i_eff,ngend]
% 1) Zeff : effective (i.e. time averaged) displacement of Z [m]
% 2) Veff: effective membrane potential [mV]
% 3) alpha_i_eff : struct with effective rate constants alpha [/s]
% 4) alphapbeta_i_eff: struct with effective rate constants alpha+beta [/s]  
% 5) ngend: gas constant at end of ultrasonic cycle [mol]

%% Based on:
% SONIC (Lemaire et al., 2018), 
% Piezoelectric solver (see Plaksin et al.,2014; Plaksin et al.,2016;
% Krasovitski et al., 2011, PhD thesis Plaksin, 2016, Tarnaud et al.,2018).

% Remark:
% The constant ks (areal modulus) is taken as the ratio of the total
% tension to the areal strain. The tension T of one leaflet is taken as
% half the total tension.

% Remark 2: Speed considerations
% SpeedUp -> NICE-model, UpSpeed -> BLS-model
SpeedUp = 5;
UpSpeed = 4;
CORRTres = 0.99;     % Threshold for normalized unbiased autocorrelation in
PerCheckMax = 4;    % At certain unpredictable parameters-values, presumed due to non-linearities, the BLS equations
                     % will never reach a periodic regime with f=USfreq. In this case, the solution can be periodic with a harmonic undertone
                     % of f=USfreq (e.g. with f2=USfreq/2) or never be truely periodic (e.g. some rough periodicity with f=USfreq can be observed, but each
                     % period has a different amplitude). 
                     % PerCheckMax is the maximal number of considered
                     % periods, in the autocorrelation (PerCheckMax should
                     % be >= 2). EventFcn will terminate BLS1Q, if there is periodicity with T in [1/USfreq,(PerCheckMax/2)*1/USfreq]
                     % OR f in [2*USfreq/PerCheckMax,USfreq]. 
NonPerEffP = 10; % If (Z, na) is never fully periodic, effective parameters are calculated by integration over [t(end)-NonPerEffP/USfreq,t(end)]
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
% Upspeed = 4 Lennard-Jones approximation of intramolecular integral (see Lemaire et al., 2018)
if SpeedUp < 1
error('SpeedUp has to be at least 1 in the SONIC implementation');
end

DISPLAY = 0;
% Display level. Note: higher display level will give more runtime information but will slow the program 
% DISPLAY = 0 -> No information displayed (use this option for HPC simulations)
% DISPLAY = 1 -> display progress level: low
% DISPLAY = 2 -> display progress level: high

ERROR = 1; % If 1, show error for unphysical situation Z<=-delta/2, if 0 extrapolate value for Pa-r
% (Use 0 for HPC simulations) 
if DISPLAY
if ERROR == 0
    disp('Note 2: solver will not check for unphysical molecular pressure');
    disp(' ');
elseif ERROR == 1
    disp('Note 2: solver will check for unphysical molecular pressure');
    disp(' ');
end
end

tic;
DiffusionApproximation = 1; % Bool: if 1, apply approximation (*) on diffusion of air.
%*Krasovitski 2011 doesn't use this approximation, Plaksin 2014 does use it. 

% 1. Parameters
% Delta0 values are calculated by:
% Pec(Vm0,0)/Ar = (deltax/delta)^x-(deltax/delta)^y
% FunDelta0 = @(delta)(deltax./delta).^x-(deltax./delta).^y-(0.01*0.001*Vm0)^2/(2*epsr*eps0*Ar);
% plot((1e-15:1e-11:50e-9),FunDelta0(1e-15:1e-11:50e-9));
% fzero(FunDelta0,[1e-15,50*10^(-9)])
xfs = 1;                % Active area fraction (-) within LU-table
ksi = 0.5*10^(-9);		% Boundary layer length leaflet-medium (m)
deltaR = ksi;           % discretisation step outside the BLS (m)
Rsim = 100*ksi;         % Total length simulation domain outside BLS (m)
eps0 = 8.854187817*10^(-12);% absolute vacuumpermittivity (F/m) 
kappa = 1;              % Polytropic constant (-)
epsr = 1;				% Relative permittivity intramembrane cavity (-)
rhol = 1028;			% Density surrounding medium (kg/m^3)
Po = 10^5;				% Static pressure surrounding medium (Pa)
omega = 2*pi*USfreq;  	% Radial frequency (Hz)
ks = 0.24; 				% Areal modulus of bilayer membrane (N/m)
delta0=2*10^(-9);		% Thickness leaflet (m)
mus = 0.035;			% Dynamic viscosity leaflets (Pa*s)
mul = 0.7*10^(-3);		% Dynamic viscosity surrounding medium (Pa*s)	
Cm0 = 0.01;				% Rest capacitance (F/m^2)	 
Rg = 8.314;             % Universal gas constant (J/(K*mol))
%Far = 96485.3329;      % Faraday constant (C/mol) 
%Temp = 309.15; 		% Surrounding medium temperature (K)	
%Temp = 306.15;
Temp = 36+273.15;
Da = 3*10^(-9);			% Diffusion coefficient of air in surrounding medium (m^2/s)
Ca = 0.62;				% Molar air concentration external medium (mol/m^3)
ka = 1.613*10^5;	    % Henry constant dissolved air external medium (Pa*m^3/mol) 
Ar = 10^5;				% Attraction/repulsion pressure coefficient (Pa) 
deltax = 1.4*10^(-9); 	% Initial gap between leaflets (no charge) (m) 
x = 5;					% Repulsion exponent
y = 3.3;				% Attraction exponent
Nout = floor(Rsim/deltaR);% Number of steps outside of BLS
if MODEL==1             % Regular spiking neuron
VT = -56.2;             % Spike treshold adjustment parameter (mV)
taumax = 608*10^(-3);   % Decay time constant slow non-inactivating K+ (s)
Vm0 = -71.9;            % rest membrane potential (mV)
delta = 1.26*10^(-9);   % Initial gap between leaflets (charge) (m)
elseif MODEL == 2       % Fast spiking neuron
VT = -57.9;             % Spike treshold adjustment parameter (mV)
taumax = 502*10^(-3);   % Decay time constant slow non-inactivating K+ (s)
Vm0 = -71.4;            % Rest membrane potential (mV)
delta = 1.26*10^(-9);   % Initial gap between leaflet (charge) (m)  
elseif MODEL == 3       % Low threshold spiking neuron
VT = -50;               % Spike treshold adjustment parameter (mV)
taumax = 4;             % Decay time constant slow non-inactivating K+ (s)
Vm0 = -54;              % Rest membrane potential (mV)
Vx = -7;				% Shift Ca2+ voltage (mV)
delta = 1.3*10^(-9);    % Initial gap between leaflets (charge) (m)
elseif MODEL == 4       % Thalamocortical neuron
VT = -52;               % Spike treshold adjustment parameter (mV)
Vm0 = -63.4;            % Rest membrane potential (mV)
Vx = 0;                 % Shift Ca2+ voltage (mV)
delta = 1.28*10^(-9);   % Initial gap between leaflet (charge) (m)  
elseif MODEL == 5       % Reticular thalamus neuron
VT = -67;               % Spike treshold adjustment parameter (mV)
Vm0 = -89.5;            % Rest membrane potential (mV)
delta = 1.21*10^(-9);   % Initial gap between leaflet (charge) (m)
elseif MODEL == 6       % Regular spiking ferret visual cortex neuron
VT = -55;               % Spike treshold adjustment parameter (mV)
taumax = 1;             % Decay time constant slow non-inactivating K+ (s)
Vm0 = -70.4;            % Rest membrane potential (mV)
delta = 1.26*10^(-9);   % Initial gap between leaflet (charge) (m)
elseif MODEL == 7       % Fast spiking ferret visual cortex neuron
VT = -55;               % Spike treshold adjustment parameter (mV)
taumax = 1;             % Decay time constant slow non-inactivating K+ (s)
Vm0 = -70;              % Rest membrane potential (mV)
delta = 1.26*10^(-9);   % Initial gap between leaflet (charge) (m)    
elseif MODEL == 8       % Low threshold spiking cat association cortex neuron
VT = -55;               % Spike treshold adjustment parameter (mV)
taumax = 1;             % Decay time constant slow non-inactivating K+ (s)
Vm0 = -84.6;            % Rest membrane potential (mV)
Vx = -2;				% Shift Ca2+ voltage (mV)
delta = 1.3*10^(-9);    % Initial gap between leaflet (charge) (m)
elseif MODEL == 9       % Subthalamic nucleus neurons -> Otsuka model
Vm0=-58;                % Resting potential (mV)
delta=1.2923e-09;       % Initial gap between leaflets (charge) (m) (calculated from Vm0 defined supra)
elseif MODEL == 10      % Thalamic Rubin-Terman based neuron
VT=0;                   % Nernst potential of T-type low threshold Ca2+ channel (mV)
Vm0=-65;                % Resting potential (mV)
delta=1.2744e-09;       % Initial gap betweel leaflets (m)
elseif MODEL == 11      % Globus pallidus internus neuron
VT=0;                   % Nernst potential of T-type Ca2+ channel (mV) 
Vm0 = -65;              % Resting potential (mV)
delta=1.2713e-09;       % Initial gap between leaflets (charge) (m)
elseif MODEL == 12      % Globus pallidus externus neuron
VT=0;                   % Nernst potential of T-type Ca2+ channel (mV) 
Vm0 = -65;              % Resting potential (mV)
delta=1.2713e-09;       % Initial gap between leaflets (charge) (m)
elseif MODEL == 13      % Striatum medium spiny neuron
Vm0=-87;                % Resting potential (mV)
delta=1.2177e-09;       % Initial gap between leaflets (charge) (m)
elseif MODEL == 14      % Hodgking-Huxley model 
Vm0 = -70;              % Resting potential (mV)
Qam = 3; Qbm = 3; Qan = 3; Qbn = 3; Qah = 3; Qbh = 3;
Qap = 3; Qbp = 3;       %#ok<NASGU> % Temperature factors (-)
Temp0 = 6.3+273.15;     % Default temperature (K)
delta = 1.26e-9;        % Initial gap between leaflets (charge) (m)
end
global Hmin;
Hmin = 0.001*1*10^(-12);     % Linear approximation Rayleigh-Plesset
Hmax = 1*10^(-9);        % Hmax for fzero

a = aBLS;
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
    ((V-VT-13)./(double((V-VT-13)~=0).*exp(-((V-VT-13)/4))-1)));% Rate constant alpha_m [1/s]
bm = @(V) 1000*0.28*(5*double((V-VT-40)==0)+...
    (V-VT-40)./(double((V-VT-40)~=0).*exp(((V-VT-40)/5))-1)); % Rate constant beta_m [1/s]
an = @(V) -1000*0.032*((-5)*double((V-VT-15)==0)+...
    (V-VT-15)./(double((V-VT-15)~=0).*exp(-((V-VT-15)/5))-1));% Rate constant alpha_n [1/s]
bn = @(V) 1000*0.5*exp(-(V-VT-10)/40);                    % Rate constant beta_n [1/s]
ah = @(V) 1000*0.128*exp(-((V-VT-17)/18));                % Rate constant alpha_h [1/s]
bh = @(V) (1000*4)./(1+exp(-((V-VT-40)/5)));               % Rate constant beta_h [1/s]
end
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 6 || MODEL == 7 || MODEL == 8 
pinf = @(V) 1./(1+exp(-((V+35)/10)));                      % Rest p-value [-]
taup = @(V) taumax./(3.3*exp((V+35)/20)+exp(-(V+35)/20));  % Time-constant [s] for p 
end
if MODEL == 3 || MODEL == 4 || MODEL == 8
sinf = @(V) 1./(1+exp(-(V+Vx+57)/6.2));					  % Rest s-value [-]
uinf = @(V) 1./(1+exp((V+Vx+81)/4));				          % Rest u-value [-] 
taus = @(V) 10^(-3)*((1/3.7)*(0.612+1./(exp(-((V+Vx+132)/16.7))+exp((V+Vx+16.8)/18.2))));  % Time-constant [s] for s
tauu = @(V) 10^(-3)*funtauu(V,Vx);                        % Time-constant [s] for u
end
if MODEL == 5
sinf = @(V) 1./(1+exp(-(V+52)/7.4));                       % Rest s-value [-]
uinf = @(V) 1./(1+exp((V+80)/5));                          % rest u-value [-]
taus = @(V) 10^(-3)*(1+0.33./(exp((V+27)/10)+exp(-(V+102)/15))); % Time constant [s] for s
tauu = @(V) 10^(-3)*(28.3+0.33./(exp((V+48)/4)+exp(-(V+407)/50))); % Time constant [s] for u
end
if MODEL == 4
winf = @(V) 1./(1+exp((V+75)/5.5));                        % Rest w-value [-]    
tauw = @(V) 10^(-3)./(exp(-14.59-0.086*V)+exp(-1.87+0.0701*V));  % Time-constant [s] for w
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
%d2inf = @(kcCai) interp1(kcCaiLU',d2infLU,kcCai);

taum = @(V) interp1(VrangeLU',taumLU,V); tauh = @(V) interp1(VrangeLU',tauhLU,V); 
taun = @(V) interp1(VrangeLU',taunLU,V); taup = @(V) interp1(VrangeLU',taupLU,V);
tauq = @(V) interp1(VrangeLU',tauqLU,V); taur = @(kcCai) interp1(kcCaiLU',taurLU,kcCai);
tauc = @(V) interp1(VrangeLU',taucLU,V); taua = @(V) interp1(VrangeLU',tauaLU,V);
taub = @(V) interp1(VrangeLU',taubLU,V); taud1 = @(V) interp1(VrangeLU',taud1LU,V);
%taud2 = @(kcCai) interp1(kcCaiLU',taud2LU,kcCai);
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
%d2inf = @(kcCai) interp1qr(kcCaiLU',d2infLU,kcCai);

taum = @(V) interp1qr(VrangeLU',taumLU,V); tauh = @(V) interp1qr(VrangeLU',tauhLU,V); 
taun = @(V) interp1qr(VrangeLU',taunLU,V); taup = @(V) interp1qr(VrangeLU',taupLU,V);
tauq = @(V) interp1qr(VrangeLU',tauqLU,V); taur = @(kcCai) interp1qr(kcCaiLU',taurLU,kcCai);
tauc = @(V) interp1qr(VrangeLU',taucLU,V); taua = @(V) interp1qr(VrangeLU',tauaLU,V);
taub = @(V) interp1qr(VrangeLU',taubLU,V); taud1 = @(V) interp1qr(VrangeLU',taud1LU,V);
%taud2 = @(kcCai) interp1qr(kcCaiLU',taud2LU,kcCai);
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
%d2inf = @(kcCai) nakeinterp1(kcCaiLU',d2infLU,kcCai);

taum = @(V) nakeinterp1(VrangeLU',taumLU,V); tauh = @(V) nakeinterp1(VrangeLU',tauhLU,V); 
taun = @(V) nakeinterp1(VrangeLU',taunLU,V); taup = @(V) nakeinterp1(VrangeLU',taupLU,V);
tauq = @(V) nakeinterp1(VrangeLU',tauqLU,V); taur = @(kcCai) nakeinterp1(kcCaiLU',taurLU,kcCai);
tauc = @(V) nakeinterp1(VrangeLU',taucLU,V); taua = @(V) nakeinterp1(VrangeLU',tauaLU,V);
taub = @(V) nakeinterp1(VrangeLU',taubLU,V); taud1 = @(V) nakeinterp1(VrangeLU',taud1LU,V);
%taud2 = @(kcCai) nakeinterp1(kcCaiLU',taud2LU,kcCai);
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
USPaT = @(t) USPa;
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
if UpSpeed == 4
[Arx,deltaxx,xx,yx,~,~] = fitLennardJones(a,delta);
Pm = @(Z) Arx*((deltaxx./(2*Z+delta)).^xx-(deltaxx./(2*Z+delta)).^yx); % Lennard-Jones approximation [Pa]
end
% Pm(Z) is not properly defined for Z < -delta/2. Z should be > -delta/2!
if UpSpeed == 1 || UpSpeed == 2 || UpSpeed == 3
if DISPLAY, disp('Initialising molecular force vector...'); end %#ok<*UNRCH>
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
dCmdZ = @(Z) FundCmdZ(Z,Cm0,a,delta); %#ok<NASGU> % derivative dCm/dZ [F/m]
PinLIN = @(Z) Po*(1+(Z./(3*delta)).*(Z.^2/a^2+3)).^(-kappa); %Linear approximation of Pin 
PmdPin = @(Z) Pm(Z)*(1+(Z./(3*delta))*(Z.^2/a^2+3))^(kappa)/Po; % Pm divided by Pinlin
Pin = @(na,Z) (na*Rg*Temp)./(Va(Z));% Internal pressure [Pa]
S = @(Z) pi*(a^2+Z.^2);% Surface of a single leaflet [m^2]
%fVCa = @(cCai) 10^(3)*((Rg*Temp)/(2*Far))*log(cCae./cCai); % Nernst equation for Ca-potential [mV] (if not assumed constant) 

% 3. Initial conditions and timespan
Z0 = 0;
dZ0 = 0;
Pin0 = Po;
na0 = (Pin0*Va(Z0))/(Rg*Temp);  % Mole air in BLS
Ca0 = Ca.*ones(1,Nout);         % External molar air concentration (mole/m^3)
Tmax = 50/USfreq;               % Max runtime before which periodicity has to be established (s) 
if Tmax <= PerCheckMax/USfreq || Tmax <= NonPerEffP/USfreq
    disp('Warning: insufficient number of periods in Tmax to check for periodicity or to calculate effective parameters over NonPerEffP');
end
dtUS = 0.025/USfreq;            % Discretisation time (s)
%dtUS = 0.1/USfreq;
%dtUS = 0.005/USfreq;
dtlin = 0.05*dtUS;

% non charge parameters:
% MODEL 1,2,6,7,13: m,n,p,h
% MODEL 3,8: m,n,p,h,s,u 
% MODEL 4: m, n, h, s, u, w, wLock, hProtein, cCai
% MODEL 5: m, n, h, s, u, cCai
% MODEL 9: m, n, p, h, q, r, a, b, c, d1, d2, cCai
% MODEL 10: h, r
% MODEL 11, 12: n,h,r,CA0
% MODEL 14: m, n, h
X0 = [Z0 dZ0 Pin0 na0 Ca0]';
%Tvalueslin = (max(0,USpstart):dtlin:ceil(min(Tsim,USpstart+USpd)/dtlin)*dtlin);
Tvalueslin = @(it) (it-1)*dtlin;

BLSlin =  @(t,Z,V) PmdPin(Z)+1-(1+(Z./(3*delta))*(Z.^2/a^2+3))^(kappa)+...
    (USPaT(t)/Po)*sin(omega*t)*(1+(Z./(3*delta))*(Z.^2/a^2+3))^(kappa)+...
    (Pec(V,Z)/Po)*(1+(Z./(3*delta))*(Z.^2/a^2+3))^(kappa);  %#ok<NASGU> % Differential-algebraic equation
BLSlinQ =  @(t,Z,Q) PmdPin(Z)+1-(1+(Z./(3*delta))*(Z.^2/a^2+3))^(kappa)+...
    (USPaT(t)/Po)*sin(omega*t)*(1+(Z./(3*delta))*(Z.^2/a^2+3))^(kappa)+...
    (PecQ(Q,Z)/Po)*(1+(Z./(3*delta))*(Z.^2/a^2+3))^(kappa);      % Quasistatic approximation assumes ng = na0 and pressure changes quasistatically. <-> Steady state where P=P0=kh*Ca
BLSlinP0 =  @(Z,Q) Pm(Z)+PecQ(Q,Z)-PS(Z);             % Steady state expression (Pin = P0, Z = 0, USPaT = 0),

if (USPa == 0)              % Solver not required. Also, would take unreasonably long in the quasistatic part. 
    infinit = 1*10^(-12);           % Small number for fzero range
    Zeff = fzero(@(Z) BLSlinP0(Z,Qm),[-delta/2+infinit,Hmax]);      % (m)
    Cmeff = Cm(Zeff);  % (F/m^2)
    Veff = (10^3)*Qm/Cmeff;      % (mV)
    ngend = (Pin0*Va(Zeff))/(Rg*Temp);                    % (mol)
    a_i_eff = struct; apb_i_eff = struct;           % (/s)
	if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 4 || MODEL == 5 || MODEL == 6 || MODEL == 7 || MODEL == 8 || MODEL == 13 || MODEL == 14
    if SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
	apb_i_eff.('m') = ampbm(Veff); apb_i_eff.('n') = anpbn(Veff); apb_i_eff.('h') = ahpbh(Veff);
    else
	apb_i_eff.('m') = am(Veff)+bm(Veff); apb_i_eff.('n') = an(Veff)+bn(Veff); apb_i_eff.('h') = ah(Veff)+bh(Veff);  
    end
	a_i_eff.('m') = am(Veff); a_i_eff.('n') = an(Veff); a_i_eff.('h') = ah(Veff);
	end
	if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 6 || MODEL == 7 || MODEL == 8
	a_i_eff.('p') = pinf(Veff)/taup(Veff); apb_i_eff.('p') = 1/taup(Veff); 
	end
	if MODEL == 4
	a_i_eff.('w') = winf(Veff)/tauw(Veff); apb_i_eff.('w') = 1/tauw(Veff); 
	end
	if MODEL == 3 || MODEL == 4 || MODEL == 5 || MODEL == 8
    a_i_eff.('s') = sinf(Veff)/taus(Veff); a_i_eff.('u') = uinf(Veff)/tauu(Veff);
    apb_i_eff.('s') = 1/taus(Veff); apb_i_eff.('u') = 1/tauu(Veff);
	end
	if MODEL == 9
    a_i_eff.('m') = minf(Veff)/taum(Veff); a_i_eff.('n') = ninf(Veff)/taun(Veff);  a_i_eff.('h') = hinf(Veff)/tauh(Veff);  
	a_i_eff.('p') = pinf(Veff)/taup(Veff); a_i_eff.('q') = qinf(Veff)/tauq(Veff);  
	a_i_eff.('a') = ainf(Veff)/taua(Veff);  a_i_eff.('b') = binf(Veff)/taub(Veff); a_i_eff.('c') = cinf(Veff)/tauc(Veff); a_i_eff.('d1') = d1inf(Veff)/taud1(Veff);
    apb_i_eff.('m') = 1/taum(Veff); apb_i_eff.('n') = 1/taun(Veff); apb_i_eff.('h') = 1/tauh(Veff); apb_i_eff.('p') = 1/taup(Veff); 
	apb_i_eff.('q') = 1/tauq(Veff);	apb_i_eff.('a') = 1/taua(Veff); apb_i_eff.('b') = 1/taub(Veff); apb_i_eff.('c') = 1/tauc(Veff); apb_i_eff.('d1') = 1/taud1(Veff);
	end
	if MODEL == 10 || MODEL == 11 || MODEL == 12
	a_i_eff.('h') = hinf(Veff)/tauh(Veff); a_i_eff.('r') = rinf(Veff)/taur(Veff);
	apb_i_eff.('h') = 1/tauh(Veff); apb_i_eff.('r') = 1/taur(Veff);
	end
	if MODEL == 11 || MODEL == 12
	a_i_eff.('n') = ninf(Veff)/taun(Veff); apb_i_eff.('n') = 1/taun(Veff);
	end
	if MODEL == 13
    if SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
	apb_i_eff.('p') = appbp(Veff);
    else
	apb_i_eff.('p') = ap(Veff)+bp(Veff);
    end
	a_i_eff.('p') = ap(Veff); 
	end
    return;
end

% 4. Solver
global reverseStr; 
global temp; global tempi; global maxlags;
if DISPLAY, disp('Solver started'); end
if SpeedUp == 1 || SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
    temp = cell(ceil(2/USfreq/dtUS),3);               % Estimated preallocation for speed
    tempi = 0;
end
it=2; itB=it; X = X0';
Event = abs(Z0)<abs(Hmin); 
TvaluesX = Tvalueslin(1); 

    % A. Calculate Z,dZ,na,Pin
    % A.1. Linear initialization
tspan = [0,Tmax];
infinit = 1*10^(-12);           % Small number for fzero range
Zold=Z0;
while Tvalueslin(it)<=tspan(2) && Event
Z = fzero(@(Z) BLSlinQ(Tvalueslin(it),Z,Qm),[-delta/2+infinit,Hmax]);
if Z <= -delta/2
    error('Unphysical solution: Z<=-delta/2');
end
dZ = (Z-Zold)/dtlin;       % Backward Euler (O(dt))
X(it,1)=Z; X(it,2)=dZ;
it=it+1; Zold=Z; Event=abs(Z)<abs(Hmin);
end
itE = it-1; 
TvaluesX = horzcat(TvaluesX,Tvalueslin(itB:itE)); %#ok<NASGU>
X(itB:itE,3) = PinLIN(X(itB:itE,1));
X(itB:itE,4) = (X(itB:itE,3).*Va(X(itB:itE,1)))./(Rg*Temp);
% Solve diffusion equation 
if ~DiffusionApproximation
diffLin = @(t,Ca) A*Ca+B(interp1(Tvalueslin(itB-1:tE)',X(itB-1:itE,3)',t));
[~, Ca] = ode113(diffLin,Tvalueslin(itB-1:itE),Ca0');
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

if (Tvalueslin(it)>tspan(2))
    error('Quasistatic solver not finished before Tmax. - Consider increasing Tmax');
end
reverseStr = '';    % String to tracks progress.
% A.2 Non-linear BLS equation
if DISPLAY==2, disp('Solve Rayleigh-Plesset equation'); end
tBLS=[Tvalueslin(itE) tspan(2)];
W0 = [X(end,1:2), X(end,4:4+Nout)]'; % W = [Z dZ na Ca]
if DiffusionApproximation
%BLS = @(t,Z,dZ,na) BLS1Q(DISPLAY,tBLS,t,Z,dZ,na,R,rhol,PecQ,Qm,Pin,Pm,Po,PaT,omega,PS,delta0,mus,mul,S,Da,Ca,ka,ksi); %#ok<*UNRCH>
if SpeedUp == 1 || SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
EventFcn = @(t,W) EventFcn1(t,W(1),W(2),W(3),USfreq,dtUS,CORRTres,PerCheckMax);
OdeOpts = odeset('MaxStep',dtUS,'Events',EventFcn,'AbsTol',[1e-13;1e-9;1e-24],'RelTol',1e-4);
%OdeOpts = odeset('MaxStep',dtUS,'Events',EventFcn,'AbsTol',1e-3,'RelTol',1e-1);
end
% For large displacements (i.e. high intensity, low frequency...) a stiff
% solver is significantly faster (e.g. ode23t) (Large Pm gradient)
% For small displacements, a non-stiff solver (e.g. ode113) is faster. Note
% that in Lemaire et al. a hybrid stiff-nonstiff python solver is used.
% Because we have to choose, the stiff solver (ode23t) is better, because
% simulations for large displacements are the bottlenecks (are more time
% intensive than small displacements)
ME = ''; MEb = '';          
try [t,W] = ode23t(@(t,W) BLS1Q(DISPLAY,tBLS,t,W(1),W(2),W(3),R,rhol,PecQ,Qm,Pin,Pm,Po,USPaT,omega,PS,delta0,mus,mul,S,Da,Ca,ka,ksi),tBLS,W0(1:3),OdeOpts);
catch ME
end
if ~isempty(ME)
if strcmp(ME.identifier,'FunPm:UNPHYS')     % Unphysical solution might be caused by instability in the ode23t solver
try [t,W] = ode113(@(t,W) BLS1Q(DISPLAY,tBLS,t,W(1),W(2),W(3),R,rhol,PecQ,Qm,Pin,Pm,Po,USPaT,omega,PS,delta0,mus,mul,S,Da,Ca,ka,ksi),tBLS,W0(1:3),OdeOpts);
catch MEb   % We retry with ode113
end
if ~isempty(MEb)
error(MEb.message);
end
else
error(ME.message);
end
end
SONICPer = maxlags;
else
%BLS = @(t,Z,dZ,na,Ca) BLS2Q(DISPLAY,tBLS,t,Z,dZ,na,Ca,R,rhol,PecQ,Qm,Pin,Pm,Po,PaT,omega,PS,delta0,mus,mul,S,Da,ka,A,B,deltaR);
if SpeedUp == 1 || SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6 
EventFcn = @(t,W) EventFcn2(t,W(1),W(2),W(3),W(4:3+Nout),USfreq,dtUS,CORRTres,PerCheckMax);
OdeOpts = odeset('MaxStep',dtUS,'Events',EventFcn);
end
[t,W] = ode113(@(t,W) BLS2Q(DISPLAY,tBLS,t,W(1),W(2),W(3),W(4:3+Nout),R,rhol,PecQ,Qm,Pin,Pm,Po,USPaT,omega,PS,delta0,mus,mul,S,Da,ka,A,B,deltaR),tBLS,W0,OdeOpts);
SONICPer = maxlags;
end
if t(end) == tBLS(end)
disp(['Note: no periodicity is found with (CORRtres,Tmax,Qm,USPa,USfreq,aBLS): ' num2str(CORRTres) ',' num2str(Tmax) ',' num2str(Qm) ',' num2str(USPa) ',' num2str(USfreq) ',' num2str(aBLS)]);
SONICPer = NonPerEffP/USfreq;
end
if DISPLAY==2
disp(' ');
disp(['Periodicity found in solution with autocorrelation >' num2str(CORRTres)]);
disp('Ending Rayleigh-Plesset solver...');
end
% We now extract the periodic part 
% 1. W is (approximately with CORRTres) periodic with T=1/freq
% -> It is important to replicate EXACTLY the period, as small errors will
%    accumulate quickly
% -> Very important: the endconditions on X (or W) have to be conserved
%    over an integer times the period (otherwise discontinuities will
%    crash the program!)
% 2. length Dt = t(end)-t(firstIn) (approx<)= 1/freq

firstIn = find((t(end)-t)<SONICPer,1);
tPeriod = [(t(end)-SONICPer), t(firstIn:end)'];
Wperiod = vertcat(W(end,:),W(firstIn:end,:));

% Calculating effective parameters 
% Note: trapz (trapezoidal integration) is the same as the integral of the linear interpolation. 
% So actually trapz is more accurate than integral(nakeinterp1(...)), because in the latter twice higher order errors
if DISPLAY == 2, disp('Calculating effective parameters'); end
calcEff = @(f) (1/(tPeriod(end)-tPeriod(1)))*trapz(tPeriod,f(Wperiod(:,1)));
calcEffrate = @(f) (1/(tPeriod(end)-tPeriod(1)))*trapz(tPeriod,f(1000*Qm./Cm(Wperiod(:,1))));
Zeff = calcEff(@(X) X);             % (m)
Veff = 1000*calcEff(@(X) Qm./Cm(X));   % (mV)
Cmeff = calcEff(Cm);        % (F/m^2)
ngend = Wperiod(end,3);

% ------------------------------Rate constants-----------------------------
a_i_eff = struct; apb_i_eff = struct;
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 4 || MODEL == 5 || MODEL == 6 || MODEL == 7 || MODEL == 8 || MODEL == 13 || MODEL == 14
	ameff = calcEffrate(am);
	aneff = calcEffrate(an);
	aheff = calcEffrate(ah);
    if SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
	ampbmeff = calcEffrate(ampbm);
	anpbneff = calcEffrate(anpbn);
	ahpbheff = calcEffrate(ahpbh);
    else
	bmeff = calcEffrate(bm);
	bneff = calcEffrate(bn);
	bheff = calcEffrate(bh);
	
	ampbmeff = ameff+bmeff;
	anpbneff = aneff+bneff;
	ahpbheff = aheff+bheff;    
    end
	a_i_eff.('m') = ameff; a_i_eff.('n') = aneff; a_i_eff.('h') = aheff;
	apb_i_eff.('m') = ampbmeff; apb_i_eff.('n') = anpbneff; apb_i_eff.('h') = ahpbheff;
end
if MODEL == 1 || MODEL == 2 || MODEL == 3 || MODEL == 6 || MODEL == 7 || MODEL == 8
    apeff = calcEffrate(@(X) pinf(X)./taup(X));
	appbpeff = calcEffrate(@(X) 1./taup(X));
	a_i_eff.('p') = apeff; apb_i_eff.('p') = appbpeff; 
end
if MODEL == 4
	aweff = calcEffrate(@(X) winf(X)./tauw(X));
	awpbweff = calcEffrate(@(X) 1./tauw(X));
	a_i_eff.('w') = aweff; apb_i_eff.('w') = awpbweff; 
end
if MODEL == 3 || MODEL == 4 || MODEL == 5 || MODEL == 8
    aseff = calcEffrate(@(X) sinf(X)./taus(X));
	aspbseff = calcEffrate(@(X) 1./taus(X));
    aueff = calcEffrate(@(X) uinf(X)./tauu(X));
    aupbueff = calcEffrate(@(X) 1./tauu(X));
    
    a_i_eff.('s') = aseff; a_i_eff.('u') = aueff;
    apb_i_eff.('s') = aspbseff; apb_i_eff.('u') = aupbueff;
end
if MODEL == 9
	ameff = calcEffrate(@(X) minf(X)./taum(X));
	ampbmeff = calcEffrate(@(X) 1./taum(X));
	aneff = calcEffrate(@(X) ninf(X)./taun(X));
	anpbneff = calcEffrate(@(X) 1./taun(X));
	aheff = calcEffrate(@(X) hinf(X)./tauh(X));
	ahpbheff = calcEffrate(@(X) 1./tauh(X));
	apeff = calcEffrate(@(X) pinf(X)./taup(X));
	appbpeff = calcEffrate(@(X) 1./taup(X));
	aqeff = calcEffrate(@(X) qinf(X)./tauq(X));
	aqpbqeff = calcEffrate(@(X) 1./tauq(X));
	aaeff = calcEffrate(@(X) ainf(X)./taua(X));
	aapbaeff = calcEffrate(@(X) 1./taua(X));
	abeff = calcEffrate(@(X) binf(X)./taub(X));
	abpbbeff = calcEffrate(@(X) 1./taub(X));
	aceff = calcEffrate(@(X) cinf(X)./tauc(X));
	acpbceff = calcEffrate(@(X) 1./tauc(X));
	ad1eff = calcEffrate(@(X) d1inf(X)./taud1(X));
	ad1pbd1eff = calcEffrate(@(X) 1./taud1(X));

    
    a_i_eff.('m') = ameff; a_i_eff.('n') = aneff;  a_i_eff.('h') = aheff;  a_i_eff.('p') = apeff;  a_i_eff.('q') = aqeff;  
	a_i_eff.('a') = aaeff;  a_i_eff.('b') = abeff;  a_i_eff.('c') = aceff;  a_i_eff.('d1') = ad1eff;
    apb_i_eff.('m') = ampbmeff; apb_i_eff.('n') = anpbneff; apb_i_eff.('h') = ahpbheff; apb_i_eff.('p') = appbpeff; apb_i_eff.('q') = aqpbqeff;
	apb_i_eff.('a') = aapbaeff; apb_i_eff.('b') = abpbbeff; apb_i_eff.('c') = acpbceff; apb_i_eff.('d1') = ad1pbd1eff;
end
if MODEL == 10 || MODEL == 11 || MODEL == 12
	aheff = calcEffrate(@(X) hinf(X)./tauh(X));
	ahpbheff = calcEffrate(@(X) 1./tauh(X));
	areff = calcEffrate(@(X) rinf(X)./taur(X));
	arpbreff = calcEffrate(@(X) 1./taur(X));
	
	a_i_eff.('h') = aheff; a_i_eff.('r') = areff;
	apb_i_eff.('h') = ahpbheff; apb_i_eff.('r') = arpbreff;
end
if MODEL == 11 || MODEL == 12
	aneff = calcEffrate(@(X) ninf(X)./taun(X));
	anpbneff = calcEffrate(@(X) 1./taun(X));
	
	a_i_eff.('n') = aneff; apb_i_eff.('n') = anpbneff;
end
if MODEL == 13
	apeff = calcEffrate(ap);
    if SpeedUp == 2 || SpeedUp == 3 || SpeedUp == 4 || SpeedUp == 5 || SpeedUp == 6
	appbpeff = calcEffrate(appbp);
    else
	bpeff = calcEffrate(bp);
	appbpeff = apeff+bpeff;  
    end
	a_i_eff.('p') = apeff; 
	apb_i_eff.('p') = appbpeff; 
end

TTime = toc; %#ok<NASGU>
if DISPLAY
disp(' ');
disp(['Program finished in ' num2str(round(TTime,1)) 's']);
disp('Post-processing...');
end

end
