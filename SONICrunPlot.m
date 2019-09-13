% Plot results in order to reproduce results in Lemaire et al.
c = 1515;				% Speed of sound surrounding medium (m/s)
rhol = 1028;			% Density surrounding medium (kg/m^3)
Pa2I = @(Pa) Pa.^2/(2*rhol*c);
I2Pa = @(I) sqrt(2*rhol*c*I);
NICEpath = 'D:\users\ttarnaud\8. Piezoelectric solver\Parallellized functions for HPC calculations';
SONICpath = 'D:\users\ttarnaud\8. Piezoelectric solver\8.4. Lemaire et al. (2018) - SONIC solver';
debugSwitch = nan;            % Number: only run part of the program. nan : run everything 
FigurePlot = 7;             % Number: plot this figure number. nan: plot everything

%% FIGURE 5 Lemaire et al. (2018)
if isnan(FigurePlot) || FigurePlot == 5
MODEL = 1;
switch MODEL
    case 1, MODELstr = 'RS';
    case 2, MODELstr = 'FS';
    case 3, MODELstr = 'LTS';
end

figure; set(gcf,'color','w'); set(gca,'box','off');
% Fig.1(a-c) (top)
Tsim = 0.3;     % (s)
USpd = 0.150;   % (s)
USps = 0.05;  % (s)
USdc = 1; USprf = 0;   % (-) , (Hz)

USfreq = {(500e3),[20e3,4e6],(500e3)}; % (Hz)
aBLS = {(32e-9),(32e-9),[16e-9,64e-9]};  % (m)
RelPa = {[-5e3,0,20e3],(20e3),(20e3)};  % relative to threshold (Pa)

SONICcolors = {'b','g','r','y','m',[0.4660 0.6740 0.1880],[0.6350 0.0780 0.1840]};
iNICEcolor = 0; iSONICcolor = 0;

legendStr = {{'A_T - 5 kPa','A_T','A_T + 20 kPa'},{'20 kHz','4 MHz'},{'16 nm','64 nm'}};

cd(SONICpath); 
fprintf('Continuous wave stimulation: charge trace plots \n');
if isnan(debugSwitch) || debugSwitch == 1
for i = 1:3
plotnr = 0;
fprintf('Calculation of subplot (%u/3) \n',i);
subplot(7,3,i);
hold on;
for iFreq = 1:length(USfreq{i})
for iaBLS = 1:length(aBLS{i})
SONICrun(num2str(Tsim),'1',num2str(USps),num2str(USpd),num2str(USfreq{i}(iFreq)),num2str(USdc),num2str(USprf),...
'0','0','0','1','0','0','0',num2str(MODEL),'0','0','0',num2str(aBLS{i}(iaBLS)),'1');
lt = load(['Thresh(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq{i}(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(0) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS{i}(iaBLS)) ').mat']);
PaThreshSONIC = I2Pa(lt.IIpa);
PaThRange = PaThreshSONIC+RelPa{i};
delete(['Thresh(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq{i}(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(0) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS{i}(iaBLS)) ').mat']);

% NICE around threshold
cd(NICEpath);
for iPa = 1:length(PaThRange)
funPESa(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq{i}(iFreq)),num2str(USdc),num2str(USprf),...
num2str(Pa2I(PaThRange(iPa))),'0','0','1','0','0','0',num2str(MODEL),'0','0','0',num2str(aBLS{i}(iaBLS)));
end
for iPa = 1:length(PaThRange)
plotnr = plotnr+1;
lt2 = load(['Chargevt(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq{i}(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaThRange(iPa))) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS{i}(iaBLS)) ').mat']);
Qvt = lt2.saveChargeSample(:,2); timeline = lt2.saveChargeSample(:,1);
iNICEcolor = iNICEcolor+1;
g{i}(plotnr) = plot(10^3*timeline,Qvt,'linestyle','-','color',SONICcolors{iNICEcolor},'linewidth',2); %#ok<SAGROW>
delete(['Chargevt(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq{i}(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaThRange(iPa))) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS{i}(iaBLS)) ').mat']);
delete(['APtimes(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq{i}(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaThRange(iPa))) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS{i}(iaBLS)) ').mat']);
end
% SONIC around threshold
cd(SONICpath);
for iPa = 1:length(PaThRange)
SONICrun(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq{i}(iFreq)),num2str(USdc),num2str(USprf),...
num2str(Pa2I(PaThRange(iPa))),'0','0','1','0','0','0',num2str(MODEL),'0','0','0',num2str(aBLS{i}(iaBLS)),'1');
end
for iPa = 1:length(PaThRange)
lt3 = load(['Chargevt(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq{i}(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaThRange(iPa))) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS{i}(iaBLS)) ').mat']);
Qvt = lt3.saveChargeSample(:,2); timeline = lt3.saveChargeSample(:,1); 
iSONICcolor = iSONICcolor+1;
plot(10^3*timeline,Qvt,'linestyle','--','color',SONICcolors{iSONICcolor});
delete(['Chargevt(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq{i}(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaThRange(iPa))) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS{i}(iaBLS)) ').mat']);
delete(['APtimes(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq{i}(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaThRange(iPa))) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS{i}(iaBLS)) ').mat']);
end
end
end
if (i==1)
    ylabel('Q [nC/cm^2]');
end
xlabel('Time [ms]');
legend(g{i},legendStr{i});
end
end

% Fig. 1(a-c) bottom
USfreq = {(500e3), logspace(log10(20e3),log10(4e6),7),(500e3)}; % (Hz)
aBLS = {(32e-9),(32e-9),logspace(log10(16e-9),log10(64e-9),5)};  % (m)
RelPa = {'',(20e3),(20e3)};  % relative to threshold (Pa)
PaR = {1e3*logspace(log10(50),log10(150),4),'',''};   % Pressure range (Pa) 

fprintf('Continuous wave stimulation: latency, firing rate and spike amplitude plots \n');
if isnan(debugSwitch) || debugSwitch == 2
for i = 1:3
fprintf('Calculation of subplot (%u/3) \n',i);

USfreqRange = USfreq{i}; aBLSRange = aBLS{i};
RelPaRange = RelPa{i}; PaRRange = PaR{i};

ThreshPa = zeros(length(USfreqRange),length(aBLSRange));

if ~isempty(RelPaRange)          % Threshold is required
cd(SONICpath);
for iFreq = 1:length(USfreqRange)
    for iaBLS = 1:length(aBLSRange)
    SONICrun(num2str(Tsim),'1',num2str(USps),num2str(USpd),num2str(USfreqRange(iFreq)),num2str(USdc),num2str(USprf),...
    '0','0','0','1','0','0','0',num2str(MODEL),'0','0','0',num2str(aBLSRange(iaBLS)),'1');
    ltf = load(['Thresh(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreqRange(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(0) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLSRange(iaBLS)) ').mat']);
    ThreshPa(iFreq,iaBLS) = I2Pa(ltf.IIpa);
    delete(['Thresh(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreqRange(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(0) ')-ES(0,0,1,0,0)-aBLS=(' ... 
    num2str(aBLSRange(iaBLS)) ').mat']);
    end
end
else
RelPaRange = PaRRange;      % Now relative to TreshPa = 0
end

Latency = zeros(2,length(USfreqRange),length(aBLSRange),length(RelPaRange));            % Latency (ms)
FR =  zeros(2,length(USfreqRange),length(aBLSRange),length(RelPaRange));                % Firing rate (Hz)
spA =  zeros(2,length(USfreqRange),length(aBLSRange),length(RelPaRange));               % spike amplitude (nC/cm^2)
for pathNum = 1:2
if pathNum == 1
    cd(SONICpath);
elseif pathNum == 2
    cd(NICEpath);
end
    for iFreq = 1:length(USfreqRange)
        for iaBLS = 1:length(aBLSRange)
            for iPaR = 1:length(RelPaRange)
            appPa = ThreshPa(iFreq,iaBLS)+RelPaRange(iPaR);
            if pathNum == 1
            SONICrun(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreqRange(iFreq)),num2str(USdc),num2str(USprf),...
            num2str(Pa2I(appPa)),'0','0','1','0','0','0',num2str(MODEL),'0','0','0',num2str(aBLSRange(iaBLS)),'1');
            elseif pathNum == 2
            funPESa(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreqRange(iFreq)),num2str(USdc),num2str(USprf),...
            num2str(Pa2I(appPa)),'0','0','1','0','0','0',num2str(MODEL),'0','0','0',num2str(aBLSRange(iaBLS)));
            end
            end
        end
    end
    
    for iFreq = 1:length(USfreqRange)
        for iaBLS = 1:length(aBLSRange)
            for iPaR = 1:length(RelPaRange)
            appPa = ThreshPa(iFreq,iaBLS)+RelPaRange(iPaR);

            ll = load(['APtimes(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
            num2str(USfreqRange(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(appPa)) ')-ES(0,0,1,0,0)-aBLS=('...
            num2str(aBLSRange(iaBLS)) ').mat']);
            lr = load(['Chargevt(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
            num2str(USfreqRange(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(appPa)) ')-ES(0,0,1,0,0)-aBLS=(' ...
            num2str(aBLSRange(iaBLS)) ').mat']);
            APtimes = ll.APtimes(ll.APtimes>=USps&ll.APtimes<=(USps+USpd)); %(s)
            
            if isempty(APtimes)
            Latency(pathNum,iFreq,iaBLS,iPaR) = inf;
            else
            Latency(pathNum,iFreq,iaBLS,iPaR) = 10^3*(APtimes(1)-USps);    % (ms)
            end
            FR(pathNum,iFreq,iaBLS,iPaR) = length(APtimes)./(USpd);        % (Hz)
            
            if length(APtimes)>=4
            Qvt = lr.saveChargeSample(:,2); timeline = lr.saveChargeSample(:,1);  % (nC/cm^2) and (s)
            QvtCrop = Qvt(timeline>=APtimes(2)&timeline<=APtimes(end-1));  % Crop outer lying action potentials to determine amplitude
            timelineCrop = timeline(timeline>=APtimes(2)&timeline<=APtimes(end-1));
            spA(pathNum,iFreq,iaBLS,iPaR) = (max(QvtCrop)-min(QvtCrop));       % (nC/cm^2)
            else
            spA(pathNum,iFreq,iaBLS,iPaR) = nan;
            end
            
            delete(['APtimes(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
            num2str(USfreqRange(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(appPa)) ')-ES(0,0,1,0,0)-aBLS=('...
            num2str(aBLSRange(iaBLS)) ').mat']);
            delete(['Chargevt(' MODELstr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
            num2str(USfreqRange(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(appPa)) ')-ES(0,0,1,0,0)-aBLS=(' ...
            num2str(aBLSRange(iaBLS)) ').mat']);            
            
            end
        end
    end
end

% Plotting results
LatencySONIC = Latency(1,:); LatencyNICE = Latency(2,:);
FRSONIC = FR(1,:); FRNICE = FR(2,:);
spASONIC = spA(1,:); spANICE = spA(2,:);

switch i
    case 1, xaxis = 1e-3*(ThreshPa(1,1)+RelPaRange(:)); xlab = 'Amplitude (kPa)';  
    case 2, xaxis = 1e-3*USfreq{2}; xlab = 'Frequency (kHz)';
    case 3, xaxis = 1e9*aBLS{3}; xlab = 'Sonophore radius (nm)';

end
subplot(7,3,3+i);
hold on;
plot(xaxis,LatencyNICE,'linestyle','-','color',[0.5 0.5 0.5],'linewidth',2);
plot(xaxis,LatencySONIC,'linestyle','--','color','k','marker','o');
xlim([xaxis(1),xaxis(end)]);
ylabel('Latency [ms]');
hold off;
subplot(7,3,6+i);
hold on;
plot(xaxis,FRNICE,'linestyle','-','color',[0.5 0.5 0.5],'linewidth',2);
plot(xaxis,FRSONIC,'linestyle','--','color','k','marker','o');
xlim([xaxis(1),xaxis(end)]);
ylabel('Firing rate [Hz]');
hold off;
subplot(7,3,9+i);
hold on;
plot(xaxis,spANICE,'linestyle','-','color',[0.5 0.5 0.5],'linewidth',2);
plot(xaxis,spASONIC,'linestyle','--','color','k','marker','o');
xlim([xaxis(1),xaxis(end)]);
ylabel('Spike amp. [nC/cm^2]');
xlabel(xlab);
hold off;
end
end

% Fig. 1(d-e) top
aBLS =  32e-9; % (m)
USdc = 0.05; % (-)
USfreq = 500e3; % (Hz)
PaR = 100e3;  % relative to threshold (Pa)

USprf = {(100),[10,100,1000,10000]};   % (Hz)
MODELnr = {[1,3],(3)};
MODELstr = {'RS','FS','LTS'};

iNICEcolor = 0; iSONICcolor = 0;
SONICcolors = {'b','r','r',[0.9290 0.6940 0.1250],'m',[0.6350 0.0780 0.1840]};

legendStr = {{'RS neuron','FS neuron'},{'10 Hz','100 Hz','1 kHz','10 kHz'}};

fprintf('Pulsed wave stimulation: charge trace plots \n');
if isnan(debugSwitch) || debugSwitch == 3
for i = 1:2
plotnr = 0;
fprintf('Calculation of subplot (%u/2) \n',i);
subplot(7,3,12+i);
hold on;
for iprf = 1:length(USprf{i})
for imodel = 1:length(MODELnr{i})
plotnr = plotnr+1;
% NICE around threshold
cd(NICEpath);
funPESa(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf{i}(iprf)),...
num2str(Pa2I(PaR)),'0','0','1','0','0','0',num2str(MODELnr{i}(imodel)),'0','0','0',num2str(aBLS));

lt2 = load(['Chargevt(' MODELstr{MODELnr{i}(imodel)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf{i}(iprf)) ',' num2str(Pa2I(PaR)) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS) ').mat']);
Qvt = lt2.saveChargeSample(:,2); timeline = lt2.saveChargeSample(:,1);
iNICEcolor = iNICEcolor+1;
h(plotnr) = plot(10^3*timeline,Qvt,'linestyle','-','color',SONICcolors{iNICEcolor},'linewidth',2); %#ok<SAGROW>
delete(['Chargevt(' MODELstr{MODELnr{i}(imodel)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf{i}(iprf)) ',' num2str(Pa2I(PaR)) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS) ').mat']);
delete(['APtimes(' MODELstr{MODELnr{i}(imodel)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf{i}(iprf)) ',' num2str(Pa2I(PaR)) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS) ').mat']);
% SONIC around threshold
cd(SONICpath);

SONICrun(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf{i}(iprf)),...
num2str(Pa2I(PaR)),'0','0','1','0','0','0',num2str(MODELnr{i}(imodel)),'0','0','0',num2str(aBLS),'1');

lt3 = load(['Chargevt(' MODELstr{MODELnr{i}(imodel)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf{i}(iprf)) ',' num2str(Pa2I(PaR)) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS) ').mat']);
Qvt = lt3.saveChargeSample(:,2); timeline = lt3.saveChargeSample(:,1); 
iSONICcolor = iSONICcolor+1;
plot(10^3*timeline,Qvt,'linestyle','--','color',SONICcolors{iSONICcolor});
delete(['Chargevt(' MODELstr{MODELnr{i}(imodel)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf{i}(iprf)) ',' num2str(Pa2I(PaR)) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS) ').mat']);
delete(['APtimes(' MODELstr{MODELnr{i}(imodel)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf{i}(iprf)) ',' num2str(Pa2I(PaR)) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS) ').mat']);
end
end
legend(h,legendStr{i});
end
if (i==1)
    ylabel('Q [nC/cm^2]');
end
xlabel('Time [ms]');
end

%Fig. 1(d-e) bottom
USfreq = 500e3; % (Hz) 
PaR = 100e3;  % relative to threshold (Pa) 
iNICEcolor = 0; iSONICcolor = 0;

USdc = {linspace(0.05,1,5),0.05}; % (-) 
USprf = {(100),logspace(log10(10),log10(10e3),10)};   % (Hz)
MODELnr = {[1,3],(3)};

NICEcolor =  {'b','r',[0.5 0.5 0.5]};
SONICcolor = {'b','r','k'};

fprintf('Pulsed wave stimulation: latency, firing rate and spike amplitude plots \n');
if isnan(debugSwitch) || debugSwitch == 4
for i = 1:2
fprintf('Calculation of subplot (%u/2) \n',i);

USdcRange = USdc{i}; USprfRange = USprf{i};
MODELnrRange = MODELnr{i};

Latency = zeros(2,length(USdcRange),length(USprfRange),length(MODELnrRange));            % Latency (ms)
FR =  zeros(2,length(USdcRange),length(USprfRange),length(MODELnrRange));                % Firing rate (Hz)
spA =  zeros(2,length(USdcRange),length(USprfRange),length(MODELnrRange));               % spike amplitude (nC/cm^2)
for pathNum = 1:2
if pathNum == 1
    cd(SONICpath);
elseif pathNum == 2
    cd(NICEpath);
end
    for idc = 1:length(USdcRange)
        for iprf = 1:length(USprfRange)
            for iMODEL = 1:length(MODELnrRange)
            if pathNum == 1
            SONICrun(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdcRange(idc)),num2str(USprfRange(iprf)),...
            num2str(Pa2I(PaR)),'0','0','1','0','0','0',num2str(MODELnrRange(iMODEL)),'0','0','0',num2str(aBLS),'1');
            elseif pathNum == 2
            funPESa(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdcRange(idc)),num2str(USprfRange(iprf)),...
            num2str(Pa2I(PaR)),'0','0','1','0','0','0',num2str(MODELnrRange(iMODEL)),'0','0','0',num2str(aBLS));
            end
            end
        end
    end
    
    for idc = 1:length(USdcRange)
        for iprf = 1:length(USprfRange)
            for iMODEL = 1:length(MODELnrRange)
            ll = load(['APtimes(' MODELstr{MODELnrRange(iMODEL)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
            num2str(USfreq) ',' num2str(USdcRange(idc)) ',' num2str(USprfRange(iprf)) ',' num2str(Pa2I(PaR)) ')-ES(0,0,1,0,0)-aBLS=('...
            num2str(aBLS) ').mat']);
            lr = load(['Chargevt(' MODELstr{MODELnrRange(iMODEL)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
            num2str(USfreq) ',' num2str(USdcRange(idc)) ',' num2str(USprfRange(iprf)) ',' num2str(Pa2I(PaR)) ')-ES(0,0,1,0,0)-aBLS=(' ...
            num2str(aBLS) ').mat']);
            APtimes = ll.APtimes(ll.APtimes>=USps&ll.APtimes<=(USps+USpd)); %(s)
            
            if isempty(APtimes)
            Latency(pathNum,idc,iprf,iMODEL) = inf;
            else
            Latency(pathNum,idc,iprf,iMODEL) = 10^3*(APtimes(1)-USps);    % (ms)
            end
            FR(pathNum,idc,iprf,iMODEL) = length(APtimes)./(USpd);        % (Hz)
            
            if length(APtimes) >= 4
            Qvt = lr.saveChargeSample(:,2); timeline = lr.saveChargeSample(:,1);  % (nC/cm^2) and (s)
            QvtCrop = Qvt(timeline>=APtimes(2)&timeline<=APtimes(end-1));  % Crop outer lying action potentials to determine amplitude
            timelineCrop = timeline(timeline>=APtimes(2)&timeline<=APtimes(end-1));
            spA(pathNum,idc,iprf,iMODEL) = (max(QvtCrop)-min(QvtCrop));       % (nC/cm^2)
            else
            spA(pathNum,idc,iprf,iMODEL) = nan;
            end
            
            delete(['APtimes(' MODELstr{MODELnrRange(iMODEL)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
            num2str(USfreq) ',' num2str(USdcRange(idc)) ',' num2str(USprfRange(iprf)) ',' num2str(Pa2I(PaR)) ')-ES(0,0,1,0,0)-aBLS=('...
            num2str(aBLS) ').mat']);
            delete(['Chargevt(' MODELstr{MODELnrRange(iMODEL)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
            num2str(USfreq) ',' num2str(USdcRange(idc)) ',' num2str(USprfRange(iprf)) ',' num2str(Pa2I(PaR)) ')-ES(0,0,1,0,0)-aBLS=(' ...
            num2str(aBLS) ').mat']);            
            
            end
        end
    end
end
 
% Plotting results
for iMODEL = 1:length(MODELnrRange)
iSONICcolor = iSONICcolor+1; iNICEcolor = iNICEcolor+1;

LatencySONIC = Latency(1,:,:,iMODEL); LatencyNICE = Latency(2,:,:,iMODEL);
FRSONIC = FR(1,:,:,iMODEL); FRNICE = FR(2,:,:,iMODEL);
spASONIC = spA(1,:,:,iMODEL); spANICE = spA(2,:,:,iMODEL);

switch i
    case 1, xaxis = 100*USdc{i}; xlab = 'Duty cycle (Hz)';  
    case 2, xaxis =  USprf{i}; xlab = 'Pulse repetition frequency (Hz)';
end

subplot(7,3,15+i);
hold on;
plot(xaxis,LatencyNICE(:),'linestyle','-','color',NICEcolor(iNICEcolor),'linewidth',2);
plot(xaxis,LatencySONIC(:),'linestyle','--','color',SONICcolor(iSONICcolor),'marker','o');
xlim([xaxis(1),xaxis(end)]);
ylabel('Latency [ms]');
hold off;
subplot(7,3,18+i);
hold on;
plot(xaxis,FRNICE(:),'linestyle','-','color',NICEcolor(iNICEcolor),'linewidth',2);
plot(xaxis,FRSONIC(:),'linestyle','--','color',SONICcolor(iSONICcolor),'marker','o');
xlim([xaxis(1),xaxis(end)]);
ylabel('Firing rate [Hz]');
hold off;
% subplot(7,3,21+i);
% hold on;
% plot(xaxis,spANICE(:),'linestyle','-','color',NICEcolor(iNICEcolor),'linewidth',2);
% plot(xaxis,spASONIC(:),'linestyle','--','color',SONICcolor(iSONICcolor),'marker','o');
% xlim([xaxis(1),xaxis(end)]);
% ylabel('Spike amp. [nC/cm^2]');
% xlabel(xlab);
% hold off;
end

end
end
set(findobj('type','axes'),'fontsize',18);
set(findobj('type','axes'),'box','off');
end

%% FIGURE 7 Lemaire et al. (2018)
% (a) Color maps 
if isnan(FigurePlot) || FigurePlot == 7
MODELstr = {'RS','FS','LTS'};

Tsim = 2;     % (s)
USpd = 1;   % (s)
USps = 0.5;  % (s)
USfreq = 500e3; % (Hz)
aBLS = 32e-9;  % (m)

USdc = linspace(0.01,1,100); % (-)
USprf = [10,100,1000];          % Hz   
PaR = logspace(log10(10e3),log10(600e3),30); % (Pa)
Modelnr = [1,3];
[mModelnr,mUSprf,mUSdc,mPaR] = ndgrid(Modelnr,USprf,USdc,PaR);
mModelnr = permute(mModelnr,[2 1 3 4]); mUSprf = permute(mUSprf,[2 1 3 4]); 
mUSdc = permute(mUSdc,[2 1 3 4]); mPaR = permute(mPaR,[2 1 3 4]);

FR = zeros(length(Modelnr),length(USprf),length(USdc),length(PaR));
p = gcp;
fprintf('Calculating LIFUS behaviour maps \n');
fprintf('Progress:');
fprintf(['\n' repmat('.',1,numel(FR)) '\n\n']);
parfor ii = 1:numel(FR)          
            fprintf('\b|\n');
            SONICrun(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq),num2str(mUSdc(ii)),num2str(mUSprf(ii)),...
            num2str(Pa2I(mPaR(ii))),'0','0','1','0','0','0',num2str(mModelnr(ii)),'0','0','0',num2str(aBLS),'1');
end

for iMODEL = 1:length(Modelnr)
    for iprf = 1:length(USprf)
        for idc = 1:length(USdc)
            for iPa = 1:length(PaR) 
            ll = load(['APtimes(' MODELstr{Modelnr(iMODEL)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
            num2str(USfreq) ',' num2str(USdc(idc)) ',' num2str(USprf(iprf)) ',' num2str(Pa2I(PaR(iPa))) ')-ES(0,0,1,0,0)-aBLS=('...
            num2str(aBLS) ').mat']);
            APtimes = ll.APtimes(ll.APtimes>=USps&ll.APtimes<=(USps+USpd)); %(s)
            
            FR(iMODEL,iprf,idc,iPa) = length(APtimes)./(USpd);        % (Hz)
            
            delete(['APtimes(' MODELstr{Modelnr(iMODEL)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
            num2str(USfreq) ',' num2str(USdc(idc)) ',' num2str(USprf(iprf)) ',' num2str(Pa2I(PaR(iPa))) ')-ES(0,0,1,0,0)-aBLS=('...
            num2str(aBLS) ').mat']);   
            delete(['Chargevt(' MODELstr{Modelnr(iMODEL)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
            num2str(USfreq) ',' num2str(USdc(idc)) ',' num2str(USprf(iprf)) ',' num2str(Pa2I(PaR(iPa))) ')-ES(0,0,1,0,0)-aBLS=('...
            num2str(aBLS) ').mat']); 
            end
        end
    end
end

figure; set(gcf,'color','w');
for iMODEL = 1:length(Modelnr)
    for iprf = 1:length(USprf)
        cFR(:,:) = FR(iMODEL,iprf,:,:); 
        subplot(2,3,(iMODEL-1)*3+iprf);
        pcolor(100*USdc,10^(-3)*PaR,cFR'); 
        if iprf == 1, ylabel('Amplitude (kPa)');
        elseif iprf == 2, title([MODELstr{Modelnr(iMODEL)} ' neuron']); 
        elseif iprf == 3, cb = colorbar('location','eastoutside'); ylabel(cb,'Firing rate (Hz)');
        end
        xlabel('Duty cycle (%)');
        set(gca,'yscale','log','ydir','normal','colorscale','log');
        colormap(gca,'winter');
        shading flat
    end
end
end
% (b) individual sims
SONICRS = load('SONIC-RS.mat'); SONICLTS = load('SONIC-LTS.mat');
SONICtableRS = SONICRS.SONICtable; SONICtableLTS = SONICLTS.SONICtable;
QmRangeRS = SONICtableRS.QmRange; USPaRangeRS = SONICtableRS.USPaRange; 
USfreqRangeRS = SONICtableRS.USfreqRange; aBLSRangeRS = SONICtableRS.aBLSRange;
QmRangeLTS = SONICtableLTS.QmRange; USPaRangeLTS = SONICtableLTS.USPaRange; 
USfreqRangeLTS = SONICtableLTS.USfreqRange; aBLSRangeLTS = SONICtableLTS.aBLSRange;

Veff4DRS = SONICtableRS.Veff; Veff4DLTS = SONICtableLTS.Veff;

f4VeffRS = @(Qm,USPa,USfreq,aBLS) interpn(QmRangeRS,USPaRangeRS,USfreqRangeRS,aBLSRangeRS,Veff4DRS,Qm,USPa,USfreq,aBLS,'linear');
f2VeffRS = @(Qm,USPa) f4VeffRS(Qm,USPa,USfreq,aBLS); 
f4VeffLTS = @(Qm,USPa,USfreq,aBLS) interpn(QmRangeLTS,USPaRangeLTS,USfreqRangeLTS,aBLSRangeLTS,Veff4DLTS,Qm,USPa,USfreq,aBLS,'linear');
f2VeffLTS = @(Qm,USPa) f4VeffLTS(Qm,USPa,USfreq,aBLS); 
f2Veff = {f2VeffRS,f2VeffLTS};
QmRange = {QmRangeRS,QmRangeLTS};

figure; set(gcf,'color','w');

USdc = {[0.24,0.35],[0.52,0.59],[0.42,0.71];[0.08,0.2],[0.17,0.54],[0.09,0.57]};
PaR = {[200e3,200e3],[392.8e3,392.8e3],[127e3,223.3e3];[54.42e3,127e3],[223.3e3,257.2e3],[41.04e3,146.2e3]};
fprintf('\nCalculating cases in behaviour maps \n'); 
for iMODEL = 1:2
    for iprf = 1:3
    fprintf('Case progress: (%d/6) \n',(iMODEL-1)*3+iprf)
    for icase = 1:2
    subplot(5,3,9*(iMODEL-1)+iprf+(icase-1)*3);
    SONICrun(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc{iMODEL,iprf}(icase)),num2str(USprf(iprf)),...
            num2str(Pa2I(PaR{iMODEL,iprf}(icase))),'0','0','1','0','0','0',num2str(Modelnr(iMODEL)),'0','0','0',num2str(aBLS),'1');
    
    ll = load(['Chargevt(' MODELstr{Modelnr(iMODEL)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
         num2str(USfreq) ',' num2str(USdc{iMODEL,iprf}(icase)) ',' num2str(USprf(iprf)) ',' num2str(Pa2I(PaR{iMODEL,iprf}(icase))) ')-ES(0,0,1,0,0)-aBLS=('...
         num2str(aBLS) ').mat']); 
    
    USprp = (1/USprf(iprf));      % Pulse repetition period (s)
    USstep = @(t) double(mod(t-USps,USprp)<=USdc{iMODEL,iprf}(icase)*USprp).*double(t>=USps&t<=USpd+USps);
    
    PaLine = (PaR{iMODEL,iprf}(icase).*USstep(ll.saveChargeSample(:,1)));
    VeffSample = zeros(size(PaLine));
    for i = 1:length(VeffSample)
    VeffSample(i) = f2Veff{iMODEL}(10^(-5)*ll.saveChargeSample(i,2),PaLine(i));
    end

    yyaxis('left');
    hold on;
    plot(ll.saveChargeSample(:,1),VeffSample); ylim([-100,50]);
    if iprf == 1, ylabel('V_{eff} [mV]'); end
    hold off;
    yyaxis('right');
    hold on;
    plot(ll.saveChargeSample(:,1),ll.saveChargeSample(:,2)); ylim([-100,50]);
    if iprf == 1, ylabel('Q [nC/cm^2]'); end
    if icase == 2, xlabel('Time [s]'); end
    if icase == 1 && iprf == 2, title([MODELstr{Modelnr(iMODEL)} ' neuron']); end
    hold off;

     
    delete(['APtimes(' MODELstr{Modelnr(iMODEL)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc{iMODEL,iprf}(icase)) ',' num2str(USprf(iprf)) ',' num2str(Pa2I(PaR{iMODEL,iprf}(icase))) ')-ES(0,0,1,0,0)-aBLS=('...
    num2str(aBLS) ').mat']);   
    delete(['Chargevt(' MODELstr{Modelnr(iMODEL)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc{iMODEL,iprf}(icase)) ',' num2str(USprf(iprf)) ',' num2str(Pa2I(PaR{iMODEL,iprf}(icase))) ')-ES(0,0,1,0,0)-aBLS=('...
    num2str(aBLS) ').mat']);
    end 
    end
end

%% FIGURE 8 Lemaire et al. (2018)
if isnan(FigurePlot) || FigurePlot == 8
MODELstr = {'RS','FS','LTS'}; lineStyles = {'-','--'}; iColor = 0;
lineColors = {[0.4568 0.4568 0.9906],[0.1092 0.1092 0.9978],[0 0 0.5931],...
    [0.4568 0.9906 0.4568],[0.1092 0.9978 0.1092],[0 0.5931 0]}; 
Legend = {{'16 nm','32 nm','64 nm'},{'20 kHz','500 kHz','4 MHz'}};

Tsim = 2;     % (s)
USpd = 1;   % (s)
USps = 0.5;  % (s)
USprf = 100;          % Hz   

USfreq = {(500e3),[20e3,500e3,4e6]}; % (Hz)
aBLS = {[16e-9,32e-9,64e-9],(32e-9)};  % (m)
USdc = {linspace(0.01,1,20),linspace(0.01,1,20)}; % (-)
Modelnr = {[1,3],[1,3]};

figure; set(gcf,'color','w'); 
fprintf('Calculating sensitivity of threshold-DC plots to frequency and aBLS \n');
for subPlot = 1:2
    iUP = 0; reverseStr = '';
    fprintf('\n Subplot (%d/2) \n',subPlot);
    subplot(1,2,subPlot);
    USfreqRange = USfreq{subPlot}; aBLSRange = aBLS{subPlot};
    USdcRange = USdc{subPlot}; ModelnrRange = Modelnr{subPlot};
    
    ThreshPa = zeros(length(ModelnrRange),length(aBLSRange),length(USfreqRange),length(USdcRange));
for iMODEL = 1:length(ModelnrRange)
    for iaBLS = 1:length(aBLSRange)
        for iFreq = 1:length(USfreqRange)
            for idc = 1:length(USdcRange)
            iUP = iUP+1;    
            Progress = 100*iUP/numel(ThreshPa);  %#ok<*NASGU>
            msg = sprintf('Progress: %3.3f', Progress); 
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));     
                
            SONICrun(num2str(Tsim),'1',num2str(USps),num2str(USpd),num2str(USfreqRange(iFreq)),num2str(USdcRange(idc)),num2str(USprf),...
            '0','0','0','1','0','0','0',num2str(ModelnrRange(iMODEL)),'0',num2str(Pa2I(600e3)),'1',num2str(aBLSRange(iaBLS)),'1');
            ltf = load(['Thresh(' MODELstr{ModelnrRange(iMODEL)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
            num2str(USfreqRange(iFreq)) ',' num2str(USdcRange(idc)) ',' num2str(USprf) ',' num2str(0) ')-ES(0,0,1,0,0)-aBLS=(' ...
            num2str(aBLSRange(iaBLS)) ').mat']);
            ThreshPa(iMODEL,iaBLS,iFreq,idc) = I2Pa(ltf.IIpa);
            delete(['Thresh(' MODELstr{ModelnrRange(iMODEL)} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
            num2str(USfreqRange(iFreq)) ',' num2str(USdcRange(idc)) ',' num2str(USprf) ',' num2str(0) ')-ES(0,0,1,0,0)-aBLS=(' ... 
            num2str(aBLSRange(iaBLS)) ').mat']);            
            end
        end
    end
end
     
     hold on;
     for iaBLS = 1:length(aBLSRange)
         for iFreq = 1:length(USfreqRange)
             iColor = iColor+1;
             for iMODEL = 1:length(ModelnrRange)
             plot(100*USdcRange,permute(10^(-3)*ThreshPa(iMODEL,iaBLS,iFreq,:),[4 1 2 3]),'linestyle',lineStyles{iMODEL},'color',lineColors{iColor});
             end   
         end
     end
     set(gca,'yscale','log','ylim',[10 600]);
     Lines = findobj(gca,'type','line');
     legend(Lines(6:-2:2),Legend{subPlot});
     xlabel('Duty cycle (%)');
     ylabel('Pressure amplitude (kPa)');
     hold off;

end
end

%% FIGURE 9 Lemaire et al. (2018)
if isnan(FigurePlot) || FigurePlot == 9
Tsim = 2.5;     % (s)
USpd = 1;   % (s)
USps = 0.5;  % (s)
USprf = 0; USdc = 1; % (Hz), (-)  
USfreq = 500e3; % (Hz)
aBLS = 32e-9;  % (m)
Modelnr = 9;    
PaR = {I2Pa([(10:10:100) (105:5:120) (121:1:140)]),[I2Pa(70),I2Pa(125),I2Pa(135)]};
figure; set(gcf,'color','w'); 
MODELstr = {};
MODELstr{9} = 'STN';
fprintf('Calculating subthalamic nucleus plots \n');
hf_sub = nan(2,1); hp = nan(2,1); 
load('redCmap.mat'); redCmap = redCmap(1:round(size(redCmap,1)*(0.8)),:);
for subPlot = 1:2
reverseStr = ''; 
iUP = 0;
fprintf('Subplot (%d/2) \n',subPlot);
PaRRange = PaR{subPlot};
hf_sub(subPlot) = figure(subPlot);
hp(subPlot) = uipanel('Parent',hf_sub(subPlot),'Position',[0 0 1 1]);
if subPlot == 1 
subplot(1,1,1,'Parent',hp(subPlot));
end
    
    for iPa = 1:length(PaRRange)
    if subPlot == 2
    subplot(length(PaRRange),1,iPa,'Parent',hp(subPlot));
    end
    iUP = iUP+1;    
    Progress = 100*iUP/numel(PaRRange);  %#ok<*NASGU>
    msg = sprintf('Progress: %3.3f', Progress); 
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));     

    SONICrun(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),...
    num2str(Pa2I(PaRRange(iPa))),'0','0','1','0','0','0',num2str(Modelnr),'0','0','0',num2str(aBLS),'1');

    if subPlot == 1
    hold on;
    ll = load(['APtimes(' MODELstr{Modelnr} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaRRange(iPa))) ')-ES(0,0,1,0,0)-aBLS=('...
    num2str(aBLS) ').mat']);
    APtimes = ll.APtimes(ll.APtimes>=USps&ll.APtimes<=(USps+USpd)); %(s)
    FR = 1./diff(APtimes);        % (Hz)
    plot(10^3*APtimes(1:end-1),FR,'color',redCmap(1+round((size(redCmap,1)-1)*(PaRRange(iPa)/PaRRange(end))),:));
    xlabel('Time [ms]'); ylabel('Firing rate [Hz]');
    rc = colorbar('location','eastoutside','colormap',redCmap);
    set(rc,'TickLabels',10^(-3)*PaRRange(end)*rc.Ticks);
    ylabel(rc,'Pressure [kPa]');
    hold off;
    else
    lr = load(['Chargevt(' MODELstr{Modelnr} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaRRange(iPa))) ')-ES(0,0,1,0,0)-aBLS=('...
    num2str(aBLS) ').mat']);
    
    Qvt = lr.saveChargeSample(:,2); timeline = lr.saveChargeSample(:,1); 
    plot(timeline*10^3,Qvt*10^5);
    title({[num2str(PaRRange(iPa)*10^(-3)) 'kPa'],['(' num2str(Pa2I(PaRRange(iPa))) ' W/m^2)']},'Units','Normalized','Position',[1.05,0.5,1])
    ylabel('Q [nC/cm^2]');
    if iPa == length(PaRRange)
    xlabel('Time [ms]');    
    end
    end
    delete(['APtimes(' MODELstr{Modelnr} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaRRange(iPa))) ')-ES(0,0,1,0,0)-aBLS=('...
    num2str(aBLS) ').mat']);   
    delete(['Chargevt(' MODELstr{Modelnr} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaRRange(iPa))) ')-ES(0,0,1,0,0)-aBLS=('...
    num2str(aBLS) ').mat']); 
    end
end
hf_main = figure(3);
npanels = numel(hp);
hp_sub = nan(1,npanels);
for idx = 1:npanels
hp_sub(idx) = copyobj(hp(idx),hf_main);
set(hp_sub(idx),'Position',[(idx-1)/npanels,0,1/npanels,1]);
close(figure(idx));
end
set(gcf,'color','w');
end