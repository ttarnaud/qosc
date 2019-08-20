% Plot results in order to reproduce results in Lemaire et al.
c = 1515;				% Speed of sound surrounding medium (m/s)
rhol = 1028;			% Density surrounding medium (kg/m^3)
Pa2I = @(Pa) Pa.^2/(2*rhol*c);
I2Pa = @(I) sqrt(2*rhol*c*I);
NICEpath = 'D:\users\ttarnaud\8. Piezoelectric solver\Parallellized functions for HPC calculations';
SONICpath = 'D:\users\ttarnaud\8. Piezoelectric solver\8.4. Lemaire et al. (2018) - SONIC solver';

%% FIGURE 1
MODEL = 1;
switch MODEL
    case 1, ModelStr = 'RS';
    case 2, ModelStr = 'FS';
    case 3, ModelStr = 'LTS';
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

cd(SONICpath); 
fprintf('Continuous wave stimulation: charge trace plots \n');
for i = 1:3
fprintf('Calculation of subplot (%u/3) \n',i);
subplot(7,3,i);
hold on;
for iFreq = 1:length(USfreq{i})
for iaBLS = 1:length(aBLS{i})
SONICrun(num2str(Tsim),'1',num2str(USps),num2str(USpd),num2str(USfreq{i}(iFreq)),num2str(USdc),num2str(USprf),...
'0','0','0','1','0','0','0',num2str(MODEL),'0','0','0',num2str(aBLS{i}(iaBLS)));
lt = load(['Thresh(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq{i}(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(0) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS{i}(iaBLS)) ').mat']);
PaThreshSONIC = I2Pa(lt.IIpa);
PaThRange = PaThreshSONIC+RelPa{i};
delete(['Thresh(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq{i}(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(0) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS{i}(iaBLS)) ').mat']);

% NICE around threshold
cd(NICEpath);
for iPa = 1:length(PaThRange)
funPESa(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq{i}(iFreq)),num2str(USdc),num2str(USprf),...
num2str(Pa2I(PaThRange(iPa))),'0','0','1','0','0','0',num2str(MODEL),'0','0','0',num2str(aBLS{i}(iaBLS)));
end
for iPa = 1:length(PaThRange)
lt2 = load(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq{i}(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaThRange(iPa))) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS{i}(iaBLS)) ').mat']);
Qvt = lt2.saveChargeSample(:,2); timeline = lt2.saveChargeSample(:,1);
iNICEcolor = iNICEcolor+1;
plot(10^3*timeline,Qvt,'linestyle','-','color',SONICcolors{iNICEcolor},'linewidth',2);
delete(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq{i}(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaThRange(iPa))) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS{i}(iaBLS)) ').mat']);
delete(['APtimes(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq{i}(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaThRange(iPa))) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS{i}(iaBLS)) ').mat']);
end
% SONIC around threshold
cd(SONICpath);
for iPa = 1:length(PaThRange)
SONICrun(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq{i}(iFreq)),num2str(USdc),num2str(USprf),...
num2str(Pa2I(PaThRange(iPa))),'0','0','1','0','0','0',num2str(MODEL),'0','0','0',num2str(aBLS{i}(iaBLS)));
end
for iPa = 1:length(PaThRange)
lt3 = load(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq{i}(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaThRange(iPa))) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS{i}(iaBLS)) ').mat']);
Qvt = lt3.saveChargeSample(:,2); timeline = lt3.saveChargeSample(:,1); 
iSONICcolor = iSONICcolor+1;
plot(10^3*timeline,Qvt,'linestyle','--','color',SONICcolors{iSONICcolor});
delete(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq{i}(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaThRange(iPa))) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS{i}(iaBLS)) ').mat']);
delete(['APtimes(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq{i}(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaThRange(iPa))) ')-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS{i}(iaBLS)) ').mat']);
end
end
end
if (i==1)
    ylabel('Q [nC/cm^2]');
end
xlabel('Time [ms]');
end

% Fig. 1(a-c) bottom
USfreq = {(500e3), logspace(log10(20e3),log10(4e6),7),(500e3)}; % (Hz)
aBLS = {(32e-9),(32e-9),linspace(16e-9,64e-9,5)};  % (m)
RelPa = {'',(20e3),(20e3)};  % relative to threshold (Pa)
PaR = {1e3*logspace(log10(50),log10(150),4),'',''};   % Pressure range (Pa) 

fprintf('Continuous wave stimulation: latency, firing rate and spike amplitude plots \n');
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
    '0','0','0','1','0','0','0',num2str(MODEL),'0','0','0',aBLSRange(iaBLS));
    ltf = load(['Thresh(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreqRange(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(0) ')-ES(0,0,1,0,0).mat-aBLS=(' ...
    num2str(aBLSRange(iaBLS)) ').mat']);
    ThreshPa(iFreq,iaBLS) = I2Pa(ltf.IIpa);
    delete(['Thresh(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreqRange(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(0) ')-ES(0,0,1,0,0).mat-aBLS=(' ... 
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
    for iFreq = length(USfreqRange)
        for iaBLS = length(aBLSRange)
            for iPaR = 1:length(RelPaRange)
            appPa = ThreshPa(iFreq,iaBLS)+RelPaRange(iPaR);
            if pathNum == 1
            SONICrun(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreqRange(iFreq)),num2str(USdc),num2str(USprf),...
            num2str(Pa2I(appPa)),'0','0','1','0','0','0',num2str(MODEL),'0','0','0',aBLSRange(iaBLS));
            elseif pathNum == 2
            funPESa(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreqRange(iFreq)),num2str(USdc),num2str(USprf),...
            num2str(Pa2I(appPa)),'0','0','1','0','0','0',num2str(MODEL),'0','0','0',aBLSRange(iaBLS));
            end
            end
        end
    end
    
    for iFreq = length(USfreqRange)
        for iaBLS = length(aBLSRange)
            for iPaR = 1:length(RelPaRange)
            appPa = ThreshPa(iFreq,iaBLS)+RelPaRange(iPaR);

            ll = load(['APtimes(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
            num2str(USfreqRange(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(appPa)) ')-ES(0,0,1,0,0)-aBLS=('...
            num2str(aBLSRange(iaBLS)) ').mat']);
            lr = load(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
            num2str(USfreqRange(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(appPa)) ')-ES(0,0,1,0,0)-aBLS=(' ...
            num2str(aBLSRange(iaBLS)) ').mat']);
            APtimes = ll.APtimes(ll.APtimes>=USps&ll.APtimes<=(USps+USpd)); %(s)
            
            Latency(pathNum,USfreqRange(iFreq),aBLSRange(iaBLS),RelPaRange(iPaR)) = 10^3*(APtimes(1)-USps);    % (ms)
            FR(pathNum,USfreqRange(iFreq),aBLSRange(iaBLS),RelPaRange(iPaR)) = length(APtimes)./(USpd);        % (Hz)
            
            Qvt = lr.saveChargeSample(:,2); timeline = lr.saveChargeSample(:,1);  % (nC/cm^2) and (s)
            QvtCrop = Qvt(timeline>=APtimes(2)&timeline<=APtimes(end-1));  % Crop outer lying action potentials to determine amplitude
            timelineCrop = timeline(timeline>=APtimes(2)&timeline<=APtimes(end-1));
            spA(pathNum,USfreqRange(iFreq),aBLSRange(iaBLS),RelPaRange(iPaR)) = (max(QvtCrop)-min(QvtCrop));       % (nC/cm^2)
            
            delete(['APtimes(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
            num2str(USfreqRange(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(appPa)) ')-ES(0,0,1,0,0)-aBLS=('...
            num2str(aBLSRange(iaBLS)) ').mat']);
            delete(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
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
    case 3, xaxis = 1e-9*aBLS{3}; xlab = 'Sonophore radius (nm)';

end
subplot(7,3,3+i);
hold on;
plot(xaxis,LatencyNICE,'linestyle','-','color',[0.5 0.5 0.5],'linewidth',2,'xlim',[xaxis(1),xaxis(end)]);
plot(xaxis,LatencySONIC,'linestyle','--','color','k','marker','o','xlim',[xaxis(1),xaxis(end)]);
ylabel('Latency [ms]');
hold off;
subplot(7,3,6+i);
hold on;
plot(xaxis,FRNICE,'linestyle','-','color',[0.5 0.5 0.5],'linewidth',2,'xlim',[xaxis(1),xaxis(end)]);
plot(xaxis,FRSONIC,'linestyle','--','color','k','marker','o','xlim',[xaxis(1),xaxis(end)]);
ylabel('Firing rate [Hz]');
hold off;
subplot(7,3,9+i);
hold on;
plot(xaxis,spANICE,'linestyle','-','color',[0.5 0.5 0.5],'linewidth',2,'xlim',[xaxis(1),xaxis(end)]);
plot(xaxis,spASONIC,'linestyle','--','color','k','marker','o','xlim',[xaxis(1),xaxis(end)]);
ylabel('Spike amp. [nC/cm^2]');
xlabel(xlab);
hold off;
end
set(findobj('type','axes'),'fontsize',18);
set(findobj('type','axes'),'box','off');
