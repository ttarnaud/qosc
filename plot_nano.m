% Script will calculate and plot different nanoscale models of the BLS
clear all;
c = 1515;				% Speed of sound surrounding medium (m/s)
rhol = 1028;			% Density surrounding medium (kg/m^3)
Pa2I = @(Pa) Pa.^2/(2*rhol*c);
I2Pa = @(I) sqrt(2*rhol*c*I);
NICEpath = 'D:\users\ttarnaud\8. Piezoelectric solver\Parallellized functions for HPC calculations';
SONICpath = 'D:\users\ttarnaud\8. Piezoelectric solver\8.4. Lemaire et al. (2018) - SONIC solver';
LoadPathNICE = 'D:\users\ttarnaud\8. Piezoelectric solver\8.4. Lemaire et al. (2018) - SONIC solver\ThreshSims\NICE (HPC)';
LoadPathSONIC = 'D:\users\ttarnaud\8. Piezoelectric solver\8.4. Lemaire et al. (2018) - SONIC solver\ThreshSims\SONIC (LOCAL)';

customLoadPath = 1;         % Bool (if 1, look in loadPath) to load the files 
keepTemps = 1;              % Bool (if 1, keep temporary files, if 0; delete)
performSimuls = 0;          % Bool (if 1, do simulations first)
plotMode = 2;
% Model 1: (SONICrun_nanoMC,proteinMode=0,gateMultip=1) Multi-compartmental SONIC model: protein channels full coverage
% Model 2: (SONICrun,proteinMode=0,gateMultip=1) Point-like SONIC model: protein channels full coverage
% Model 3: (SONICrun_nanoMC,proteinMode=1,gateMultip=1) Multi-compartmental SONIC model: no protein
% channels in BLS compartment. Leakage current full coverage
% Model 4: (SONICrun,proteinMode=1,gateMultip=1) Point-like SONIC model: protein channels partial
% coverage. Leakage current full coverage.
% Model 5: (SONICrun_nanoMC, proteinMode = 1, gateMultip=4) Multi-compartmental SONIC model: no protein channels
% in BLS compartment. Leakage current full coverage. Gains multiplied by gateMultip
% Model 6: (SONICrun, proteinMode = 1, gateMultip=4) Point-like SONIC model: protein channels partial coverage.
% Gains multiplied by gateMultip.
% Model 7-12 : idem for funPES_nanoMC/funPES

Tsim = 0.3;   % (s)
USpd = 0.1;   % (s)
USps = 0.1; % (s)
USfreq = 500e3; % (Hz)
USdc = 1; USprf = 0; % (-,Hz)
aBLS = 32e-9;       % (m)
modelnr = 1;

fBLSRange = (0.01:0.01:0.99);  % (-)
fBLSRange = (0.05:0.1:0.95); 

ThreshPaSONIC = zeros(length(fBLSRange),6,2);     % [Pa] (fBLS x Model type x threshMode)
ThreshPaNICE = zeros(length(fBLSRange),6,2);
TempThreshPa = zeros(length(fBLSRange),6,2);

for i = 1:2
    switch i
        case 1, cd(SONICpath);
        case 2, cd(NICEpath);
    end
for proteinMode = 0:1
for threshMode = 0:1
if proteinMode == 0, MultipLoop = 1; elseif proteinMode == 1, MultipLoop = [1,4]; end
if (performSimuls)
for gateMultip = MultipLoop %#ok<UNRCH>
parfor ifBLS=1:length(fBLSRange)
switch i
    case 1
SONICrun_nanoMC(num2str(Tsim),'1',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),'0','0','0','1','0','0','0',num2str(modelnr),'0',num2str(Pa2I(600e3)),'1',num2str(aBLS),num2str(fBLSRange(ifBLS)),num2str(proteinMode),num2str(threshMode),num2str(gateMultip));
SONICrun(num2str(Tsim),'1',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),'0','0','0','1','0','0','0',num2str(modelnr),'0',num2str(Pa2I(600e3)),'1',num2str(aBLS),num2str(fBLSRange(ifBLS)),num2str(proteinMode),num2str(threshMode),num2str(gateMultip));
    case 2
funPES_nanoMC(num2str(Tsim),'1',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),'0','0','0','1','0','0','0',num2str(modelnr),'0',num2str(Pa2I(600e3)),'1',num2str(aBLS),num2str(fBLSRange(ifBLS)),num2str(proteinMode),num2str(threshMode),num2str(gateMultip));
funPES(num2str(Tsim),'1',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),'0','0','0','1','0','0','0',num2str(modelnr),'0',num2str(Pa2I(600e3)),'1',num2str(aBLS),num2str(fBLSRange(ifBLS)),num2str(proteinMode),num2str(threshMode),num2str(gateMultip));      
end
end
end
end
if customLoadPath
switch i
    case 1, cd(LoadPathSONIC);
    case 2, cd(LoadPathNICE);
end
end
for gateMultip = MultipLoop
for ifBLS=1:length(fBLSRange)
LoadStr=['nanoMC-Thresh(RS)-Tsim=' num2str(Tsim) '-US(' num2str(USps) ','...
    num2str(USpd) ',' num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',0)-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS) ')-fBLS=(' num2str(fBLSRange(ifBLS)) ')-(proteinMode,threshMode,gateMultip)=(' num2str(proteinMode) ',' num2str(threshMode) ',' num2str(gateMultip) ').mat']; 
if ~exist(LoadStr)
TempThreshPa(ifBLS,2*proteinMode+1+2*(gateMultip==4),threshMode+1) = nan;
else
ll = load(LoadStr);
TempThreshPa(ifBLS,2*proteinMode+1+2*(gateMultip==4),threshMode+1) = I2Pa(ll.IIpa);
if (~keepTemps)
delete(LoadStr);
end
end

if ~exist(LoadStr)
TempThreshPa(ifBLS,2*proteinMode+2+2*(gateMultip==4),threshMode+1) = nan;
else
LoadStr2=['Thresh(RS)-Tsim=' num2str(Tsim) '-US(' num2str(USps) ','...
    num2str(USpd) ',' num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',0)-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS) ')-fBLS=(' num2str(fBLSRange(ifBLS)) ')-(proteinMode,threshMode,gateMultip)=(' num2str(proteinMode) ',' num2str(threshMode) ',' num2str(gateMultip) ').mat']; 
lr = load(LoadStr2);
TempThreshPa(ifBLS,2*proteinMode+2+2*(gateMultip==4),threshMode+1) = I2Pa(lr.IIpa);
if (~keepTemps)
delete(LoadStr2);
end
end

end
end
end
end
    switch i
        case 1, ThreshPaSONIC = TempThreshPa;
        case 2, ThreshPaNICE = TempThreshPa;
    end
end

paramsStr = ['(Tsim,USps,USpd,USfreq,USdc,USprf,aBLS,modelnr) = (' num2str(Tsim) ',' num2str(USps) ',' ...
    num2str(USpd) ',' num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(aBLS) ',' num2str(modelnr) ')'];
if (plotMode==1)
for i = 1:2
switch i
    case 1, ThreshPa = ThreshPaSONIC; titleStr = ['SONIC model - ' paramsStr];
    case 2, ThreshPa = ThreshPaNICE; titleStr = ['NICE model - ' paramsStr];
end
markerStyle = {'o','none'};
figure;
hold on;
title(titleStr);
for threshMode = 1:2
plot(100*fBLSRange,ThreshPa(:,1,threshMode)*10^(-3),'linestyle','-','color','b','linewidth',2,'marker',markerStyle{threshMode});
plot(100*fBLSRange,ThreshPa(:,2,threshMode)*10^(-3),'linestyle','-','color',[0.5 0.5 0.5],'linewidth',2,'marker',markerStyle{threshMode});
plot(100*fBLSRange,ThreshPa(:,3,threshMode)*10^(-3),'linestyle','--','color','b','linewidth',2,'marker',markerStyle{threshMode});
plot(100*fBLSRange,ThreshPa(:,4,threshMode)*10^(-3),'linestyle','--','color',[0.5 0.5 0.5],'linewidth',2,'marker',markerStyle{threshMode});
plot(100*fBLSRange,ThreshPa(:,5,threshMode)*10^(-3),'linestyle','-.','color','b','linewidth',2,'marker',markerStyle{threshMode});
plot(100*fBLSRange,ThreshPa(:,6,threshMode)*10^(-3),'linestyle','-.','color',[0.5 0.5 0.5],'linewidth',2,'marker',markerStyle{threshMode});
end
xlabel('Sonophore coverage (%)');
ylabel('Amplitude (kPa)');

LegendPlot(1) = plot(nan,'linestyle','-','color','b');
LegendPlot(2) = plot(nan,'linestyle','-','color',[0.5 0.5 0.5]);
LegendPlot(3) = plot(nan,'linestyle','-','color','k');
LegendPlot(4) = plot(nan,'linestyle','--','color','k');
LegendPlot(5) = plot(nan,'linestyle','-.','color','k');
LegendPlot(6) = plot(nan,'marker','o','color','k');
legendNames = {'Nanoscale multicompartmental SONIC';'Point-like SONIC';'Full protein coverage';'Partial protein coverage';'Partial coverage - modified channel gain';'Excitation at break'};

hold off;
ylim([10,600]);
set(gca,'yscale','log');
set(gcf,'color','w');
legend(LegendPlot,legendNames);
set(findobj('type','axes'),'fontsize',18);
end
end

if (plotMode == 2)
legPos = {[0.0116 0.9334 0.2198 0.0580];[0.2345 0.9062 0.2276 0.0852]};
figure('units','normalized','position',[0 0 1 1]);
for i = 1:2
switch i
case 2, ThreshPa = ThreshPaSONIC; titleStr = ['SONIC model'];
case 1, ThreshPa = ThreshPaNICE; titleStr = ['NICE model'];
end
markerStyle = {'none','none'};
for threshMode = 1:2
subplot(2,2,2*(i-1)+threshMode);
hold on;
plot(100*fBLSRange,ThreshPa(:,1,threshMode)*10^(-3),'linestyle','-','color','b','linewidth',2,'marker',markerStyle{threshMode});
plot(100*fBLSRange,ThreshPa(:,2,threshMode)*10^(-3),'linestyle','-','color',[0.5 0.5 0.5],'linewidth',2,'marker',markerStyle{threshMode});
plot(100*fBLSRange,ThreshPa(:,3,threshMode)*10^(-3),'linestyle','--','color','b','linewidth',2,'marker',markerStyle{threshMode});
plot(100*fBLSRange,ThreshPa(:,4,threshMode)*10^(-3),'linestyle','--','color',[0.5 0.5 0.5],'linewidth',2,'marker',markerStyle{threshMode});
plot(100*fBLSRange,ThreshPa(:,5,threshMode)*10^(-3),'linestyle','-.','color','b','linewidth',2,'marker',markerStyle{threshMode});
plot(100*fBLSRange,ThreshPa(:,6,threshMode)*10^(-3),'linestyle','-.','color',[0.5 0.5 0.5],'linewidth',2,'marker',markerStyle{threshMode});
if (i==2)
if threshMode==1
LegendPlot(1) = plot(nan,'linestyle','-','color','b');
LegendPlot(2) = plot(nan,'linestyle','-','color',[0.5 0.5 0.5]);
legendNames = {'Nanoscale multicompartmental SONIC';'Point-like SONIC'};
else
LegendPlot(1) = plot(nan,'linestyle','-','color','k');
LegendPlot(2) = plot(nan,'linestyle','--','color','k');
LegendPlot(3) = plot(nan,'linestyle','-.','color','k');
legendNames = {'Full protein coverage';'Partial protein coverage';'Partial coverage - modified channel gain'};
end
leg(threshMode) = legend(LegendPlot,legendNames,'position',legPos{threshMode});
end
if (threshMode == 1)
ylabel('Amplitude (kPa)');
end
if (i==2)
xlabel('Sonophore coverage (%)');
end
hold off;
ylim([10,600]);
set(gca,'yscale','log');
set(gcf,'color','w');
end
tit(i) = title(titleStr,'Position',[60.7271 741.3792 0]);
set(findobj('type','axes'),'fontsize',18);
end
set(tit(2),'Position',[-9.8794 815.4413 0]);

end