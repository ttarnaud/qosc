% Script will calculate and plot different nanoscale models of the BLS
c = 1515;				% Speed of sound surrounding medium (m/s)
rhol = 1028;			% Density surrounding medium (kg/m^3)
Pa2I = @(Pa) Pa.^2/(2*rhol*c);
I2Pa = @(I) sqrt(2*rhol*c*I);
NICEpath = 'D:\users\ttarnaud\8. Piezoelectric solver\Parallellized functions for HPC calculations';
SONICpath = 'D:\users\ttarnaud\8. Piezoelectric solver\8.4. Lemaire et al. (2018) - SONIC solver';
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
% Model 7: (funPES_nanoMC,proteinMode=0,gateMultip=1): Multi-compartmental NICE model: protein channels full coverage
% Model 8: (funPES_nanoMC,proteinMode=1,gateMultip=1): Multi-compartmental NICE model: no protein
% channels in BLS compartment. Leakage current full coverage. 
% Model 9: (funPES_nanoMC,proteinMode=1,gateMultip=4): Multi-compartmental NICE model: no protein
% channels in BLS compartment. Leakage current full coverage. Gains
% multiplied by gateMultip

Tsim = 2;   % (s)
USpd = 1;   % (s)
USps = 0.5; % (s)
USfreq = 500e3; % (Hz)
USdc = 1; USprf = 0; % (-,Hz)
aBLS = 32e-9;       % (m)
modelnr = 1;

fBLSRange = (0.01:0.01:0.99);  % (-)
fBLSRange = (0.05:0.1:0.95); fprintf('debug mode');

ThreshPa = zeros(length(fBLSRange),9,2);     % [Pa] (fBLS x Model type x threshMode)

for proteinMode = 0:1
for threshMode = 0:1
cd(SONICpath);
if proteinMode == 0, MultipLoop = 1; elseif proteinMode == 1, MultipLoop = [1,4]; end
for gateMultip = MultipLoop
parfor ifBLS=1:length(fBLSRange)
SONICrun_nanoMC(num2str(Tsim),'1',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),'0','0','0','1','0','0','0',num2str(modelnr),'0',num2str(Pa2I(600e3)),'1',num2str(aBLS),num2str(fBLSRange(ifBLS)),num2str(proteinMode),num2str(threshMode),num2str(gateMultip));
SONICrun(num2str(Tsim),'1',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),'0','0','0','1','0','0','0',num2str(modelnr),'0',num2str(Pa2I(600e3)),'1',num2str(aBLS),num2str(fBLSRange(ifBLS)),num2str(proteinMode),num2str(threshMode),num2str(gateMultip));
end
end
for gateMultip = MultipLoop
for ifBLS=1:length(fBLSRange)
LoadStr=['nanoMC-Thresh(RS)-Tsim=' num2str(Tsim) '-US(' num2str(USps) ','...
    num2str(USpd) ',' num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',0)-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS) ')-fBLS=(' num2str(fBLSRange(ifBLS)) ')-(proteinMode,threshMode,gateMultip)=(' num2str(proteinMode) ',' num2str(threshMode) ',' num2str(gateMultip) ').mat']; 
ll = load(LoadStr);
ThreshPa(ifBLS,2*proteinMode+1+2*(gateMultip==4),threshMode+1) = I2Pa(ll.IIpa);
delete(LoadStr);

LoadStr2=['Thresh(RS)-Tsim=' num2str(Tsim) '-US(' num2str(USps) ','...
    num2str(USpd) ',' num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',0)-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS) ')-fBLS=(' num2str(fBLSRange(ifBLS)) ')-(proteinMode,threshMode,gateMultip)=(' num2str(proteinMode) ',' num2str(threshMode) ',' num2str(gateMultip) ').mat']; 
lr = load(LoadStr2);
ThreshPa(ifBLS,2*proteinMode+2+2*(gateMultip==4),threshMode+1) = I2Pa(lr.IIpa);
delete(LoadStr2);
end
end
end
end

markerStyle = {'none','o'};
figure;
hold on;
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