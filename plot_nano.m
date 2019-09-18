% Script will calculate and plot different nanoscale models of the BLS

% Model 1: (SONICrun_nanoMC) Multi-compartmental SONIC model: protein channels full coverage
% Model 2: (SONICrun) Point-like SONIC model: protein channels full coverage
% Model 3: (SONICrun_nanoMC) Multi-compartmental SONIC model: no protein
% channels in BLS compartment. Leakage current full coverage
% Model 4: (SONICrun) Point-like SONIC model: protein channels partial
% coverage. Leakage current full coverage.
% Model 5: (funPES_nanoMC): Multi-compartmental NICE model: protein channels full coverage
% Model 6: (funPES_nanoMC): Multi-compartmental NICE model: no protein
% channels in BLS compartment. Leakage current full coverage. 
Tsim = 2;   % (s)
USpd = 1;   % (s)
USps = 0.5; % (s)

fBLSRange = (0.01:0.01:0.99);  % (-)
ThreshPa = zeros(length(fBLSRange),2);

parfor ifBLS=1:length(fBLSRange)
SONICrun_nanoMC(num2str(Tsim),'1',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),'0','0','0','1','0','0','0',num2str(modelnr),'0',num2str(Pa2I(600e3)),'1',num2str(aBLS),num2str(fBLSRange(ifBLS)));
SONICrun(num2str(Tsim),'1',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),'0','0','0','1','0','0','0',num2str(modelnr),'0',num2str(Pa2I(600e3)),'1',num2str(aBLS),num2str(fBLSRange(ifBLS)));
end
for ifBLS=1:length(fBLSRange)
LoadStr=['nanoMC-Thresh(RS)-Tsim=' num2str(Tsim) '-US(' num2str(USps) ','...
    num2str(USpd) ',' num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',0)-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS) ')-fBLS=(' num2str(fBLSRange(ifBLS)) ').mat']; 
ll = load(LoadStr);
ThreshPa(ifBLS,1) = I2Pa(ll.IIpa);
delete(LoadStr);

LoadStr2=['Thresh(RS)-Tsim=' num2str(Tsim) '-US(' num2str(USps) ','...
    num2str(USpd) ',' num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',0)-ES(0,0,1,0,0)-aBLS=(' ...
    num2str(aBLS) ')-fBLS=(' num2str(fBLSRange(ifBLS)) ').mat']; 
lr = load(LoadStr2);
ThreshPa(ifBLS,2) = I2Pa(lr.IIpa);
delete(LoadStr2);
end
figure;
hold on;
plot(100*fBLSRange,ThreshPa(:,1)*10^(-3),'linestyle','-','color','b','linewidth',2);
plot(100*fBLSRange,ThreshPa(:,2)*10^(-3),'linestyle','-','color',[0.5 0.5 0.5],'linewidth',2);
xlabel('Sonophore coverage (%)');
ylabel('Amplitude (kPa)');
hold off;
ylim([10,600]);
set(gca,'yscale','log');
set(gcf,'color','w');