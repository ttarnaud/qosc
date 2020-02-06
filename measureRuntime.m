clear all;
c = 1515;				% Speed of sound surrounding medium (m/s)
rhol = 1028;			% Density surrounding medium (kg/m^3)
Pa2I = @(Pa) Pa.^2/(2*rhol*c);
I2Pa = @(I) sqrt(2*rhol*c*I);
NICEpath = 'D:\users\ttarnaud\8. Piezoelectric solver\Parallellized functions for HPC calculations';
NICEloadpath = 'D:\users\ttarnaud\8. Piezoelectric solver\8.4. Lemaire et al. (2018) - SONIC solver\measureRuntime\NICE';
SONICpath = 'D:\users\ttarnaud\8. Piezoelectric solver\8.4. Lemaire et al. (2018) - SONIC solver';
SONICloadpath = 'D:\users\ttarnaud\8. Piezoelectric solver\8.4. Lemaire et al. (2018) - SONIC solver\measureRuntime\SONIC';
addpath(genpath(pwd));

Tsim = 0.1;   % (s)
USpd = 0.1;   % (s)
USps = 0; % (s)
USfreq = 500e3; % (Hz)
USdc = 1; USprf = 0; % (-,Hz)
aBLS = 32e-9;       % (m)
modelnr = 1;
fBLS = 0.75;
proteinMode = 0; threshMode = 0; gateMultip = 1;

USPaRange = logspace(log10(30e3),log10(600000),5);  % Ultrasonic pressure (Pa)
fBLSRange = linspace(0.5,0.95,5);

SONICorNICE = 1; compiled = 0;
noPlot = 0; noSim = 1;
%% --- SONIC 
if SONICorNICE == 1 || SONICorNICE == 0
if ~compiled
cd(SONICpath);
end
if noSim
cd(SONICloadpath);
load('SimTimeSONIC.mat');
TTimeSONIC = SimTimeSONIC.pointN; TTimeSONICnanoMC = SimTimeSONIC.nanoMC;
else
p = gcp;
nLoops = 100;
TTimeSONICtemp = zeros(length(USPaRange),length(fBLSRange)); TTimeSONICnanoMCtemp = zeros(length(USPaRange),length(fBLSRange));
TTimeSONICnanoMC = zeros(length(USPaRange),length(fBLSRange),nLoops);
TTimeSONIC = zeros(length(USPaRange),length(fBLSRange),nLoops);
for ll = 1:nLoops
for i = 1:length(USPaRange)
    for j = 1:length(fBLSRange)
USPa = USPaRange(i);  fBLS = fBLSRange(j);
futureTTimeSONIC(i,j) = parfeval(p,@SONICrun,1,num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),num2str(Pa2I(USPa)),'0','0','1','0','0','0',num2str(modelnr),'0',num2str(Pa2I(600e3)),'1',num2str(aBLS),num2str(fBLS),num2str(proteinMode),num2str(threshMode),num2str(gateMultip)); %#ok<SAGROW>
    end
end
for i = 1:length(USPaRange)
    for j = 1:length(fBLSRange)
[futInd,TTimeFut] = fetchNext(futureTTimeSONIC);
TTimeSONICtemp(ind2sub([length(USPaRange),length(fBLSRange)],futInd)) = TTimeFut;
    end
end
TTimeSONIC(:,:,ll) = TTimeSONICtemp;
for i = 1:length(USPaRange)
    for j = 1:length(fBLSRange)
USPa = USPaRange(i);  fBLS = fBLSRange(j);
futureTTimeSONICnanoMC(i,j) = parfeval(p,@SONICrun_nanoMC,1,num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),num2str(Pa2I(USPa)),'0','0','1','0','0','0',num2str(modelnr),'0',num2str(Pa2I(600e3)),'1',num2str(aBLS),num2str(fBLS),num2str(proteinMode),num2str(threshMode),num2str(gateMultip)); %#ok<SAGROW>
    end
end
for i = 1:length(USPaRange)
    for j = 1:length(fBLSRange)
[futInd,TTimeFut] = fetchNext(futureTTimeSONICnanoMC); %#ok<SAGROW>
TTimeSONICnanoMCtemp(ind2sub([length(USPaRange),length(fBLSRange)],futInd)) = TTimeFut;
    end
end
TTimeSONICnanoMC(:,:,ll) = TTimeSONICnanoMCtemp; 
end
SimTimeSONIC.pointN = TTimeSONIC; SimTimeSONIC.nanoMC = TTimeSONICnanoMC;
save('SimTimeSONIC.mat','SimTimeSONIC');
end
if ~noPlot
% Plotting the results
preStr = {'','nanoMC-'}; %#ok<*UNRCH>
for ii = 1:2
switch ii
    case 1, TTime = mean(TTimeSONIC,3); titleStr = 'SONIC';         % [s]
    case 2, TTime = mean(TTimeSONICnanoMC,3); titleStr = 'SONIC (MC)'; % [s]
end
    
figure('units','normalized','position',[0 0 1 1]); nwdth = 0.7; nhght = 0.7; % Normalized width and height
hold on;
for i = 1:length(USPaRange)
    for j = 1:length(fBLSRange)
    USPa = USPaRange(end+1-i); fBLS = fBLSRange(j);
    h(i,j) = subplot(length(USPaRange)+3,2*length(fBLSRange)+2,(i)*(2*length(fBLSRange)+2)+j+1);
    loadStr = [preStr{ii} 'Chargevt(RS)-Tsim=0.1-US(0,0.1,500000,1,0,'  num2str(Pa2I(USPa)) ')-ES(0,0,1,0,0)-aBLS=(3.2e-08)-fBLS=(' num2str(fBLS) ')-(proteinMode,threshMode,gateMultip)=(0,0,1).mat'];
    ll = load(loadStr); saveChargeSample = ll.saveChargeSample;
    switch ii
        case 1, plot(1e3*saveChargeSample(:,1),saveChargeSample(:,2));
        case 2, plot(1e3*saveChargeSample(:,1),fBLS*saveChargeSample(:,2)+(1-fBLS)*saveChargeSample(:,3));
    end
    if i ~= 5 || j ~= 1 
    set(gca,'visible','off');
    else
    box off; xlabel('Time [ms]','position',[56.2500 -149.0909 -1]); ylabel('Charge [nC/cm^2]');
    end
    set(gca,'ylim',[-100 50]);
    end
    set(gcf,'color','w');
end
for i = 1:length(USPaRange)
    for j = 1:length(fBLSRange)
    set(h(i,j),'position',[(1-nwdth)/(2*length(fBLSRange)+2)/2 + (j)/(2*length(fBLSRange)+1), (1-nhght)/(length(USPaRange)+3)/2 + 1-(i+1)/(length(USPaRange)+3), nwdth/(2*length(fBLSRange)+2), nhght/(length(USPaRange)+3)]);
    end
end

colormap(cm_viridis(100));
subindc = ((length(fBLSRange)+3:2*length(fBLSRange)+2)+(1:length(USPaRange))'*(2*length(fBLSRange)+2)); subindc = sort(subindc(:));
subplot(length(USPaRange)+3,2*length(fBLSRange)+2,subindc); surf(fBLSRange,1e-3*USPaRange,TTime,'facecolor','interp','edgecolor','interp'); ylabel('P_{US} [kPa]'); xlabel('f_{BLS} [-]');
ylim(1e-3*[USPaRange(1),USPaRange(end)]); xlim([fBLSRange(1),fBLSRange(end)]); view([0 90]);
cb = colorbar('position',[0.5911 0.8386 0.3141 0.0198],'orientation','horizontal');
set(get(cb,'title'),'String','Simulation time [s]');

supX = [fBLSRange(1)-fBLSRange(2)+fBLSRange(1),fBLSRange,fBLSRange(end)+fBLSRange(2)-fBLSRange(1)];
supXlbl = num2cell(supX);
supXlbl{1} = []; supXlbl{end} = [];
supY = USPaRange;
supY = [supY(1)/(supY(3)/supY(2)), supY, supY(end)*supY(3)/supY(2)];
supYlbl = num2cell(1e-3*supY);
supYlbl{1} = []; supYlbl{end} = [];
supL = suplabel('f_{BLS} [-]','x',[0.0479 0.1972 0.5437 0.7513]);
set(get(supL,'ylabel'),'String','P_{US} [kPa]'); set(get(supL,'ylabel'),'rotation',0,'position',[0.3853 14.0240]);
set(supL,'yscale','lin','ytick',log(supY),'ylim',[log(supY(1)), log(supY(end))]);
set(supL,'yticklabel',cellfun(@(X) num2str(X,3),supYlbl,'UniformOutput',false))
set(supL,'visible','on');
set(supL,'xscale','lin','xtick',supX,'xlim',[supX(1),supX(end)]);
set(supL,'xticklabel',cellfun(@(X) num2str(X,3),supXlbl,'UniformOutput',false))

ttl = text(681.7500,235.7273,titleStr,'fontsize',25,'fontweight','bold');

set(findobj('type','axes'),'fontsize',18,'fontweight','bold');
end
end
end
%% --- NICE
if SONICorNICE == 2 || SONICorNICE == 0
if ~compiled
cd(NICEpath);
end
if noSim
cd(NICEloadpath);
load('SimTimeNICE');
TTimeNICE = SimTimeNICE.pointN;
TTimeNICEnanoMC = SimTimeNICE.nanoMC;
else
p = gcp;
nLoops = 1;
TTimeNICEtemp =  zeros(length(USPaRange),length(fBLSRange));
TTimeNICEnanoMCtemp =  zeros(length(USPaRange),length(fBLSRange));
TTimeNICEnanoMC = zeros(length(USPaRange),length(fBLSRange),nLoops);
TTimeNICE = zeros(length(USPaRange),length(fBLSRange),nLoops);
for ll = 1:nLoops
for i = 1:length(USPaRange)
    for j = 1:length(fBLSRange)
USPa = USPaRange(i);  fBLS = fBLSRange(j);
futureTTimeNICE(i,j) = parfeval(p,@funPES,1,num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),num2str(Pa2I(USPa)),'0','0','1','0','0','0',num2str(modelnr),'0',num2str(Pa2I(600e3)),'1',num2str(aBLS),num2str(fBLS),num2str(proteinMode),num2str(threshMode),num2str(gateMultip)); %#ok<SAGROW>
    end
end
for i = 1:length(USPaRange)
    for j = 1:length(fBLSRange)
[futInd,TTimeFut] = fetchNext(futureTTimeNICE);
TTimeNICEtemp(ind2sub([length(USPaRange),length(fBLSRange)],futInd)) = TTimeFut;
    end
end
TTimeNICE(:,:,ll) = TTimeNICEtemp;
for i = 1:length(USPaRange)
    for j = 1:length(fBLSRange)
USPa = USPaRange(i);  fBLS = fBLSRange(j);
futureTTimeNICEnanoMC(i,j) = parfeval(p,@funPES_nanoMC,1,num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),num2str(Pa2I(USPa)),'0','0','1','0','0','0',num2str(modelnr),'0',num2str(Pa2I(600e3)),'1',num2str(aBLS),num2str(fBLS),num2str(proteinMode),num2str(threshMode),num2str(gateMultip)); %#ok<SAGROW>
    end
end
for i = 1:length(USPaRange)
    for j = 1:length(fBLSRange)
[futInd,TTimeFut] = fetchNext(futureTTimeNICEnanoMC); %#ok<SAGROW>
TTimeNICEnanoMCtemp(ind2sub([length(USPaRange),length(fBLSRange)],futInd)) = TTimeFut;
    end
end
TTimeNICEnanoMC(:,:,ll) = TTimeNICEnanoMCtemp;
end
SimTimeNICE.pointN = TTimeNICE; SimTimeNICE.nanoMC = TTimeNICEnanoMC;
save('SimTimeNICE.mat','SimTimeNICE');
end
if ~noPlot
% Plotting the results
preStr = {'','nanoMC-'}; %#ok<UNRCH>
for ii = 1:2
switch ii
    case 1, TTime = mean(TTimeNICE,3)/60^2; titleStr = 'NICE';              % [h]              
    case 2, TTime = mean(TTimeNICEnanoMC,3)/60^2; titleStr = 'NICE (MC)';   % [h]     
end
    
figure('units','normalized','position',[0 0 1 1]); nwdth = 0.7; nhght = 0.7; % Normalized width and height
hold on;
for i = 1:length(USPaRange)
    for j = 1:length(fBLSRange)
    USPa = USPaRange(end+1-i); fBLS = fBLSRange(j);
    h(i,j) = subplot(length(USPaRange)+3,2*length(fBLSRange)+2,(i)*(2*length(fBLSRange)+2)+j+1);
    loadStr = [preStr{ii} 'Chargevt(RS)-Tsim=0.1-US(0,0.1,500000,1,0,'  num2str(Pa2I(USPa)) ')-ES(0,0,1,0,0)-aBLS=(3.2e-08)-fBLS=(' num2str(fBLS) ')-(proteinMode,threshMode,gateMultip)=(0,0,1).mat'];
    ll = load(loadStr); saveChargeSample = ll.saveChargeSample;
    switch ii
        case 1,  plot(1e3*saveChargeSample(:,1),saveChargeSample(:,2));
        case 2,  plot(1e3*saveChargeSample(:,1),fBLS*saveChargeSample(:,2)+(1-fBLS)*saveChargeSample(:,3));
    end
    if i ~= 5 || j ~= 1 
    set(gca,'visible','off');
    else
    box off; xlabel('Time [ms]','position',[56.2500 -149.0909 -1]); ylabel('Charge [nC/cm^2]');
    end
    set(gca,'ylim',[-100 50]);
    end
    set(gcf,'color','w');
end
for i = 1:length(USPaRange)
    for j = 1:length(fBLSRange)
    set(h(i,j),'position',[(1-nwdth)/(2*length(fBLSRange)+2)/2 + (j)/(2*length(fBLSRange)+1), (1-nhght)/(length(USPaRange)+3)/2 + 1-(i+1)/(length(USPaRange)+3), nwdth/(2*length(fBLSRange)+2), nhght/(length(USPaRange)+3)]);
    end
end

colormap(cm_viridis(100));
subindc = ((length(fBLSRange)+3:2*length(fBLSRange)+2)+(1:length(USPaRange))'*(2*length(fBLSRange)+2)); subindc = sort(subindc(:));
subplot(length(USPaRange)+3,2*length(fBLSRange)+2,subindc); surf(fBLSRange,1e-3*USPaRange,TTime,'facecolor','interp','edgecolor','interp'); ylabel('P_{US} [kPa]'); xlabel('f_{BLS} [-]');
ylim(1e-3*[USPaRange(1),USPaRange(end)]); xlim([fBLSRange(1),fBLSRange(end)]); view([0 90]);
cb = colorbar('location','westoutside');
set(get(cb,'title'),'String','Simulation time [h]');

supX = [fBLSRange(1)-fBLSRange(2)+fBLSRange(1),fBLSRange,fBLSRange(end)+fBLSRange(2)-fBLSRange(1)];
supXlbl = num2cell(supX);
supXlbl{1} = []; supXlbl{end} = [];
supY = USPaRange;
supY = [supY(1)/(supY(3)/supY(2)), supY, supY(end)*supY(3)/supY(2)];
supYlbl = num2cell(1e-3*supY);
supYlbl{1} = []; supYlbl{end} = [];
supL = suplabel('f_{BLS} [-]','x',[0.0479 0.1972 0.5437 0.7513]);
set(get(supL,'ylabel'),'String','P_{US} [kPa]'); set(get(supL,'ylabel'),'rotation',0,'position',[0.3853 14.0240]);
set(supL,'yscale','lin','ytick',log(supY),'ylim',[log(supY(1)), log(supY(end))]);
set(supL,'yticklabel',cellfun(@(X) num2str(X,3),supYlbl,'UniformOutput',false))
set(supL,'visible','on');
set(supL,'xscale','lin','xtick',supX,'xlim',[supX(1),supX(end)]);
set(supL,'xticklabel',cellfun(@(X) num2str(X,3),supXlbl,'UniformOutput',false))

ttl = text(681.7500,235.7273,titleStr,'fontsize',25,'fontweight','bold');

set(findobj('type','axes'),'fontsize',18,'fontweight','bold');
end
end
end