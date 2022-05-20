c = 1515;				% Speed of sound surrounding medium (m/s)
rhol = 1028;			% Density surrounding medium (kg/m^3)
Pa2I = @(Pa) Pa.^2/(2*rhol*c);
I2Pa = @(I) sqrt(2*rhol*c*I);
NICEpath = 'D:\users\ttarnaud\8. Piezoelectric solver\Parallellized functions for HPC calculations';
SONICpath = 'D:\users\ttarnaud\8. Piezoelectric solver\8.4. Lemaire et al. (2018) - SONIC solver';

PlotN = 2;
Simul = 0;
keepTemps = 1;
postProc = 0;

if (PlotN == 1)
% Plot 1
Tsim = 2.5;     % (s)
USpd = 1;   % (s)
USps = 0.5;  % (s)
USprf = 0; USdc = 1; % (Hz), (-)  
USfreq = 500e3; % (Hz)
aBLS = 32e-9;  % (m)
Modelnr = 9;    
PaR = {I2Pa([(10:10:100) (105:5:120) (121:1:140)]),[I2Pa(70),I2Pa(125),I2Pa(135)]};
ESps = 1; 
ESpd = 20e-3;
ESdc = 1; ESprf = 0;
ESi = -0.2;
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
    
    if ~~(Simul)
    SONICrun(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),...
    num2str(Pa2I(PaRRange(iPa))), num2str(ESps),num2str(ESpd),num2str(ESdc),num2str(ESprf),num2str(ESi),'0',num2str(Modelnr),'0','0','0',num2str(aBLS),'1');
    end
    
    if subPlot == 1
    hold on;
    ll = load(['APtimes(' MODELstr{Modelnr} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaRRange(iPa))) ')-ES(' num2str(num2str(ESps)) ',' ...
    num2str(ESpd) ',' num2str(ESdc) ',' num2str(ESprf) ',' num2str(ESi) ')-aBLS=('...
    num2str(aBLS) ')-fBLS=(1).mat']);
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
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaRRange(iPa))) ')-ES(' num2str(ESps) ',' ...
    num2str(ESpd) ',' num2str(ESdc) ',' num2str(ESprf) ',' num2str(ESi) ')-aBLS=('...
    num2str(aBLS) ')-fBLS=(1).mat']);
    
    Qvt = lr.saveChargeSample(:,2); timeline = lr.saveChargeSample(:,1); 
    plot(timeline*10^3,Qvt);
    title({[num2str(PaRRange(iPa)*10^(-3)) 'kPa'],['(' num2str(Pa2I(PaRRange(iPa))) ' W/m^2)']},'Units','Normalized','Position',[1.05,0.5,1])
    ylabel('Q [nC/cm^2]');
    if iPa == length(PaRRange)
    xlabel('Time [ms]');    
    end
    end
    if ~(keepTemps)
    delete(['APtimes(' MODELstr{Modelnr} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaRRange(iPa))) ')-ES(' num2str(ESps) ',' ...
    num2str(ESpd) ',' num2str(ESdc) ',' num2str(ESprf) ',' num2str(ESi) ')-aBLS=('...
    num2str(aBLS) ')-fBLS=(1).mat']);   
    delete(['Chargevt(' MODELstr{Modelnr} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaRRange(iPa))) ')-ES(' num2str(ESps) ',' ...
    num2str(ESpd) ',' num2str(ESdc) ',' num2str(ESprf) ',' num2str(ESi) ')-aBLS=('...
    num2str(aBLS) ')-fBLS=(1).mat']); 
    end
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

% Plot 2
if (PlotN == 2)
Tsim = 2.5;     % (s)
USpd = 1;   % (s)
USps = 0.5;  % (s)
USprf = 100; USdc = 0.05; % (Hz), (-)  
USfreq = 500e3; % (Hz)
aBLS = 32e-9;  % (m)
Modelnr = 9;    
%PaR = {I2Pa([(10:10:100) (105:5:120) (121:1:140)]),[I2Pa(20),I2Pa(132),I2Pa(140)]};
%PaR = {I2Pa([(10:10:100) (105:5:120) (121:1:140)]),I2Pa([(10:10:100) (105:5:120) (121:1:140)])};
%PaR = {I2Pa((0:50:500));[]};
PaR = {I2Pa((300));I2Pa((300))};
ESps = 1; 
ESpd = 20e-3;
ESdc = 1; ESprf = 0;
ESi = 0; 
MODELstr = {};
MODELstr{9} = 'STN';
fprintf('Calculating subthalamic nucleus plots \n');
 
load('redCmap.mat'); redCmap = redCmap(1:round(size(redCmap,1)*(0.8)),:);
for subPlot = 1:2
fprintf('Subplot (%d/2) \n',subPlot);
PaRRange = PaR{subPlot};
    
    if ~~(Simul)
    parfor iPa = 1:length(PaRRange)    
    SONICrun(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),...
    num2str(Pa2I(PaRRange(iPa))), num2str(ESps),num2str(ESpd),num2str(ESdc),num2str(ESprf),num2str(ESi),'0',num2str(Modelnr),'0','0','0',num2str(aBLS),'1');
    end
    end
    
    figure; set(gcf,'color','w');
    for iPa = 1:length(PaRRange)
    if subPlot == 1
    hold on;
    ll = load(['APtimes(' MODELstr{Modelnr} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaRRange(iPa))) ')-ES(' num2str(num2str(ESps)) ',' ...
    num2str(ESpd) ',' num2str(ESdc) ',' num2str(ESprf) ',' num2str(ESi) ')-aBLS=('...
    num2str(aBLS) ')-fBLS=(1).mat']);
    APtimes = ll.APtimes(ll.APtimes>=USps&ll.APtimes<=(USps+USpd)); %(s)
    FR = 1./diff(APtimes);        % (Hz)
    plot(10^3*APtimes(1:end-1),FR,'color',redCmap(1+round((size(redCmap,1)-1)*(PaRRange(iPa)/PaRRange(end))),:));
    xlabel('Time [ms]'); ylabel('Firing rate [Hz]');
    rc = colorbar('location','eastoutside','colormap',redCmap);
    set(rc,'TickLabels',10^(-3)*PaRRange(end)*rc.Ticks);
    ylabel(rc,'Pressure [kPa]');
    hold off;
    else
    subplot(ceil(sqrt(length(PaRRange))),ceil(sqrt(length(PaRRange))),iPa);    
    %subplot(1,length(PaRRange),iPa);
    lr = load(['Chargevt(' MODELstr{Modelnr} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaRRange(iPa))) ')-ES(' num2str(ESps) ',' ...
    num2str(ESpd) ',' num2str(ESdc) ',' num2str(ESprf) ',' num2str(ESi) ')-aBLS=('...
    num2str(aBLS) ')-fBLS=(1).mat']);
    
    Qvt = lr.saveChargeSample(:,2); timeline = lr.saveChargeSample(:,1); 
    plot(timeline*10^3,Qvt);
    title({[num2str(PaRRange(iPa)*10^(-3)) 'kPa'],['(' num2str(Pa2I(PaRRange(iPa))) ' W/m^2)']},'Units','Normalized','Position',[1.05,0.5,1])
    title({[num2str(PaRRange(iPa)*10^(-3)) 'kPa'],['(' num2str(Pa2I(PaRRange(iPa))) ' W/m^2)']},'Units','Normalized')
    ylabel('Q [nC/cm^2]');
    if iPa == length(PaRRange)
    xlabel('Time [ms]');    
    end
    end
    if ~(keepTemps)
    delete(['APtimes(' MODELstr{Modelnr} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaRRange(iPa))) ')-ES(' num2str(ESps) ',' ...
    num2str(ESpd) ',' num2str(ESdc) ',' num2str(ESprf) ',' num2str(ESi) ')-aBLS=('...
    num2str(aBLS) ')-fBLS=(1).mat']);   
    delete(['Chargevt(' MODELstr{Modelnr} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaRRange(iPa))) ')-ES(' num2str(ESps) ',' ...
    num2str(ESpd) ',' num2str(ESdc) ',' num2str(ESprf) ',' num2str(ESi) ')-aBLS=('...
    num2str(aBLS) ')-fBLS=(1).mat']); 
    end
    end
end
end

if ~~(postProc) 
% Plot post-processed results:
% a. CW-stim
ESiRange = [0,-0.1,-0.2,-0.3,-0.4,-0.5];   % Simulated inhibitory perturbations
ESpdRange = 20e-3;                               % Pulse durations of inhibitory perturbations
USi2SP = [135,135,135,132,127,120];            % [W/m^2] Corresponding ultrasonic intensity threshold to obtain  a silenced plateau (from visual inspection)
FRwindow = 100e-3;             % Window in which the FR pre and post the GABAergic inhibition is calculated
DeltaFR = nan(length(PaRRange),length(ESiRange));

for jESI = 1:length(ESiRange) 
for iPa = 1:length(PaRRange)
if Pa2I(PaRRange(iPa)) < USi2SP(jESI)

ll = load(['APtimes(' MODELstr{Modelnr} ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaRRange(iPa))) ')-ES(' num2str(num2str(ESps)) ',' ...
num2str(ESpd) ',' num2str(ESdc) ',' num2str(ESprf) ',' num2str(ESiRange(jESI)) ')-aBLS=('...
num2str(aBLS) ')-fBLS=(1).mat']);

APtimes = ll.APtimes(ll.APtimes>=USps&ll.APtimes<=(USps+USpd)); %(s)
FRpre = sum(APtimes(APtimes>=(ESps-FRwindow)&(APtimes<=ESps)))/FRwindow;
FRpost = sum(APtimes(APtimes>=(ESps+ESpd)&(APtimes<=(ESps+ESpd+FRwindow))))/FRwindow;
DeltaFR(iPa,jESI) = FRpost-FRpre;
end
end
end

end