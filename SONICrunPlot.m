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

figure;
% Fig.1(a-c) (top)
Tsim = 0.3;     % (s)
USpd = 0.150;   % (s)
USps = 0.05;  % (s)
USdc = 1; USprf = 0;   % (-) , (Hz)

USfreq = {(500e3),[20e3,4e6],(500e3)}; % (Hz)
aBLS = {(32e-9),(32e-9),[16e-9,64e-9]};  % (m)
RelPa = {[-5e3,0,20e3],(20e3),(20e3)};  % relative to threshold (Pa)


% (A) Determine threshold
cd(SONICpath);

for i = 1:3
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
    num2str(USfreq{i}(iFreq)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(0) ')-ES(0,0,1,0,0)-aBLS(' ...
    num2str(aBLS{i}(iaBLS)) ').mat']);
end
end
end

for i = 1:length(PaThRange)
SONICrun(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),...
num2str(Pa2I(PaThRange(i))),'0','0','1','0','0','0',num2str(MODEL),'0','0','0');
end
figure;
subplot(7,3,1);
hold on;
for i = 1:length(PaThRange)
lt2 = load(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaThRange(i))) ')-ES(0,0,1,0,0).mat']);
Qvt = lt2.saveChargeSample(:,2); timeline = lt2.saveChargeSample(:,1); 
plot(10^3*timeline,Qvt);
delete(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaThRange(i))) ')-ES(0,0,1,0,0).mat']);
end
cd(NICEpath);
for i = 1:length(PaThRange)
funPES(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),...
num2str(Pa2I(PaThRange(i))),'0','0','1','0','0','0',num2str(MODEL),'0','0','0');
end
for i = 1:length(PaThRange)
lt3 = load(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaThRange(i))) ')-ES(0,0,1,0,0).mat']);
Qvt = lt3.saveChargeSample(:,2); timeline = lt3.saveChargeSample(:,1); 
plot(10^3*timeline,Qvt);
delete(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaThRange(i))) ')-ES(0,0,1,0,0).mat']);
end
xlabel('Time [ms]');
ylabel('Charge [nC/cm^2]');
set(gcf,'color','w');
hold off;

% (B) Sweep over pressure
Latency = cell(2,1); FR = cell(2,1); spA = cell(2,1);
for pathNum = 1:2
if pathNum == 1
    cd(SONICpath);
    for i = 1:length(PaR)
    SONICrun(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),...
    num2str(Pa2I(PaR(i))),'0','0','1','0','0','0',num2str(MODEL),'0','0','0');
    end
else
    cd(NICEpath);
    for i = 1:length(PaR)
    funPES(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq),num2str(USdc),num2str(USprf),...
    num2str(Pa2I(PaR(i))),'0','0','1','0','0','0',num2str(MODEL),'0','0','0');
    end
end

Latency{pathNum,1} = zeros(length(PaR),1);  % Latency (ms) 
FR{pathNum,1} = zeros(length(PaR),1);          % Firing rate (Hz)
spA{pathNum,1} = zeros(length(PaR),1);      % spiking amplitude (nC/cm^2)
for i = 1:length(PaR)
ll = load(['APtimes(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaR(i))) ')-ES(0,0,1,0,0).mat']);
lr = load(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaR(i))) ')-ES(0,0,1,0,0).mat']);
APtimes = ll.APtimes(ll.APtimes>=USps&ll.APtimes<=(USps+USpd)); %(s)
Latency{pathNum,1}(i) = 10^3*(APtimes(1)-USps);    % (ms)
FR{pathNum,1}(i) = length(APtimes)./(USpd);        % (Hz)
Qvt = lr.saveChargeSample(:,2); timeline = lr.saveChargeSample(:,1);  % (nC/cm^2) and (s)
QvtCrop = Qvt(timeline>=APtimes(2)&timeline<=APtimes(end-1));  % Crop outer lying action potentials to determine amplitude
timelineCrop = timeline(timeline>=APtimes(2)&timeline<=APtimes(end-1));
spA{pathNum,1}(i) = (max(QvtCrop)-min(QvtCrop));       % (nC/cm^2)
end
for i = 1:length(PaR)
    delete(['APtimes(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaR(i))) ')-ES(0,0,1,0,0).mat']);
    delete(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(PaR(i))) ')-ES(0,0,1,0,0).mat']);
end
end

subplot(7,3,2);
hold on;
plot(1e-3*PaR,Latency{2,1},'linestyle','-','color',[0.5 0.5 0.5],'linewidth',2);
plot(1e-3*PaR,Latency{1,1},'linestyle','--','color','k','marker','o');
ylabel('Latency [ms]');
hold off;
subplot(7,3,3);
hold on;
plot(1e-3*PaR,FR{2,1},'linestyle','-','color',[0.5 0.5 0.5],'linewidth',2);
plot(1e-3*PaR,FR{1,1},'linestyle','--','color','k','marker','o');
ylabel('Firing rate [Hz]');
hold off;
subplot(7,3,4);
hold on;
plot(1e-3*PaR,spA{2,1},'linestyle','-','color',[0.5 0.5 0.5],'linewidth',2);
plot(1e-3*PaR,spA{1,1},'linestyle','--','color','k','marker','o');
ylabel('Spike amp. [nC/cm^2]');
xlabel('Amplitude [kPa]');
hold off;
set(gcf,'color','white');
set(findobj('type','axes'),'fontsize',18);
set(findobj('type','axes'),'xlim',1e-3*[PaR(1),PaR(end)]);
set(findobj('type','axes'),'box','off');

% Fig. 1(b)
USfreq = logspace(log10(20e3),log10(4e6),7);        % (Hz)
ThreshPavFreq = zeros(1,7);             % (Pa)

cd(SONICpath);
for i = 1:length(USfreq)
SONICrun(num2str(Tsim),'1',num2str(USps),num2str(USpd),num2str(USfreq(i)),num2str(USdc),num2str(USprf),...
'0','0','0','1','0','0','0',num2str(MODEL),'0','0','0');
ltf = load(['Thresh(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq(i)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(0) ')-ES(0,0,1,0,0).mat']);
ThreshPavFreq(i) = I2Pa(ltf.IIpa);
delete(['Thresh(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq(i)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(0) ')-ES(0,0,1,0,0).mat']);
end

for i = 1:length(ThreshPavFreq)
SONICrun(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq(1)),num2str(USdc),num2str(USprf),...
num2str(Pa2I(ThreshPavFreq(i))),'0','0','1','0','0','0',num2str(MODEL),'0','0','0');
SONICrun(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq(end)),num2str(USdc),num2str(USprf),...
num2str(Pa2I(ThreshPavFreq(i))),'0','0','1','0','0','0',num2str(MODEL),'0','0','0');
end
figure;
subplot(7,3,7);
hold on;
for i = 1:length(ThreshPavFreq)
ltf2 = load(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq(1)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(ThreshPavFreq(i))) ')-ES(0,0,1,0,0).mat']);
ltf3 = load(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq(end)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(ThreshPavFreq(i))) ')-ES(0,0,1,0,0).mat']);
Qvt2 = lt2.saveChargeSample(:,2); timeline2 = lt2.saveChargeSample(:,1);
Qvt3 = lt3.saveChargeSample(:,2); timeline3 = lt3.saveChargeSample(:,1); 
plot(10^3*timeline2,Qvt2); plot(10^3*timeline3,Qvt3);
delete(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq(1)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(ThreshPavFreq(i))) ')-ES(0,0,1,0,0).mat']);
delete(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq(end)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(ThreshPavFreq(i))) ')-ES(0,0,1,0,0).mat']);
end
cd(NICEpath);
for i = 1:length(ThreshPavFreq)
funPES(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq(1)),num2str(USdc),num2str(USprf),...
num2str(Pa2I(ThreshPavFreq(i))),'0','0','1','0','0','0',num2str(MODEL),'0','0','0');
funPES(num2str(Tsim),'2',num2str(USps),num2str(USpd),num2str(USfreq(end)),num2str(USdc),num2str(USprf),...
num2str(Pa2I(ThreshPavFreq(i))),'0','0','1','0','0','0',num2str(MODEL),'0','0','0');
end
for i = 1:length(ThreshPavFreq)
lt4 = load(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq(1)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(ThreshPavFreq(i))) ')-ES(0,0,1,0,0).mat']);
lt5 = load(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq(end)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(ThreshPavFreq(i))) ')-ES(0,0,1,0,0).mat']);
Qvt4 = lt4.saveChargeSample(:,2); timeline4 = lt4.saveChargeSample(:,1); 
Qvt5 = lt5.saveChargeSample(:,2); timeline5 = lt5.saveChargeSample(:,1); 
plot(10^3*timeline4,Qvt4); plot(10^3*timeline5,Qvt5);
delete(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq(1)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(ThreshPavFreq(i))) ')-ES(0,0,1,0,0).mat']);
delete(['Chargevt(' ModelStr ')-Tsim=' num2str(Tsim) '-US(' num2str(USps) ',' num2str(USpd) ',' ...
    num2str(USfreq(end)) ',' num2str(USdc) ',' num2str(USprf) ',' num2str(Pa2I(ThreshPavFreq(i))) ')-ES(0,0,1,0,0).mat']);
end
xlabel('Time [ms]');
ylabel('Charge [nC/cm^2]');
set(gcf,'color','w');
hold off;

PaR = 1e3*logspace(log10(50),log10(150),4);   % Pressure range (Pa) 



