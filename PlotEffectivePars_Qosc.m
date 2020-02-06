MODEL = 1;
QoscStr = {'-QoscFourier2-FourierIn1','-QoscFourier8-FourierIn1'};

switch MODEL
    case 1, modelStrN = 'RS'; modelName = 'regular spiking';
    case 2, modelStrN = 'FS'; modelName = 'fast spiking';
    case 3, modelStrN = 'LTS'; modelName = 'low threshold spiking';
end
for iload = 1:numel(QoscStr)
modelStr = ['SONIC-' modelStrN QoscStr{iload} '.mat'];
l(iload) = load(modelStr,'SONICtable');
end
figure('units','normalized','position',[0 0 1 1]); set(gcf,'color','w');
for ii = 1:numel(QoscStr)
SONICtable = l(ii).SONICtable;

Q = SONICtable.QmRange;       % Order in the 4D matrix : (Q,Pa,Freq,aBLS)
Pa = SONICtable.USPaRange;
Freq = SONICtable.USfreqRange;
aBLS = SONICtable.aBLSRange; 
fBLS = SONICtable.fBLSRange;
DeltaQm = SONICtable.DeltaQmRange;
psiQRange = SONICtable.psiQRange;

% gof values
rsquare = cellfun(@(X) X(2),SONICtable.gof);
adjrsquare = cellfun(@(X) X(4),SONICtable.gof);
rsquare = permute(rsquare(:,:,1,1,1,:,1),[1 2 6 3 4 5]);
adjrsquare = permute(adjrsquare(:,:,1,1,1,:,1),[1 2 6 3 4 5]);
adjrsquareStr.(char('a' + ii-1)) = adjrsquare;
rsquare(:,1,1) = nan; adjrsquare(:,1,1) = nan;  % Constant fit, rsquare values are determined by numerical imprecisions while the sse\approx 0

iPa = 51; iQm = 20; iDeltaQm = 1;
colormap(cm_plasma(100));
subplot(3,3,3+1+3*(ii-1)); hold on; surf(10^(-3)*Pa,10^5*DeltaQm,permute(adjrsquare(iQm,:,:),[3 2 1]),'edgecolor','none','facecolor','interp');
if (ii==1), title(['Q_m = ' num2str(10^(5)*Q(iQm)) ' nC/cm^2']); end
plot3(600,500,adjrsquare(iQm,end,end),'marker','x','color','b','markerfacecolor','b','markersize',10,'linewidth',2);
plot3(600,0,adjrsquare(iQm,end,1),'marker','x','color','r','markerfacecolor','r','markersize',10,'linewidth',2);
plot3(31,0,adjrsquare(iQm,31,1),'marker','o','color','b','markerfacecolor','b','markersize',10,'linewidth',2);
plot3(0,400,adjrsquare(iQm,1,17),'marker','o','color','r','markerfacecolor','r','markersize',10,'linewidth',2);
ylabel('\DeltaQ_m [nC/cm^2]'); xlabel('P_{US} [kPa]'); zlabel('\boldmath$\overline{R}^2$ [-]','interpreter','latex');
view([106.5266 48.3254]); caxis manual; caxis([min(adjrsquare(:)),max(adjrsquare(:))]);
xlim(10^(-3)*[Pa(1),Pa(end)]); ylim(10^5*[DeltaQm(1),DeltaQm(end)]); hold off;
subplot(3,3,3+2+3*(ii-1));
hold on;
surf(10^(-3)*Pa,10^5*Q,adjrsquare(:,:,iDeltaQm),'edgecolor','none','facecolor','interp');
plot3(600,-77.9,adjrsquare(iQm,end,1),'marker','x','color','r','markerfacecolor','r','markersize',10,'linewidth',2);
plot3(31,-77.9,adjrsquare(iQm,31,1),'marker','o','color','b','markerfacecolor','b','markersize',10,'linewidth',2);
if (ii == 1), title(['\DeltaQ_m = ' num2str(10^(5)*DeltaQm(iDeltaQm)) ' nC/cm^2']); end
ylabel('Q_m [nC/cm^2]'); xlabel('P_{US} [kPa]');
view([106.5266 48.3254]);  caxis manual; caxis([min(adjrsquare(:)),max(adjrsquare(:))]);
xlim(10^(-3)*[Pa(1),Pa(end)]); ylim(10^5*[Q(1),Q(end)]);
hold off;
subplot(3,3,3+3+3*(ii-1)); 
hold on; surf(10^5*Q,10^5*DeltaQm,permute(adjrsquare(:,iPa,:),[3 1 2]),'edgecolor','none','facecolor','interp');
plot3(-77.9,500,adjrsquare(iQm,end,end),'marker','x','color','b','markerfacecolor','b','markersize',10,'linewidth',2);
plot3(-77.9,0,adjrsquare(iQm,end,1),'marker','x','color','r','markerfacecolor','r','markersize',10,'linewidth',2);

if (ii==1), title(['P_{US} = ' num2str(10^(-3)*Pa(iPa)) ' kPa']); end
ylabel('\DeltaQ_m [nC/cm^2]'); xlabel('Q_m [nC/cm^2]');
view([106.5266-90 48.3254]); colorbar;  caxis manual; caxis([min(adjrsquare(:)),max(adjrsquare(:))]);
xlim(10^5*[Q(1),Q(end)]); ylim(10^5*[DeltaQm(1),DeltaQm(end)]);
hold off;
end
subplot(3,3,3); hold on; histogram(adjrsquareStr.('a'),'binedges',(0:0.001:1),'Normalization','probability'); histogram(adjrsquareStr.('b'),'BinEdges',(0:0.001:1),'Normalization','Probability'); 
xlabel('\boldmath$\overline{R}^2$','interpreter','latex'); ylabel('Incidence [-]'); xlim([0.95,1]); ylim([0,0.2]);
ll = legend({'N_{FS}=2','N_{FS}=8'},'position',[0.6972 0.8284 0.0594 0.0778]); hold off; 

tempfitEx(1) = load('fitEx-RS-(Q,P,freq,a,fBLS,corrPEC,deltaQ,deltaPsi)-(-77.9e-5,600e3,500e3,32e-9,1,0,500e-5,0).mat');
tempfitEx(2) = load('fitEx-RS-(Q,P,freq,a,fBLS,corrPEC,deltaQ,deltaPsi)-(-77.9e-5,600e3,500e3,32e-9,1,0,0,0).mat');
tempfitEx(3) = load('fitEx-RS-(Q,P,freq,a,fBLS,corrPEC,deltaQ,deltaPsi)-(-77.9e-5,17.2e3,500e3,32e-9,1,0,0,0).mat');
tempfitEx(4) = load('fitEx-RS-(Q,P,freq,a,fBLS,corrPEC,deltaQ,deltaPsi)-(-77.9e-5,0,500e3,32e-9,1,0,400e-5,0).mat');

for i = 1:4
fitEx(i) = tempfitEx(i).fitExample;
end

refitN = @(X,t,Nfr) X(1)+X(2:2:end-1)*cos((1:Nfr)'*X(end)*t)+X(3:2:end-1)*sin((1:Nfr)'*X(end)*t);

subplot(3,3,1);
hold on; iff = 1;
leg(1) = plot(nan,'color','k','linewidth',2,'linestyle','--','marker','none');
leg(2) = plot(nan,'color','k','linewidth',2,'linestyle','-.','marker','none');
yyaxis left, plot(fitEx(iff).tPeriod-fitEx(iff).tPeriod(1),fitEx(iff).Vperiod,'color','b','linewidth',1.5);
cf2 = plot(fitEx(iff).tPeriod-fitEx(iff).tPeriod(1),refitN(coeffvalues(fitEx(iff).cfit2),fitEx(iff).tPeriod,2),'b--'); set(cf2,'linewidth',2);
cf8 = plot(fitEx(iff).tPeriod-fitEx(iff).tPeriod(1),refitN(coeffvalues(fitEx(iff).cfit8),fitEx(iff).tPeriod,8),'b-.'); set(cf8,'linewidth',2);
set(gca,'YColor','b')
box off;
ylabel('V [mV]'); 
hold off; 
iff=2;
yyaxis right, 
hold on; plot(fitEx(iff).tPeriod-fitEx(iff).tPeriod(1),fitEx(iff).Vperiod,'color','r','linewidth',1.5);
cf2 = plot(fitEx(iff).tPeriod-fitEx(iff).tPeriod(1),refitN(coeffvalues(fitEx(iff).cfit2),fitEx(iff).tPeriod,2),'r--'); set(cf2,'linewidth',2);
cf8 = plot(fitEx(iff).tPeriod-fitEx(iff).tPeriod(1),refitN(coeffvalues(fitEx(iff).cfit8),fitEx(iff).tPeriod,8),'r-.'); set(cf8,'linewidth',2);
set(gca,'YColor','r')
%ylabel('V [mV]');
hold off;
xlabel('time [\mus]');
legend(leg,{'N_{FS}=2';'N_{FS}=8'},'location','westoutside'); xlim([0,max(fitEx(1).tPeriod(end)-fitEx(1).tPeriod(1),fitEx(2).tPeriod(end)-fitEx(2).tPeriod(1))]);
set(gca,'xticklabels',num2cell(10^6*get(gca,'xtick')));

subplot(3,3,2);
hold on; iff = 3;
yyaxis left, plot(fitEx(iff).tPeriod-fitEx(iff).tPeriod(1),fitEx(iff).Vperiod,'color','b','linewidth',1.5);
cf2 = plot(fitEx(iff).tPeriod-fitEx(iff).tPeriod(1),refitN(coeffvalues(fitEx(iff).cfit2),fitEx(iff).tPeriod,2),'b--'); set(cf2,'linewidth',2);
cf8 = plot(fitEx(iff).tPeriod-fitEx(iff).tPeriod(1),refitN(coeffvalues(fitEx(iff).cfit8),fitEx(iff).tPeriod,8),'b-.'); set(cf8,'linewidth',2);
%ylabel('V [mV]'); 
set(gca,'YColor','b')
hold off; 
iff=4;
yyaxis right, 
hold on; plot(fitEx(iff).tPeriod-fitEx(iff).tPeriod(1),fitEx(iff).Vperiod,'color','r','linewidth',1.5);
cf2 = plot(fitEx(iff).tPeriod-fitEx(iff).tPeriod(1),refitN(coeffvalues(fitEx(iff).cfit2),fitEx(iff).tPeriod,2),'r--'); set(cf2,'linewidth',2);
cf8 = plot(fitEx(iff).tPeriod-fitEx(iff).tPeriod(1),refitN(coeffvalues(fitEx(iff).cfit8),fitEx(iff).tPeriod,8),'r-.'); set(cf8,'linewidth',2);
ylabel('V [mV]');
hold off;
xlabel('time [\mus]');
legend('off'); xlim([0,max(fitEx(3).tPeriod(end)-fitEx(3).tPeriod(1),fitEx(4).tPeriod(end)-fitEx(4).tPeriod(1))]);
set(gca,'xticklabels',num2cell(10^6*get(gca,'xtick'))); set(gca,'YColor','r')
set(findobj('type','axes'),'fontweight','bold','fontsize',15);


% cfit values
FourierComps = cellfun(@(X) X(1,2:end-1),SONICtable.cfit,'UniformOutput',0); % Extract fourier components
ConfI = cellfun(@(X) X(2:3,2:end-1),SONICtable.cfit,'UniformOutput',0); % Extract confidence intervals

for i = 1:numel(FourierComps{1})/2
    SONICtable.(['ampV_' num2str(i)]) = cellfun(@(X) sqrt(X(2*i-1)^2+X(2*i)^2),FourierComps);
    SONICtable.(['DeltaPhi_' num2str(i)]) = cellfun(@(X) -atan2(X(2*i),X(2*i-1)),FourierComps);
end

% Colours
load('redCmap.mat');
greenCmap = redCmap; greenCmap(:,[2 1]) = greenCmap(:,[1 2]);
blueCmap = redCmap; blueCmap(:,[3 1]) = blueCmap(:,[1 3]);

SONICfields = fieldnames(SONICtable);
for i = length(SONICfields):-1:1
if contains(SONICfields{i},'Range') || contains(SONICfields{i},'cfit') || contains(SONICfields{i},'gof') || strcmp(SONICfields{i},'DeltaPhi')
SONICfields(i) = [];
end
end
for i = length(SONICfields):-1:1
if ~contains(SONICfields{i},'amp') 
SONICfields(i) = [];
end
end

plotFields = cell(1,length(SONICfields)); UnitsScale = zeros(1,length(SONICfields));
for i = 1:length(SONICfields)
   switch SONICfields{i}
       case 'Zeff', plotFields{i} = 'Zeff [nm]'; UnitsScale(i) = 10^9;
       case 'Veff', plotFields{i} = 'Veff [mV]'; UnitsScale(i) = 1;
       case 'Cmeff', plotFields{i} = 'Cmeff [\muF/cm^2]'; UnitsScale(i) = 100;
       case 'ngend', plotFields{i} = 'ngend [mole]';UnitsScale(i) = 1;
       case 'ampV', plotFields{i} = 'V_{osc} [mV]'; UnitsScale(i) = 1;
       %case 'DeltaPhi', plotFields{i} = '\Delta \phi [rad]'; UnitsScale(i) = 1;
       otherwise
           if contains(SONICfields{i},'a_')
               plotFields{i} = [SONICfields{i} '[1/ms]']; UnitsScale(i) = 0.001;
           elseif contains(SONICfields{i},'apb_')
               plotFields{i} =['b_' SONICfields{i}(5:end) '[1/ms]']; UnitsScale(i) = 0.001;
           elseif contains(SONICfields{i},'ampV')
               plotFields{i} = [SONICfields{i} ' [mV]']; UnitsScale(i) = 1;
           elseif contains(SONICfields{i},'DeltaPhi')
               plotFields{i} = [SONICfields{i} ' [rad]']; UnitsScale(i) = 1;
           end
   end   
end   

%-- Surf plots of effective parameters

fs = 17;
fig = figure('position',[1 41 1920 963]);
set(gcf,'color','w');
colormap('jet');
% 1. Dependence on (Q,USPa)   (Freq=500 kHz, a = 32 nm, fBLS = 1, DeltaQm = 0)
for j = 1:length(SONICfields)
surfQP(j) = subplot(length(SONICfields),3,-2+3*j); %#ok<*SAGROW>
hold on;
FieldName = SONICfields{j};
Y = SONICtable.(FieldName);
if strcmp(FieldName(1:3),'apb')
Y = Y-SONICtable.(['a_' FieldName(5)]);   
end
maxY = max(max(Y(:,:,1,1,1,1)));
if contains(FieldName,'DeltaPhi')
surf(10^(-3)*Pa,10^5*Q,UnitsScale(j)*(mod(Y(:,:,1,1,1,1,1),2*pi)),'edgecolor','none','facecolor','interp');
else
surf(10^(-3)*Pa,10^5*Q,UnitsScale(j)*Y(:,:,1,1,1,1,1)/maxY,'edgecolor','none','facecolor','interp');
end
xlim([10^(-3)*Pa(1),10^(-3)*Pa(end)]); ylim(10^5*[Q(1),Q(end)]);
if (j==1), ylabel('Q [nC/cm^2]'); set(get(gca,'YLabel'),'Rotation',0,'Position',[-37.9462 88.4078 0]); end 
%colorbar;
if strcmp(plotFields{j}(1:3),'amp')
text(-0.3345,0.7937,plotFields{j}(4:end),'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontweight','bold','fontsize',fs);
text(0.9977,1.4594,['(' num2str(maxY,3) ')'],'Units', 'Normalized', 'VerticalAlignment', 'Top','horizontalalignment','right','Fontweight','bold','fontsize',fs);
else
text(-0.3345,0.7937,plotFields{j},'Units', 'Normalized', 'VerticalAlignment', 'Top','Fontweight','bold','fontsize',fs);
text(0.9977,1.4594,['(' num2str(maxY,3) ')'],'Units', 'Normalized', 'VerticalAlignment', 'Top','horizontalalignment','right','Fontweight','bold','fontsize',fs);
end
if j == length(SONICfields)
xlabel('P_{US} [kPa]');
else 
set(gca,'xtick',[]);   
end
yt = get(gca,'ytick');
set(gca,'ytick',[yt(1),yt(end)]);
hold off;
end

% 2. Dependence on (Q,DeltaQm)  (freq = 500 kHz, Pa = 0 kPa, a = 32 nm, fBLS = 1)
for j = 1:length(SONICfields)
surfQDeltaQm(j) = subplot(length(SONICfields),3,-1+3*j);
hold on;
FieldName = SONICfields{j};
Y = SONICtable.(FieldName);
if strcmp(FieldName(1:3),'apb')
Y = Y-SONICtable.(['a_' FieldName(5)]);   
end
maxY = max(max(Y(:,1,1,1,1,:)));
if contains(FieldName,'DeltaPhi')
surf(10^5*DeltaQm,10^5*Q,permute(UnitsScale(j)*mod(Y(:,1,1,1,1,:,1),2*pi),[1 6 2 3 4 5]),'edgecolor','none','facecolor','interp');
else
surf(10^5*DeltaQm,10^5*Q,permute(UnitsScale(j)*Y(:,1,1,1,1,:,1)./maxY,[1 6 2 3 4 5]),'edgecolor','none','facecolor','interp');
end
xlim([10^(5)*DeltaQm(1),10^(5)*DeltaQm(end)]); ylim(10^5*[Q(1),Q(end)]);
if (j==1), ylabel('Q [nC/cm^2]'); set(get(gca,'YLabel'),'Rotation',0,'Position',[-37.9462 88.4078 0]); end 
%colorbar;
if strcmp(plotFields{j}(1:3),'amp')
%title([plotFields{j}(4:end) ' (' num2str(maxY,3) ')']);
text(0.9977,1.4594,['(' num2str(maxY,3) ')'],'Units', 'Normalized', 'VerticalAlignment', 'Top','horizontalalignment','right','Fontweight','bold','fontsize',fs);
else
%title([plotFields{j} ' (' num2str(maxY,3) ')']);
text(0.9977,1.4594,['(' num2str(maxY,3) ')'],'Units', 'Normalized', 'VerticalAlignment', 'Top','horizontalalignment','right','Fontweight','bold','fontsize',fs);
end
if j == length(SONICfields)
xlabel('\Delta Q [nC/cm^2]');
else 
set(gca,'xtick',[]);   
end
yt = get(gca,'ytick');
set(gca,'ytick',[yt(1),yt(end)]);
hold off;
end
% 2. Dependence on (USPa,DeltaQm)  (freq = 500 kHz, Q = 0 nC/cm^2, a = 32 nm, fBLS = 1)
for j = 1:length(SONICfields)
surfPDeltaQm(j) = subplot(length(SONICfields),3,3*j);
hold on;
FieldName = SONICfields{j};
Y = SONICtable.(FieldName);
if strcmp(FieldName(1:3),'apb')
Y = Y-SONICtable.(['a_' FieldName(5)]);   
end
maxY = max(max(Y(1,:,1,1,1,:,1)));
if contains(FieldName,'DeltaPhi')
surf(10^5*DeltaQm,10^(-3)*Pa,permute(UnitsScale(j)*mod(Y(1,:,1,1,1,:,1),2*pi),[2 6 1 3 4 5 7]),'edgecolor','none','facecolor','interp');
else
surf(10^5*DeltaQm,10^(-3)*Pa,permute(UnitsScale(j)*Y(1,:,1,1,1,:,1)./maxY,[2 6 1 3 4 5 7]),'edgecolor','none','facecolor','interp');
end
ylim([10^(-3)*Pa(1),10^(-3)*Pa(end)]); xlim(10^5*[DeltaQm(1),DeltaQm(end)]);
if (j==1), ylabel('P_{US} [kPa]'); set(get(gca,'YLabel'),'Rotation',0,'Position',[-22.0925 596.1001 7.1054e-15]); end 
%colorbar;
if strcmp(plotFields{j}(1:3),'amp')
%title([plotFields{j}(4:end) ' (' num2str(maxY,3) ')']);
%title(plotFields{j}(4:end));
text(0.9977,1.4594,['(' num2str(maxY,3) ')'],'Units', 'Normalized', 'VerticalAlignment', 'Top','HorizontalAlignment','right','Fontweight','bold','fontsize',fs);
else
%title([plotFields{j} ' (' num2str(maxY,4) ')']);
text(0.9977,1.4594,['(' num2str(maxY,4) ')'],'Units', 'Normalized', 'VerticalAlignment', 'Top','HorizontalAlignment','right','Fontweight','bold','fontsize',fs);
end
if j == length(SONICfields)
xlabel('\Delta Q_m [nC/cm^2]'); 
else 
set(gca,'xtick',[]);   
end
yt = get(gca,'ytick');
set(gca,'ytick',[yt(1),yt(end)]);
hold off;
end
set(findobj('type','axes'),'fontsize',fs,'fontweight','bold');
cb = colorbar(surfQDeltaQm(1),'location','northoutside');

% ---- Effective params plot
SONICfields = fieldnames(SONICtable);
for i = length(SONICfields):-1:1
if contains(SONICfields{i},'Range') || contains(SONICfields{i},'cfit') || contains(SONICfields{i},'gof') || strcmp(SONICfields{i},'DeltaPhi')
SONICfields(i) = [];
end
end
for i = length(SONICfields):-1:1
if contains(SONICfields{i},'amp') || contains(SONICfields{i},'DeltaPhi')
SONICfields(i) = [];
end
end

fig = figure('position',[1 41 1920 963]);
set(gcf,'color','w');
% 1. Dependence on (Q,USPa)   (Freq=500 kHz, a = 32 nm, fBLS = 1, DeltaQm = 0)
for j = 1:length(SONICfields)
red(j) = subplot(length(SONICfields),3,-2+3*j); %#ok<*SAGROW>
hold on;
for i = 1:length(Pa)
FieldName = SONICfields{j};
Y = SONICtable.(FieldName);
if strcmp(FieldName(1:3),'apb')
Y = Y-SONICtable.(['a_' FieldName(5)]);   
end
if contains(FieldName,'DeltaPhi')
if i == 1               % No pressure
plot(10^5*Q,UnitsScale(j)*unwrap(Y(:,i,1,1,1,1)),'color','k','linestyle','--','linewidth',2);
else
plot(10^5*Q,UnitsScale(j)*unwrap(Y(:,i,1,1,1,1)),'color',redCmap(1+round((size(redCmap,1)-1)*(Pa(i)/Pa(end))),:));
end
else
if i == 1               % No pressure
plot(10^5*Q,UnitsScale(j)*Y(:,i,1,1,1,1),'color','k','linestyle','--','linewidth',2);
else
plot(10^5*Q,UnitsScale(j)*Y(:,i,1,1,1,1),'color',redCmap(1+round((size(redCmap,1)-1)*(Pa(i)/Pa(end))),:));
end
end
end
ylabel(plotFields{j},'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
if j == length(SONICfields)
xlabel('Q [nC/cm^2]');
else 
set(gca,'xtick',[]);   
end
hold off;
end

% 2. Dependence on (Q,DeltaQm)  (freq = 500 kHz, Pa = 0 kPa, a = 32 nm, fBLS = 1)
colormap(greenCmap);
for j = 1:length(SONICfields)
green(j) = subplot(length(SONICfields),3,-1+3*j);
hold on;
for i = 1:length(DeltaQm)
FieldName = SONICfields{j};
Y = SONICtable.(FieldName);
if strcmp(FieldName(1:3),'apb')
Y = Y-SONICtable.(['a_' FieldName(5)]);   
end
if contains(FieldName,'DeltaPhi')
plot(10^5*Q,UnitsScale(j)*unwrap(Y(:,1,1,1,1,i)),'color',greenCmap(1+round((size(greenCmap,1)-1)*(DeltaQm(i)/DeltaQm(end))),:));
else
plot(10^5*Q,UnitsScale(j)*Y(:,1,1,1,1,i),'color',greenCmap(1+round((size(greenCmap,1)-1)*(DeltaQm(i)/DeltaQm(end))),:));
end
end
if contains(FieldName,'DeltaPhi')
plot(10^5*Q,UnitsScale(j)*unwrap(Y(:,1,1,1,1,1)),'color','k','linestyle','--','linewidth',2);
else
plot(10^5*Q,UnitsScale(j)*Y(:,1,1,1,1,1),'color','k','linestyle','--','linewidth',2);
end
if j == length(SONICfields)
xlabel('Q [nC/cm^2]');
else 
set(gca,'xtick',[]);   
end
hold off;
end
% 2. Dependence on (Q,DeltaQm)  (freq = 500 kHz, Pa = 50 kPa, a = 32 nm, fBLS = 1)
colormap(blueCmap);
for j = 1:length(SONICfields)
blue(j) = subplot(length(SONICfields),3,3*j);
hold on;
for i = 1:length(DeltaQm)
FieldName = SONICfields{j};
Y = SONICtable.(FieldName);
if strcmp(FieldName(1:3),'apb')
Y = Y-SONICtable.(['a_' FieldName(5)]);   
end
plot(10^5*Q,UnitsScale(j)*Y(:,37,1,1,1,i),'color',blueCmap(1+round((size(blueCmap,1)-1)*(DeltaQm(i)/DeltaQm(end))),:));
end
plot(10^5*Q,UnitsScale(j)*Y(:,1,1,1,1,1),'color','k','linestyle','--','linewidth',2);
if j == length(SONICfields)
xlabel('Q [nC/cm^2]');
else 
set(gca,'xtick',[]);   
end
hold off;
end
% Colorbars
colormap(red(1),redCmap); colormap(green(1),greenCmap); colormap(blue(1),blueCmap);
rc = colorbar(red(1),'location','northoutside');
bc = colorbar(blue(1),'location','northoutside');
gc = colorbar(green(1),'location','northoutside');
set(rc,'TickLabels',10^(-3)*Pa(end)*rc.Ticks);
set(bc,'TickLabels',10^(5)*DeltaQm(end)*bc.Ticks);
set(gc,'TickLabels',10^(5)*DeltaQm(end)*gc.Ticks);
set(rc,'position',get(rc,'position')+[0 0.04 0 0]);
set(bc,'position',get(bc,'position')+[0 0.04 0 0]);
set(gc,'position',get(gc,'position')+[0 0.04 0 0]);
ylabel(rc,'P_{US} (\DeltaQ_m=0 nC/cm^2) [kPa]'); ylabel(gc,'\DeltaQ_m (P_{US} = 0 Pa) [nC/cm^2]'); ylabel(bc,'\DeltaQ_m (P_{US} = 50 kPa) [nC/cm^2]'); 
axes = get(fig,'children'); tt = title(axes(25),['Effective parameters for ' modelName ' neuron']);
set(tt,'position',get(tt,'position')+[0 5 0],'fontsize',20);