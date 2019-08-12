MODEL = 1;
switch MODEL
    case 1, modelStr = 'RS'; modelName = 'regular spiking';
    case 2, modelStr = 'FS'; modelName = 'fast spiking';
    case 3, modelStr = 'LTS'; modelName = 'low threshold spiking';
end
modelStr = ['SONIC-' modelStr '.mat'];
l = load(modelStr,'SONICtable');
SONICtable = l.SONICtable;

Q = SONICtable.QmRange;       % Order in the 4D matrix : (Q,Pa,Freq,aBLS)
Pa = SONICtable.USPaRange;
Freq = SONICtable.USfreqRange;
aBLS = SONICtable.aBLSRange; 

% Colours
load('redCmap.mat');
greenCmap = redCmap; greenCmap(:,[2 1]) = greenCmap(:,[1 2]);
blueCmap = redCmap; blueCmap(:,[3 1]) = blueCmap(:,[1 3]);

SONICfields = {'Zeff','Veff','Cmeff','a_m','apb_m','a_n','apb_n','a_h','apb_h','a_p','apb_p','ngend'};
plotFields = {'Zeff [nm]','Veff [mV]','Cmeff [\muF/cm^2]','a_m [1/ms]','b_m [1/ms]','a_n [1/ms]','b_n [1/ms]','a_h [1/ms]','b_h [1/ms]','a_p [1/ms]','b_p [1/ms]','ngend [mole]'};
UnitsScale = [10^9,1,100,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,1];

fig = figure('position',[1 41 1920 963]);
set(gcf,'color','w');
% 1. Dependence on (Q,USPa)   (Freq=500 kHz, a = 32 nm)
for j = 1:length(SONICfields)
red(j) = subplot(length(SONICfields),3,-2+3*j); %#ok<*SAGROW>
hold on;
for i = 1:length(Pa)
FieldName = SONICfields{j};
Y = SONICtable.(FieldName);
if strcmp(FieldName(1:3),'apb')
Y = Y-SONICtable.(['a_' FieldName(5)]);   
end
if i == 1               % No pressure
plot(10^5*Q,UnitsScale(j)*Y(:,i,3,1),'color','k','linestyle','--','linewidth',2);
else
plot(10^5*Q,UnitsScale(j)*Y(:,i,3,1),'color',redCmap(1+round((size(redCmap,1)-1)*(Pa(i)/Pa(end))),:));
end
end
ylabel(plotFields{j},'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
if j == length(SONICfields)
xlabel('Charge [nC/cm^2]');
else 
set(gca,'xtick',[]);   
end
hold off;
end
% 2. Dependence on (Q,Freq)  (Pa = 50 kPa, a = 32 nm)
colormap(greenCmap);
for j = 1:length(SONICfields)
green(j) = subplot(length(SONICfields),3,-1+3*j);
hold on;
for i = 1:length(Freq)
FieldName = SONICfields{j};
Y = SONICtable.(FieldName);
if strcmp(FieldName(1:3),'apb')
Y = Y-SONICtable.(['a_' FieldName(5)]);   
end
plot(10^5*Q,UnitsScale(j)*Y(:,37,i,1),'color',greenCmap(1+round((size(greenCmap,1)-1)*(Freq(i)/Freq(end))),:));
end
plot(10^5*Q,UnitsScale(j)*Y(:,1,3,1),'color','k','linestyle','--','linewidth',2);
if j == length(SONICfields)
xlabel('Charge [nC/cm^2]');
else 
set(gca,'xtick',[]);   
end
hold off;
end
% 3. Dependence on (Q,aBLS) (Freq = 500 kHz, Pa = 50 kPa); 
colormap(blueCmap);
for j = 1:length(SONICfields)
blue(j) = subplot(length(SONICfields),3,3*j);
hold on;
for i = 1:length(aBLS)
FieldName = SONICfields{j};
Y = SONICtable.(FieldName);
if strcmp(FieldName(1:3),'apb')
Y = Y-SONICtable.(['a_' FieldName(5)]);   
end
plot(10^5*Q,UnitsScale(j)*Y(:,37,3,i),'color',blueCmap(1+round((size(blueCmap,1)-1)*(aBLS(i)/aBLS(end))),:));
end
plot(10^5*Q,UnitsScale(j)*Y(:,1,3,1),'color','k','linestyle','--','linewidth',2);
if j == length(SONICfields)
xlabel('Charge [nC/cm^2]');
else 
set(gca,'xtick',[]);   
end
hold off;
end
% Colorbars
colormap(red(1),redCmap); colormap(blue(1),blueCmap); colormap(green(1),greenCmap);
rc = colorbar(red(1),'location','northoutside');
bc = colorbar(blue(1),'location','northoutside');
gc = colorbar(green(1),'location','northoutside');
set(rc,'TickLabels',10^(-3)*Pa(end)*rc.Ticks);
set(bc,'TickLabels',10^(9)*aBLS(end)*bc.Ticks);
set(gc,'TickLabels',10^(-3)*Freq(end)*gc.Ticks);
set(rc,'position',get(rc,'position')+[0 0.04 0 0]);
set(bc,'position',get(bc,'position')+[0 0.04 0 0]);
set(gc,'position',get(gc,'position')+[0 0.04 0 0]);
ylabel(rc,'Pressure [kPa]'); ylabel(gc,'Frequency [kHz]'); ylabel(bc,'BLS diameters [nm]'); 
axes = get(fig,'children'); tt = title(axes(25),['Effective parameters for ' modelName ' neuron']);
set(tt,'position',get(tt,'position')+[0 5 0],'fontsize',20);