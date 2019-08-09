MODEL = 1;
switch MODEL
    case 1, modelStr = 'RS';
    case 2, modelStr = 'FS';
    case 3, modelStr = 'LTS';     
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

SONICfields = {'Zeff','Veff','Cmeff','a_m','apb_m','a_n','apb_n','a_h','apb_h','a_p','apb_h','ngend'};
plotFields = {'Zeff [nm]','Veff [mV]','Cmeff [\muF/cm^2]','a_m [1/ms]','b_m [1/ms]','a_n [1/ms]','b_n [1/ms]','a_h [1/ms]','b_h [1/ms]','a_p [1/ms]','b_p [1/ms]','ngend [mole]'};
UnitsScale = [10^9,1,100,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,1];

figure;
% 1. Dependence on (Q,USPa)   (Freq=500 kHz, a = 32 nm)
redCmap(:,2) = 0; redCmap(:,3) = 0;
colormap(redCmap);
for j = 1:length(SONICfields)
subplot(length(SONICfields),3,-2+3*j);
hold on;
for i = 1:length(Pa)
FieldName = SONICfields{j};
Y = SONICtable.(FieldName);
if strcmp(FieldName(1:3),'apb')
Y = Y-SONICtable.(['a_' FieldName(5)]);   
end
plot(Q,UnitsScale(j)*Y(:,i,3,1),'color',redCmap(round(i*size(redCmap,1)/length(Pa)),:));
end
ylabel(plotFields{j});
hold off;
end

