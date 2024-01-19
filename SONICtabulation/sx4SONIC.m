% sx4SONIC (singleton expansion for SONIC) 
% Script will prepare SONICtable for cubic, spline, pchip or makima interpolation,
% by creating 4 points in the singleton directions

% ! As of now, the singleton dimensions are assumed to be USPa, USfreq,
% aBLS and fBLS
clear all; clc; 
tableVersion = '-v2'; 
modelName = 'RS';
NFS = 2;

%% - 

loadStr = ['SONIC-' modelName '-QoscFourier' num2str(NFS) '-FourierIn' num2str(NFS) tableVersion '.mat'];
saveStr = [loadStr(1:end-4) '-sx' '.mat'];

SONIC = load(loadStr); 
SONICtable = SONIC.SONICtable; clear SONIC;

SONICfields = fieldnames(SONICtable);
SONICrates = sort(SONICfields(cellfun(@(X) contains(X,'a_')|contains(X,'apb_'),SONICfields)));

QmRange = SONICtable.QmRange; 
USPaRange = SONICtable.USPaRange; USfreqRange = SONICtable.USfreqRange; 
aBLSRange = SONICtable.aBLSRange; fBLSRange = SONICtable.fBLSRange;
DeltaQmRange = SONICtable.DeltaQmRange; psiQRange = SONICtable.psiQRange;

USPaRange = linspace(USPaRange(1),USPaRange(end),4);
USfreqRange = linspace(USfreqRange(1),USfreqRange(end),4);
aBLSRange = linspace(aBLSRange(1),aBLSRange(end),4);
fBLSRange = linspace(fBLSRange(1),fBLSRange(end),4);

SONICtable.USPaRange = USPaRange; SONICtable.USfreqRange = USfreqRange;
SONICtable.aBLSRange = aBLSRange; SONICtable.fBLSRange = fBLSRange;

for k = 1:2
SONICtable.Veff = cat(2,SONICtable.Veff,SONICtable.Veff(:,1,:,:,:,:,:));
SONICtable.Veff = cat(3,SONICtable.Veff,SONICtable.Veff(:,:,1,:,:,:,:));
SONICtable.Veff = cat(4,SONICtable.Veff,SONICtable.Veff(:,:,:,1,:,:,:));
SONICtable.Veff = cat(5,SONICtable.Veff,SONICtable.Veff(:,:,:,:,1,:,:));

SONICtable.Zeff = cat(2,SONICtable.Zeff,SONICtable.Zeff(:,1,:,:,:,:,:));
SONICtable.Zeff = cat(3,SONICtable.Zeff,SONICtable.Zeff(:,:,1,:,:,:,:));
SONICtable.Zeff = cat(4,SONICtable.Zeff,SONICtable.Zeff(:,:,:,1,:,:,:));
SONICtable.Zeff = cat(5,SONICtable.Zeff,SONICtable.Zeff(:,:,:,:,1,:,:));

SONICtable.Cmeff = cat(2,SONICtable.Cmeff,SONICtable.Cmeff(:,1,:,:,:,:,:));
SONICtable.Cmeff = cat(3,SONICtable.Cmeff,SONICtable.Cmeff(:,:,1,:,:,:,:));
SONICtable.Cmeff = cat(4,SONICtable.Cmeff,SONICtable.Cmeff(:,:,:,1,:,:,:));
SONICtable.Cmeff = cat(5,SONICtable.Cmeff,SONICtable.Cmeff(:,:,:,:,1,:,:));

SONICtable.ngend = cat(2,SONICtable.ngend,SONICtable.ngend(:,1,:,:,:,:,:));
SONICtable.ngend = cat(3,SONICtable.ngend,SONICtable.ngend(:,:,1,:,:,:,:));
SONICtable.ngend = cat(4,SONICtable.ngend,SONICtable.ngend(:,:,:,1,:,:,:));
SONICtable.ngend = cat(5,SONICtable.ngend,SONICtable.ngend(:,:,:,:,1,:,:));

SONICtable.cfit = cat(2,SONICtable.cfit,SONICtable.cfit(:,1,:,:,:,:,:));
SONICtable.cfit = cat(3,SONICtable.cfit,SONICtable.cfit(:,:,1,:,:,:,:));
SONICtable.cfit = cat(4,SONICtable.cfit,SONICtable.cfit(:,:,:,1,:,:,:));
SONICtable.cfit = cat(5,SONICtable.cfit,SONICtable.cfit(:,:,:,:,1,:,:));

for i = 1:length(SONICrates)
SONICtable.(SONICrates{i}) = cat(2,SONICtable.(SONICrates{i}),SONICtable.(SONICrates{i})(:,1,:,:,:,:,:));
SONICtable.(SONICrates{i}) = cat(3,SONICtable.(SONICrates{i}),SONICtable.(SONICrates{i})(:,:,1,:,:,:,:));
SONICtable.(SONICrates{i}) = cat(4,SONICtable.(SONICrates{i}),SONICtable.(SONICrates{i})(:,:,:,1,:,:,:));
SONICtable.(SONICrates{i}) = cat(5,SONICtable.(SONICrates{i}),SONICtable.(SONICrates{i})(:,:,:,:,1,:,:));
end
end
save(saveStr,'SONICtable','-v7.3');