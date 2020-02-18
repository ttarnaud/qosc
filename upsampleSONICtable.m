% This script will upsample a SONIC table. E.g. upsample SONIC table with
% spline/makima interpolation, while linear nakeinterp1 interpolation is
% used at runtime of SONICrun.m
% Motivation: To improve accuracy there is the choice between increasing the effective
% table resolution (computationally intensive on the tabulation step) or to
% upsample with a higher order interpolation method

% In the current version, only one dimension will be upsampled

clear all; clc; 
tableVersion = '-v2';               % Note: spline, pchip and cubic need 4 samples in every direction: use singleton expansion (sx). (Makima only requires 2 samples)
modelName = 'RS';
NFS = 2;
UpsampleFactor = [4,2,5];           % Upsample factor of dim: [DeltaQm,phiQ,Qm]
UpsampleMethod = 'makima';
HomogeneousOutput = 1;              % To reduce memory requirements, HomogeneousOutput==1 gets rid of cells/structs/...
sredpfaf = 1;                       % If 1, pressure (p), frequency (f), bilayer radius (a) and coverage (f_BLS) 
                                    % are reduced to singleton dimensions 
                               
%% - 

loadStr = ['SONIC-' modelName '-QoscFourier' num2str(NFS) '-FourierIn' num2str(NFS) tableVersion '.mat'];
UpDimStr = '';
if any(UpsampleFactor~=1)
UpDimStr = [UpDimStr '-up'];
if UpsampleFactor(1) ~= 1, UpDimStr = [UpDimStr '_DeltaQm' num2str(UpsampleFactor(1))]; end
if UpsampleFactor(2) ~= 1, UpDimStr = [UpDimStr '_phiQ' num2str(UpsampleFactor(2))]; end
if UpsampleFactor(3) ~= 1, UpDimStr = [UpDimStr '_Qm' num2str(UpsampleFactor(3))]; end
UpDimStr = [UpDimStr '_' UpsampleMethod];
saveStr = [loadStr(1:end-4) UpDimStr];
if sredpfaf == 1
saveStr = [saveStr '_sredpfaf.mat'];
else
saveStr = [saveStr '.mat'];
end
else
if sredpfaf == 1
saveStr = [loadStr(1:end-4) '_sredpfaf.mat'];
else
error('sredpfaf = 0 and UpsampleFactor = 1. Data will not change');
end   
end

SONIC = load(loadStr); 
SONICtable = SONIC.SONICtable; clear SONIC;

QmRange = SONICtable.QmRange; 
USPaRange = SONICtable.USPaRange; USfreqRange = SONICtable.USfreqRange; 
aBLSRange = SONICtable.aBLSRange; fBLSRange = SONICtable.fBLSRange;
DeltaQmRange = SONICtable.DeltaQmRange; psiQRange = SONICtable.psiQRange;
if (sredpfaf)
USPaRange = USPaRange(1); USfreqRange = USfreqRange(1);
aBLSRange = aBLSRange(1); fBLSRange = fBLSRange(1);
SONICtable.USPaRange = USPaRange; SONICtable.USfreqRange = USfreqRange;
SONICtable.aBLSRange = aBLSRange; SONICtable.fBLSRange = fBLSRange;
end

if ~(sredpfaf)
resh = [numel(QmRange),numel(USPaRange),numel(USfreqRange),numel(aBLSRange),numel(fBLSRange),repmat(numel(DeltaQmRange),[1,NFS]),repmat(numel(psiQRange),[1,NFS])];
interpgrid = horzcat({QmRange},{USPaRange},{USfreqRange},{aBLSRange},{fBLSRange},repmat({DeltaQmRange},[1,NFS]),repmat({psiQRange},[1,NFS]));


reshUP = [UpsampleFactor(3)*numel(QmRange),numel(USPaRange),numel(USfreqRange),numel(aBLSRange),numel(fBLSRange),repmat(UpsampleFactor(1)*numel(DeltaQmRange),[1,NFS]),repmat(UpsampleFactor(2)*numel(psiQRange),[1,NFS])];
reshUPBack = [UpsampleFactor(3)*numel(QmRange),numel(USPaRange),numel(USfreqRange),numel(aBLSRange),numel(fBLSRange),(UpsampleFactor(1)*numel(DeltaQmRange))^NFS,(UpsampleFactor(2)*numel(psiQRange))^NFS];

DeltaQmRangeInterp = linspace(DeltaQmRange(1),DeltaQmRange(end),UpsampleFactor(1)*numel(DeltaQmRange));
psiQRangeInterp = linspace(psiQRange(1),psiQRange(end),UpsampleFactor(2)*numel(psiQRange));
QmRangeInterp = linspace(QmRange(1),QmRange(end),UpsampleFactor(3)*numel(QmRange));

querygrid = horzcat({QmRangeInterp},{USPaRange},{USfreqRange},{aBLSRange},{fBLSRange},repmat({DeltaQmRangeInterp},[1,NFS]),repmat({psiQRangeInterp},[1,NFS]));

else
resh = [numel(QmRange),repmat(numel(DeltaQmRange),[1,NFS]),repmat(numel(psiQRange),[1,NFS])];
interpgrid = horzcat({QmRange},repmat({DeltaQmRange},[1,NFS]),repmat({psiQRange},[1,NFS]));

reshUP = [UpsampleFactor(3)*numel(QmRange),repmat(UpsampleFactor(1)*numel(DeltaQmRange),[1,NFS]),repmat(UpsampleFactor(2)*numel(psiQRange),[1,NFS])];
reshUPBack = [UpsampleFactor(3)*numel(QmRange),(UpsampleFactor(1)*numel(DeltaQmRange))^NFS,(UpsampleFactor(2)*numel(psiQRange))^NFS];

DeltaQmRangeInterp = linspace(DeltaQmRange(1),DeltaQmRange(end),UpsampleFactor(1)*numel(DeltaQmRange));
psiQRangeInterp = linspace(psiQRange(1),psiQRange(end),UpsampleFactor(2)*numel(psiQRange));
QmRangeInterp = linspace(QmRange(1),QmRange(end),UpsampleFactor(3)*numel(QmRange));

querygrid = horzcat({QmRangeInterp},repmat({DeltaQmRangeInterp},[1,NFS]),repmat({psiQRangeInterp},[1,NFS]));
end
reshC = num2cell(resh); reshCones = cellfun(@(Z) ones(Z,1),reshC,'UniformOutput',0);
reshCUP = num2cell(reshUP); reshConesUP = cellfun(@(Z) ones(Z,1),reshCUP,'UniformOutput',0);

queryndgrid = cell(size(querygrid));
[queryndgrid{:}] = ndgrid(querygrid{:});

SONICfields = fieldnames(SONICtable);
SONICtables = SONICfields(~cellfun(@(X) contains(X,'Range'),SONICfields));

SONICtable.psiQRange = psiQRangeInterp;
SONICtable.DeltaQmRange = DeltaQmRangeInterp;
SONICtable.QmRange = QmRangeInterp;

for i = 1:numel(SONICtables)
if (sredpfaf)
SONICtable.(SONICtables{i}) = permute(SONICtable.(SONICtables{i})(:,1,1,1,1,:,:),[1 6 7 2 3 4 5]);
end
if ~strcmp(SONICtables{i},'cfit')
SONICtable.(SONICtables{i}) = reshape(interpn(interpgrid{:},reshape(SONICtable.(SONICtables{i}),resh),...
    queryndgrid{:},UpsampleMethod),reshUPBack);
else
    % Here, also fBLS (cfr. SONICrun_nanoMC_Qosc) 
    % As one line to circumvent large temp vars:
    
% X = permute(mat2cell(cell2mat(cellfun(@(X)permute(X',[2*NFS+6,(2:2*NFS+5),1]),reshape(SONICtable.(SONICtables{i}),resh),...
%     'UniformOutput',0)),reshC{:},ones(2*NFS+1,1)),[2*NFS+6,(1:2*NFS+5)]);           % Here, also fBLS (cfr. SONICrun_nanoMC_Qosc)
% XX = cellfun(@(Y) interpn(interpgrid{:},Y,queryndgrid{:},UpsampleMethod),X,'UniformOutput',0);
% 
% SONICtable.(SONICtables{i}) = reshape(cellfun(@(X) permute(X,[1, 2*NFS+6,(3:2*NFS+5),2]),...
%     mat2cell(cell2mat(permute(XX,[2*NFS+6,(2:2*NFS+5),1])),reshConesUP{:},2*NFS+1),'UniformOutput',0),reshUPBack);   
mindims = 0; if sredpfaf, mindims = 4; end
SONICtable.(SONICtables{i}) = reshape(cellfun(@(X) permute(X,[1, 2*NFS+6-mindims,(3:2*NFS+5-mindims),2]),...
    mat2cell(cell2mat(permute(cellfun(@(Y) interpn(interpgrid{:},Y,queryndgrid{:},UpsampleMethod),...
    permute(mat2cell(cell2mat(cellfun(@(X)permute(X',[2*NFS+6-mindims,(2:2*NFS+5-mindims),1]),...
    reshape(SONICtable.(SONICtables{i}),resh),'UniformOutput',0)),reshC{:},ones(2*NFS+1,1)),...
    [2*NFS+6-mindims,(1:2*NFS+5-mindims)]),'UniformOutput',0),[2*NFS+6-mindims,(2:2*NFS+5-mindims),1])),reshConesUP{:},2*NFS+1),...
    'UniformOutput',0),reshUPBack);
end
end

% Add the contraints a_X >= 0 and apb_X>=a_X (i.e. b_X >= 0)
SONICrates = sort(SONICfields(cellfun(@(X) contains(X,'a_')|contains(X,'apb_'),SONICfields)));
SONICgates = cellfun(@(X) X(3:end),SONICrates(cellfun(@(X) contains(X,'a_'),SONICrates)),'UniformOutput',0); 
for i = 1:numel(SONICgates)
SONICtable.(['a_' SONICgates{i}]) = max(SONICtable.(['a_' SONICgates{i}]),0);
SONICtable.(['apb_' SONICgates{i}]) = max(SONICtable.(['apb_' SONICgates{i}]),SONICtable.(['a_' SONICgates{i}])); 
end

if ~HomogeneousOutput, save(saveStr,'SONICtable','-v7.3'); end

if (HomogeneousOutput)
HO_SONIC = []; %#ok<UNRCH>
SONICtables = cat(1,SONICtables(~cellfun(@(X) contains(X,'cfit'),SONICtables)),{'cfit'});
for i = 1:numel(SONICtables)-1 
HO_SONIC = cat(8-mindims,HO_SONIC,SONICtable.(SONICtables{i})); 
end
for i = 1:2*NFS+1
HO_SONIC = cat(8-mindims,HO_SONIC,cellfun(@(X) X(i), SONICtable.(SONICtables{end})));
end
HO_SONIC = permute(HO_SONIC,[8-mindims,1:7-mindims]);

HO_SONICnames = SONICtables;
HO_SONICRange = struct;
for i = 1:numel(SONICfields)
if contains(SONICfields{i},'Range')
HO_SONICRange.(SONICfields{i}) = SONICtable.(SONICfields{i}); 
end
end
HO_SONICinfo = struct; HO_SONICinfo.SONICnames = HO_SONICnames; HO_SONICinfo.SONICRange = HO_SONICRange;

HO_saveStr = [saveStr(1:end-4) '_HO.mat'];
HO_saveStr_info = [saveStr(1:end-4) '_HO_info.mat'];

save(HO_saveStr,'HO_SONIC','-v7.3');
save(HO_saveStr_info,'HO_SONICinfo');
end
