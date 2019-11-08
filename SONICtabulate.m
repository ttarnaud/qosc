function SONICtabulate
ModelN = 1;
aBLSRange = [16,32,64]*10^(-9); 		% Sonophore radius (m)
USfreqRange = [0.02, 0.1, 0.5, 1, 2, 3, 4]*10^(6); % Ultrasonic frequency (Hz)
USPaRange = [0, logspace(log10(100),log10(600000),50)];  % Ultrasonic pressure (Pa)
fQmRange = @(Qm0) (Qm0-25*10^(-5):1*10^(-5):50*10^(-5)); % Membrane charge (C/m^2)
fBLSRange = 1;                                           %#ok<*NASGU> % Membrane coverage
DeltaQmRange = 0;
saveStrAdd = '';

% DeltaQmRange tables (outcomment for normal tables)
DeltaQmRange = (0*10^(-5):25*10^(-5):500*10^(-5));  
USfreqRange = 0.5*10^(6);            % Ultrasonic frequency [Hz]
aBLSRange = 32*10^(-9);              % Sonophore radius (m)  
saveStrAdd = [saveStrAdd,'-Qosc'];

% % For Fig. 10 from Lemaire et al. we need tables with multiple membrane
% % coverages -- Out-comment for normal tables.
% fBLSRange = (0.05:0.05:1);
% USfreqRange = 0.5*10^(6);            % Ultrasonic frequency [Hz]
% aBLSRange = 32*10^(-9);              % Sonophore radius (m)   
% saveStrAdd = [saveStrAdd '-xfs'];

% Correct for charge redistribution if fBLS < 1 in the electrostatic force (corrPec) 
corrPec = 0;                        % Bool: if 1, the table-name will contain "-corrPec"
% --------------------------------------------------------------------------------
ME = ''; futIND = 0;
try p = gcp;
catch ME
end
if isempty(ME)
ParallelON = 1; 
disp('Parallel computing toolbox available - running in gcp');
q = p.FevalQueue;
elseif strcmp(ME.identifier,'parallel:cluster:LicenseUnavailable')
ParallelON = 0;
disp('Parallel computing toolbox not available - running sequentially');
end
for iModel = 1:length(ModelN)
MODEL = ModelN(iModel);
Cm0 = 0.01; 				% Rest capacitance (F/m^2);   
switch MODEL
   case 1, MODELstr='RS'; Vm0 = -71.9;  % mV
   case 2, MODELstr='FS'; Vm0 = -71.4;
   case 3, MODELstr='LTS'; Vm0 = -54; 
   case 4, MODELstr='TC'; Vm0 = -63.4;          
   case 5, MODELstr='RE'; Vm0 = -89.5;
   case 6, MODELstr='RS_FVX'; Vm0 = -70.4;
   case 7, MODELstr='FS_FVX'; Vm0 = -70; 
   case 8, MODELstr='LTS_CAX'; Vm0 = -84.6; 
   case 9, MODELstr='STN'; Vm0=-58;       
   case 10, MODELstr='Th-RT'; Vm0=-65; 
   case 11, MODELstr='GPi'; Vm0 = -65;      
   case 12, MODELstr='GPe'; Vm0 = -65;
   case 13, MODELstr='Str-MSN'; Vm0=-87;
   case 14, MODELstr='HH'; Vm0 = -70;
end
QmRange = fQmRange(10^(-3)*Vm0*Cm0);
SONICInitM = zeros(numel(QmRange),numel(USPaRange),numel(USfreqRange),numel(aBLSRange),numel(fBLSRange),numel(DeltaQmRange));
SONICtable = struct('Zeff',SONICInitM,'Veff',SONICInitM,'ngend',SONICInitM,'Cmeff',SONICInitM,'ampV',SONICInitM,'DeltaPhi',SONICInitM);
ProgressIncr1 = 0; ProgressIncr2 = 0; reverseStr1 = ''; reverseStr2 = '';
if ParallelON == 1
for ifBLS = (1:length(fBLSRange))
fBLS = fBLSRange(ifBLS);
    for iaBLS = (1:length(aBLSRange))
    aBLS = aBLSRange(iaBLS);
        for iUSfreq = (1:length(USfreqRange))
        USfreq = USfreqRange(iUSfreq);
            for iUSPa = (1:length(USPaRange))
            USPa = USPaRange(iUSPa);			
                for iQm = (1:length(QmRange)) 
                Qm = QmRange(iQm);
                    for iDeltaQm = (1:length(DeltaQmRange))
                    DeltaQm = DeltaQmRange(iDeltaQm);
                    
                    ProgressIncr1 = ProgressIncr1+1;
                    Progress = 100*ProgressIncr1/numel(SONICInitM);

                    msg = sprintf('Submitting futures to gcp. Progress: %3.5f', Progress);
                    fprintf([reverseStr1, msg]);
                    reverseStr1 = repmat(sprintf('\b'), 1, length(msg)); 

                    FutureResults(iQm,iUSPa,iUSfreq,iaBLS,ifBLS,iDeltaQm) = parfeval(p,@SONICcalc,8,MODEL,Qm,USPa,USfreq,aBLS,fBLS,corrPec,DeltaQm); %#ok<AGROW> 
                    end
                end
            end
        end
    end
end
end
fprintf('\n');
for ifBLS = (1:length(fBLSRange))
    for iaBLS = (1:length(aBLSRange))
        for iUSfreq = (1:length(USfreqRange))
            for iUSPa = (1:length(USPaRange))	
                for iQm = (1:length(QmRange))
                    for iDeltaQm = (1:length(DeltaQmRange))
                    ProgressIncr2 = ProgressIncr2+1;
                    Progress = 100*ProgressIncr2/numel(SONICInitM);

                    msg = sprintf('Simulating and fetching results. Progress: %3.5f', Progress); 
                    fprintf([reverseStr2, msg]);
                    reverseStr2 = repmat(sprintf('\b'), 1, length(msg)); 

                    if ParallelON == 1
                    TIMEOUT = 60*60*2; % (s)
                    [futIND,Zeff,Veff,a_i_eff,apb_i_eff,ngend,Cmeff,ampV,DeltaPhi] = fetchNext(FutureResults,TIMEOUT);
                    if isempty(futIND)
                       fprintf('\n');
                       disp('TIME OUT ERROR. Result could not be fetched in time');
                       disp('Running futures');
                       disp(q.RunningFutures);
                       disp('Input arguments');
                       for i = 1:length(q.RunningFutures)
                       InputArgs(i,:) = q.RunningFutures(i).InputArguments; %#ok<AGROW>
                       end
                       disp(InputArgs);
                       error('---------------------------');               
                    end
                    else
                    fBLS = fBLSRange(ifBLS);
                    aBLS = aBLSRange(iaBLS);
                    USfreq = USfreqRange(iUSfreq);
                    USPa = USPaRange(iUSPa);                
                    Qm = QmRange(iQm);
                    DeltaQm = DeltaQmRange(iDeltaQm);
                    futIND = futIND+1;
                    [Zeff,Veff,a_i_eff,apb_i_eff,ngend,Cmeff,ampV,DeltaPhi] = SONICcalc(MODEL,Qm,USPa,USfreq,aBLS,fBLS,corrPec,DeltaQm);
                    end
                    SONICtable.Zeff(ind2sub(size(SONICInitM),futIND-(iModel-1)*numel(SONICInitM))) = Zeff;
                    SONICtable.Veff(ind2sub(size(SONICInitM),futIND-(iModel-1)*numel(SONICInitM))) = Veff;
                    SONICtable.ngend(ind2sub(size(SONICInitM),futIND-(iModel-1)*numel(SONICInitM))) = ngend;
                    SONICtable.Cmeff(ind2sub(size(SONICInitM),futIND-(iModel-1)*numel(SONICInitM))) = Cmeff;  
                    SONICtable.ampV(ind2sub(size(SONICInitM),futIND-(iModel-1)*numel(SONICInitM))) = ampV;  
                    SONICtable.DeltaPhi(ind2sub(size(SONICInitM),futIND-(iModel-1)*numel(SONICInitM))) = DeltaPhi;  

                    alphaFNs = fieldnames(a_i_eff); 
                    alphapbetaFNs = fieldnames(apb_i_eff);

                        for i = 1:numel(alphaFNs)
                        if ~isfield(SONICtable,['a_' alphaFNs{i}]), SONICtable.(['a_' alphaFNs{i}]) = SONICInitM; end 
                        SONICtable.(['a_' alphaFNs{i}])(ind2sub(size(SONICInitM),futIND-(iModel-1)*numel(SONICInitM))) = a_i_eff.(alphaFNs{i});
                        end
                        for i = 1:numel(alphapbetaFNs)
                        if ~isfield(SONICtable,['apb_' alphapbetaFNs{i}]), SONICtable.(['apb_' alphapbetaFNs{i}]) = SONICInitM; end 
                        SONICtable.(['apb_' alphapbetaFNs{i}])(ind2sub(size(SONICInitM),futIND-(iModel-1)*numel(SONICInitM))) = apb_i_eff.(alphapbetaFNs{i});
                        end		
                    end
                end
            end
        end
    end
end

fn = fieldnames(SONICtable);  
% If a dimension of the 5D table is a singleton, interpn will fail, because
% it needs at least 2 points in each direction. 
% To streamline calculations in SONICrun, we solve this by modifying the
% tables by replicating the matrix in this dimension.
if length(QmRange) == 1
QmRange = horzcat(QmRange,QmRange(1)+eps(QmRange(1))); %#ok<*AGROW>
for i = 1:length(fn)
SONICtable.(fn{i}) = repmat(SONICtable.(fn{i}),[2 1 1 1 1 1]);
end
end
if length(USPaRange) == 1
USPaRange = horzcat(USPaRange,USPaRange(1)+eps(USPaRange(1)));
for i = 1:length(fn)
SONICtable.(fn{i}) = repmat(SONICtable.(fn{i}),[1 2 1 1 1 1]);
end
end
if length(USfreqRange) == 1
USfreqRange = horzcat(USfreqRange,USfreqRange(1)+eps(USfreqRange(1)));
for i = 1:length(fn)
SONICtable.(fn{i}) = repmat(SONICtable.(fn{i}),[1 1 2 1 1 1]);
end
end
if length(aBLSRange) == 1
aBLSRange = horzcat(aBLSRange,aBLSRange(1)+eps(aBLSRange(1)));
for i = 1:length(fn)
SONICtable.(fn{i}) = repmat(SONICtable.(fn{i}),[1 1 1 2 1 1]);
end
end
if length(fBLSRange) == 1
fBLSRange = horzcat(fBLSRange,fBLSRange(1)+eps(fBLSRange(1)));
for i = 1:length(fn)
SONICtable.(fn{i}) = repmat(SONICtable.(fn{i}),[1 1 1 1 2 1]);
end
end
if length(DeltaQmRange) == 1
DeltaQmRange = horzcat(DeltaQmRange,DeltaQmRange(1)+eps(DeltaQmRange(1)));
for i = 1:length(fn)
SONICtable.(fn{i}) = repmat(SONICtable.(fn{i}),[1 1 1 1 1 2]);
end
end
    
SONICtable.QmRange = QmRange;
SONICtable.aBLSRange = aBLSRange;
SONICtable.USfreqRange = USfreqRange; 
SONICtable.USPaRange = USPaRange; 
SONICtable.fBLSRange = fBLSRange;
SONICtable.DeltaQmRange = DeltaQmRange;

if (corrPec)
  corrPecStr = '-corrPec';
else 
  corrPecStr = ''; %#ok<*UNRCH>
end
saveStr = ['SONIC-' MODELstr saveStrAdd corrPecStr '.mat'];
save(saveStr,'SONICtable');

end
end
