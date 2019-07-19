% This is the parameter 2 CSV-file for the PES-parallellized function 
% -> Intended for use on HPC 
% ------------------------------------------------------------------
% Input parameters are: 
% Tsim,MODE,USpstart,USpd,USfreq,USdc,USprf,USisppa,ESpstart,ESpd,...
%  ESdc,ESprf,ESisppa,PLOT,model,USibegin,USiend,SearchMode
clear all; %#ok<CLALL>

% Define parameter names for .csv file format
ParNames = {'Tsim','MODE','USpstart','USpd','USfreq','USdc','USprf','USisppa',...
    'ESpstart','ESpd','ESdc','ESprf','ESisppa','PLOT','model','USibegin',...
    'USiend','SearchMode'};
saveName = 'STN_EUS_4.csv';

% All combinations of the 18 input parameter arrays are generated.
% The Val-parameters indicate how the parameters are defined:
% Val = 1 -> each element in the array indicates a value of the parameter
% Val = 2 -> Par(2) linearly distributed elements between Par(1) and Par(3) 
% Val = 3 -> Par(3) logarithmicly distributed elements between Par(1) and Par(3)
TsimVal = 1;
MODEVal = 1;
USpstartVal = 1;
USpdVal = 1;
USfreqVal = 1;
USdcVal = 1;
USprfVal = 1;
USisppaVal = 2;
ESpstartVal = 1;
ESpdVal = 1;
ESdcVal = 1;
ESprfVal = 1;
ESisppaVal = 2;
PLOTVal = 1;
modelVal = 1;
USibeginVal = 1;
USiendVal = 1;
SearchModeVal = 1;
% ---COMBINE ALL VALUES IN PARVal-----------------------------------------
PARVal = [TsimVal,MODEVal,USpstartVal,USpdVal,USfreqVal,USdcVal,USprfVal,USisppaVal,...
    ESpstartVal,ESpdVal,ESdcVal,ESprfVal,ESisppaVal,PLOTVal,modelVal,USibeginVal,...
    USiendVal,SearchModeVal];
% ------------------------------------------------------------------------

Tsim = [1.5]; %#ok<*NBRAK>
MODE = [2];
USpstart = [500e-3];
USpd = [500e-3];
USfreq = [0.7e6];
USdc = [0.05];
USprf = [100];
USisppa = [0,11,500];
ESpstart = [500e-3];
ESpd = [500e-3];
ESdc = [1];
ESprf = [1];
ESisppa = [0,16,0.15];
PLOT = [0];
model = [9];
USibegin = [0];
USiend = [0];
SearchMode = [0];
% -------------COMBINE ALL PARAMETERS IN PAR cell--------------------------
PAR = {Tsim;MODE;USpstart;USpd;USfreq;USdc;USprf;USisppa;ESpstart;ESpd;ESdc;ESprf;...
    ESisppa;PLOT;model;USibegin;USiend;SearchMode};
% -------------------------------------------------------------------------

% CALCULATE EXPLICIT PARValues
PARValues = cell(length(PAR),1); NPARStep = zeros(length(PAR),1);
for i=1:length(PAR)
if PARVal(i) == 1
PARValues{i} = PAR{i};
else
if PAR{i}(1) == PAR{i}(3)
PARValues{i} = (PAR{i}(1));  
else
if PARVal(i) == 3
PARValues{i} =(PAR{i}(3)/PAR{i}(1)).^((0:1:(PAR{i}(2)-1))/(PAR{i}(2)-1))*PAR{i}(1);
elseif PARVal(i) == 2 
PARStep = (PAR{i}(3)-PAR{i}(1))/(PAR{i}(2)-1);
PARValues{i} = (PAR{i}(1):PARStep:PAR{i}(3));
end
end
end
NPARStep(i) = length(PARValues{i});
end

Config = zeros(prod(NPARStep),length(PAR));
digitVal = cumprod(NPARStep(end:-1:1));
digitVal = [digitVal(end-1:-1:1);1];
IndexList = cell(prod(NPARStep),1);
for i=0:prod(NPARStep)-1
    Digit=i; Index = ones(1,length(PAR));
for j=1:length(PAR)
while Digit>=digitVal(j)
Index(j) = Index(j)+1;
Digit = Digit-digitVal(j);
end
end
IndexList{i+1} = Index;
end

for i=1:size(Config,1)
Config(i,:) = cellfun(@(x,y) x(y),PARValues',num2cell(IndexList{i}));
end
Config = vertcat(ParNames,num2cell(Config)); 
xlswrite(saveName,Config);