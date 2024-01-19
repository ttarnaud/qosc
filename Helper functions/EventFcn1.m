function [value,isterminal,direction] = EventFcn1(t,Z,dZ,na,freq,dt,CORRTres,PerCheckMax) %#ok<INUSL>
isterminal = 1; direction = 1; value = -1;
global temp; global tempi; global maxlags;
if ~tempi||(t>=(temp{tempi,1}+dt))    % Save time by immediately sampling with dt (however less accurate)
    tempi = tempi+1;
    temp(tempi,:) = {t, Z, na};
end
if (t-temp{1,1})>=PerCheckMax/freq 
    tempMat = cell2mat(temp);
    % 2 periods of data -> check for periodicity (if periodic -> value=1)
firstIn = find((t-tempMat(:,1))>=(PerCheckMax/freq),1,'last');
tempor = tempMat(firstIn:end,:);     % Extract last two periods for analysis
TimeLine = tempor(:,1)-tempor(1,1); % Shift zero in timeline to first value of (Z,na)
% Resample the Z and na values, to allow calculation of autocorrelation
% Resampling is necessary, because the VSVO-ode113 solver has variable
% stepsizes, rendering the autocorrelation meaningless
% + the VSVO solver has the habit to regress locally and slightly in time,
% this behaviour vanishes after resampling
SampleT = (0:dt:PerCheckMax/freq)'; 
Ztemp = interp1(TimeLine,tempor(:,2),SampleT,'spline','extrap');
natemp = interp1(TimeLine,tempor(:,3),SampleT,'spline','extrap');
% Calculate the unbiased normalized autocorrelation
[acorrZ,~] = xcorr(Ztemp-mean(Ztemp),'unbiased');
[acorrna,lags] = xcorr(natemp-mean(natemp),'unbiased');
lags = lags*dt;
acorrZ = acorrZ/acorrZ(ceil(length(acorrZ)/2));  % Normalize autocorrelation to 1
acorrna = acorrna/acorrna(ceil(length(acorrna)/2));
lags2 = lags(lags>=(1/freq-dt)&lags<=(PerCheckMax/(2*freq)+dt));   % +- dt for the case PerCheckMax == 2
acorrZ2 = acorrZ(lags>=(1/freq-dt)&lags<=(PerCheckMax/(2*freq)+dt));
acorrna2 = acorrna(lags>=(1/freq-dt)&lags<=(PerCheckMax/(2*freq)+dt));

[maxacorrZ,Imax] = max(acorrZ2);
maxacorrna = acorrna2(Imax); maxlags = lags2(Imax);

if maxacorrZ>=CORRTres&&maxacorrna>=CORRTres
value = 1;
end
end
end