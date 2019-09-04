function Out = SONIC_ThRT(ESi,USPaT,DISPLAY,tNICE,t,Q,h,r,Gna,Vna,Gk,Vk,Gl,Vl,...
    GT,VT,minf,pinf,f1Veff0,f1VeffPa,f1rt0,f1rtPa,SONICgates)
if DISPLAY == 1
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));  
end
rate = struct;
h = h*(h<=1&h>=0)+(h>1);
r = r*(r<=1&r>=0)+(r>1);
rate.('h') = h; rate.('r') = r;

if USPaT(t) == 0
Veff = f1Veff0(Q);
Out = [ESi(t)-10^(-3)*(Gl*(Veff-Vl)+Gna*minf(Veff).^3.*h.*(Veff-Vna)+...
    Gk*(0.75*(1-h)).^4.*(Veff-Vk)+GT*pinf(Veff).^2.*r.*(Veff-VT));
cellfun(@(X) f1rt0.(['a_' X])(Q)-f1rt0.(['apb_' X])(Q)*rate.(X),SONICgates)];
else
Veff = f1VeffPa(Q);
Out = [ESi(t)-10^(-3)*(Gl*(Veff-Vl)+Gna*minf(Veff).^3.*h.*(Veff-Vna)+...
    Gk*(0.75*(1-h)).^4.*(Veff-Vk)+GT*pinf(Veff).^2.*r.*(Veff-VT));
cellfun(@(X) f1rtPa.(['a_' X])(Q)-f1rtPa.(['apb_' X])(Q)*rate.(X),SONICgates)];
end
end