function Out = SONIC_GPi(ESi,USPaT,DISPLAY,tNICE,t,Q,h,n,r,CA,Gna,Vna,Gk,Vk,Gl,Vl,...
    GT,VT,GCa,VCa,Gahp,minf,ainf,sinf,f1Veff0,f1VeffPa,f1rt0,f1rtPa,SONICgates)
if DISPLAY == 2
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));  
end
n = n*(n<=1)+(n>1);
h = h*(h<=1)+(h>1);
r = r*(r<=1)+(r>1);
rate.('n') = n; rate.('h') = h; rate.('r') = r;

if USPaT(t) == 0
Veff = f1Veff0(Q);
Out = [ESi(t)-10^(-3)*(Gl*(Veff-Vl)+Gna*minf(Veff).^3.*h.*(Veff-Vna)+...
   Gk*n.^4.*(Veff-Vk)+GT*ainf(Veff).^3.*r.*(Veff-VT)+...
   GCa.*sinf(Veff).^2.*(Veff-VCa)+Gahp*(Veff-Vk).*(CA./(CA+10)));
cellfun(@(X) f1rt0.(['a_' X])(Q)-f1rt0.(['apb_' X])(Q)*rate.(X),SONICgates)
-10^(-4)*(0.1*GT*ainf(Veff).^3.*r.*(Veff-VT)+0.1*GCa.*sinf(Veff).^2.*(Veff-VCa)+15*1000*CA)];
else
Veff = f1VeffPa(Q);
Out = [ESi(t)-10^(-3)*(Gl*(Veff-Vl)+Gna*minf(Veff).^3.*h.*(Veff-Vna)+...
   Gk*n.^4.*(Veff-Vk)+GT*ainf(Veff).^3.*r.*(Veff-VT)+...
   GCa.*sinf(Veff).^2.*(Veff-VCa)+Gahp*(Veff-Vk).*(CA./(CA+10)));
cellfun(@(X) f1rtPa.(['a_' X])(Q)-f1rtPa.(['apb_' X])(Q)*rate.(X),SONICgates)
-10^(-4)*(0.1*GT*ainf(Veff).^3.*r.*(Veff-VT)+0.1*GCa.*sinf(Veff).^2.*(Veff-VCa)+15*1000*CA)];
end
end