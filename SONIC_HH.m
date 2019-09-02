function Out = SONIC_HH(ESi,USPaT,DISPLAY,tNICE,t,Q,h,m,n,Gna,Vna,Gk,Vk,Gl,Vl,f1Veff0,f1VeffPa,f1rt0,f1rtPa,SONICgates)
if DISPLAY == 2
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));  
end
rate = struct;
m = m*(m<=1)+(m>1);
n =  n*(n<=1)+(n>1);
h = h*(h<=1)+(h>1);
rate.('m') = m; rate.('n') = n; rate.('h') = h;

if USPaT(t) == 0
Veff = f1Veff0(Q);
Out = [ESi(t)-10^(-3)*(Gl*(Veff-Vl)+Gna*m.^3.*h.*(Veff-Vna)+...
    Gk*n.^4.*(Veff-Vk));
    cellfun(@(X) f1rt0.(['a_' X])(Q)-f1rt0.(['apb_' X])(Q)*rate.(X),SONICgates)];
else
Veff = f1VeffPa(Q);   
Out = [ESi(t)-10^(-3)*(Gl*(Veff-Vl)+Gna*m.^3.*h.*(Veff-Vna)+...
    Gk*n.^4.*(Veff-Vk));
    cellfun(@(X) f1rtPa.(['a_' X])(Q)-f1rtPa.(['apb_' X])(Q)*rate.(X),SONICgates)];
end
end