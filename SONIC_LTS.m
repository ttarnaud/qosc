function Out = SONIC_LTS(ESi,USPaT,DISPLAY,tNICE,t,Q,h,m,n,p,s,u,Gna,Vna,Gk,Vk,Gm,GT,VCa,Gl,Vl,f1Veff0,f1VeffPa,f1rt0,f1rtPa,SONICgates)
if DISPLAY == 1
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg)); 
end
rate = struct;
m = m*(m<=1)+(m>1);
n = n*(n<=1)+(n>1);
p = p*(p<=1)+(p>1);
h = h*(h<=1)+(h>1);
s = s*(s<=1)+(s>1);
u = u*(u<=1)+(u>1);
rate.('m') = m; rate.('n') = n; rate.('p') = p;
rate.('h') = h; rate.('s') = s; rate.('u') = u;

if USPaT(t) == 0
Veff = f1Veff0(Q);
Out = [ESi(t)-10^(-3)*(Gna*m^3*h*(Veff-Vna)+Gk*n^4*(Veff-Vk)+Gm*p*(Veff-Vk)+GT*s^2*u*(Veff-VCa)+Gl*(Veff-Vl));
cellfun(@(X) f1rt0.(['a_' X])(Q)-f1rt0.(['apb_' X])(Q)*rate.(X),SONICgates)];
else
Veff = f1VeffPa(Q);
Out = [ESi(t)-10^(-3)*(Gna*m^3*h*(Veff-Vna)+Gk*n^4*(Veff-Vk)+Gm*p*(Veff-Vk)+GT*s^2*u*(Veff-VCa)+Gl*(Veff-Vl));
cellfun(@(X) f1rtPa.(['a_' X])(Q)-f1rtPa.(['apb_' X])(Q)*rate.(X),SONICgates)];
end
end