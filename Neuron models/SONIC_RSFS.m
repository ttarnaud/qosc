function Out = SONIC_RSFS(ESi,USPaT,DISPLAY,tNICE,t,Q,h,m,n,p,Gna,Vna,Gk,Vk,Gm,Gl,Vl,f1Veff0,f1VeffPa,f1rt0,f1rtPa,SONICgates)
if DISPLAY == 1
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));  
end
rate = struct;
m = m*(m<=1&m>=0)+(m>1);
n =  n*(n<=1&n>=0)+(n>1);
p = p*(p<=1&p>=0)+(p>1);
h = h*(h<=1&h>=0)+(h>1);
rate.('m') =  m; rate.('n') = n; rate.('p') = p; rate.('h') = h;

if USPaT(t) == 0
Veff = f1Veff0(Q);
Out = [ESi(t)-10^(-3)*(Gna*m^3*h*(Veff-Vna)+Gk*n^4*(Veff-Vk)+Gm*p*(Veff-Vk)+Gl*(Veff-Vl));
cellfun(@(X) f1rt0.(['a_' X])(Q)-f1rt0.(['apb_' X])(Q)*rate.(X),SONICgates)];
else
Veff = f1VeffPa(Q);
Out = [ESi(t)-10^(-3)*(Gna*m^3*h*(Veff-Vna)+Gk*n^4*(Veff-Vk)+Gm*p*(Veff-Vk)+Gl*(Veff-Vl));
cellfun(@(X) f1rtPa.(['a_' X])(Q)-f1rtPa.(['apb_' X])(Q)*rate.(X),SONICgates)];    
end
end