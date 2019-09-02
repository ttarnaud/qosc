function Out = SONIC_TC(ESi,USPaT,DISPLAY,tNICE,t,Q,h,m,n,s,u,w,wLock,hProtein,cCai,...
    Gna,Vna,Gk,Vk,GT,fVCa,Gl,Vl,GKL,Gh,ginc,Vh,k1,k2,k3,k4,Far,deffCa,tauCa,f1Veff0,f1VeffPa,f1rt0,f1rtPa,SONICgates)
if DISPLAY == 2
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg)); 
end
rate = struct;
m = m*(m<=1&m>=0)+(m>1);  
n = n*(n<=1&n>=0)+(n>1);
h =  h*(h<=1&h>=0)+(h>1);
s = s*(s<=1&s>=0)+(s>1);
u = u*(u<=1&u>=0)+(u>1);
w = w*(w<=1&w>=0)+(w>1);
rate.('m') = m; rate.('n') = n; rate.('h') = h;
rate.('s') = s; rate.('u') = u; rate.('w') = w;

if USPaT(t) == 0
Veff = f1Veff0(Q);
Out = [ESi(t)-10^(-3)*(Gna*m^3*h*(Veff-Vna)+Gk*n^4*(Veff-Vk)+GT*s^2*u*(Veff-fVCa(cCai))+Gl*(Veff-Vl)+...
    GKL*(Veff-Vk)+Gh*(w+ginc*wLock)*(Veff-Vh));
cellfun(@(X) f1rt0.(['a_' X])(Q)-f1rt0.(['apb_' X])(Q)*rate.(X),SONICgates);
(winf(Veff)*(1-wLock)-w)/tauw(Veff);
k3*(w*hProtein)-k4*wLock;
k1*((1-hProtein)*cCai^2)-k2*hProtein;
-(10^(-3)*GT*s^2*u*(Veff-fVCa(cCai)))/(2*Far*deffCa)-cCai/tauCa];
else
Veff = f1VeffPa(Q);
Out = [ESi(t)-10^(-3)*(Gna*m^3*h*(Veff-Vna)+Gk*n^4*(Veff-Vk)+GT*s^2*u*(Veff-fVCa(cCai))+Gl*(Veff-Vl)+...
    GKL*(Veff-Vk)+Gh*(w+ginc*wLock)*(Veff-Vh));
cellfun(@(X) f1rtPa.(['a_' X])(Q)-f1rtPa.(['apb_' X])(Q)*rate.(X),SONICgates);
(winf(Veff)*(1-wLock)-w)/tauw(Veff);
k3*(w*hProtein)-k4*wLock;
k1*((1-hProtein)*cCai^2)-k2*hProtein;
-(10^(-3)*GT*s^2*u*(Veff-fVCa(cCai)))/(2*Far*deffCa)-cCai/tauCa];
end
end