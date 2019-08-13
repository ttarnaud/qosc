function Out = SONIC_TC(ESi,USPaT,DISPLAY,tNICE,t,Q,m,n,h,s,u,w,wLock,hProtein,cCai,...
    Gna,Vna,Gk,Vk,GT,fVCa,Gl,Vl,GKL,Gh,ginc,Vh,k1,k2,k3,k4,Far,deffCa,tauCa,f1Veff0,f1VeffPa,f1rt0,f1rtPa)
if DISPLAY == 2
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg)); 
end
rate = struct;
rate.('m') = m*(m<=1)+(m>1); 
rate.('n') = n*(n<=1)+(n>1);
rate.('h') = h*(h<=1)+(h>1);
rate.('s') = s*(s<=1)+(s>1);
rate.('u') = u*(u<=1)+(u>1);
rate.('w') = w*(w<=1)+(w>1);

if USPaT(t) == 0
Veff = f1Veff0(Q);
Out = [ESi(t)-10^(-3)*(Gna*m^3*h*(Veff-Vna)+Gk*n^4*(Veff-Vk)+GT*s^2*u*(Veff-fVCa(cCai))+Gl*(Veff-Vl)+...
    GKL*(Veff-Vk)+Gh*(w+ginc*wLock)*(Veff-Vh));
am(V)-(ampbm(V))*m;
an(V)-(anpbn(V))*n;
ah(V)-(ahpbh(V))*h;
(sinf(V)-s)/taus(V);
(uinf(V)-u)/tauu(V);
(winf(V)*(1-wLock)-w)/tauw(V);
k3*(w*hProtein)-k4*wLock;
k1*((1-hProtein)*cCai^2)-k2*hProtein;
-(10^(-3)*GT*s^2*u*(V-fVCa(cCai)))/(2*Far*deffCa)-cCai/tauCa];
end