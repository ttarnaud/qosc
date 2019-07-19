function Out = SimplNICETC(ESi,DISPLAY,tNICE,t,Q,m,n,h,s,u,w,wLock,hProtein,cCai,...
CmR,Gna,Vna,Gk,Vk,GT,fVCa,Gl,Vl,GKL,Gh,ginc,Vh,am,bm,an,bn,ah,...
    bh,sinf,taus,uinf,tauu,winf,tauw,k1,k2,k3,k4,Far,deffCa,tauCa)
if DISPLAY == 2
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg)); 
end
m = m*(m<=1)+(m>1); 
n = n*(n<=1)+(n>1);
h = h*(h<=1)+(h>1);
s = s*(s<=1)+(s>1);
u = u*(u<=1)+(u>1);
w = w*(w<=1)+(w>1);

V = 10^(3)*Q/CmR(t);
Out = [ESi(t)-10^(-3)*(Gna*m^3*h*(V-Vna)+Gk*n^4*(V-Vk)+GT*s^2*u*(V-fVCa(cCai))+Gl*(V-Vl)+...
    GKL*(V-Vk)+Gh*(w+ginc*wLock)*(V-Vh));
am(V)-(am(V)+bm(V))*m;
an(V)-(an(V)+bn(V))*n;
ah(V)-(ah(V)+bh(V))*h;
(sinf(V)-s)/taus(V);
(uinf(V)-u)/tauu(V);
(winf(V)*(1-wLock)-w)/tauw(V);
k3*(w*hProtein)-k4*wLock;
k1*((1-hProtein)*cCai^2)-k2*hProtein;
-(10^(-3)*GT*s^2*u*(V-fVCa(cCai)))/(2*Far*deffCa)-cCai/tauCa];
end