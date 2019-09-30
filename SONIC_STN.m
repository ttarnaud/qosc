function Out = SONIC_STN(ESi,USPaT,DISPLAY,tNICE,t,Q,a,b,c,d1,h,m,n,p,q,r,d2,cCai,...
    Gna,Vna,Gk,Vk,Gl,Vl,GT,fVCa,GCa,GA,GL,Far,tauCa,f1Veff0,f1VeffPa,f1rt0,f1rtPa,rinf,d2inf,taur,taud2,SONICgates)
if DISPLAY == 1
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));  
end
rate = struct;
m = m*(m<=1&m>=0)+(m>1);
n = n*(n<=1&n>=0)+(n>1);
p = p*(p<=1&p>=0)+(p>1);
h = h*(h<=1&h>=0)+(h>1);
q = q*(q<=1&q>=0)+(q>1);
r = r*(r<=1&r>=0)+(r>1);
a =  a*(a<=1&a>=0)+(a>1);
b = b*(b<=1&b>=0)+(b>1);
c =  c*(c<=1&c>=0)+(c>1);
d1 = d1*(d1<=1&d1>=0)+(d1>1);
d2 = d2*(d2<=1&d2>=0)+(d2>1);
rate.('m') = m; rate.('n') = n; rate.('p') = p; rate.('h') = h;
rate.('q') = q; rate.('r') = r; rate.('a') = a; rate.('b') = b;
rate.('c') = c; rate.('d1') = d1; rate.('d2') = d2;

kcCai = 10^3*cCai;          % A trick to improve the Jacobian (and speed-up the solver)

if USPaT(t) == 0
Veff = f1Veff0(Q);
Out = [ESi(t)-10^(-3)*(Gl*(Veff-Vl)+Gna*m.^3.*h.*(Veff-Vna)+Gk*n.^4.*(Veff-Vk)+...
   GT*p.^2.*q.*(Veff-fVCa(cCai))+GCa.*r.^2.*(Veff-Vk)+GA*a.^2.*b.*(Veff-Vk)+...
   +GL*c.^2.*d1.*d2.*(Veff-fVCa(cCai)));
cellfun(@(X) f1rt0.(['a_' X])(Q)-f1rt0.(['apb_' X])(Q)*rate.(X),SONICgates);
(rinf(kcCai)-r)/taur(kcCai)
(d2inf(kcCai)-d2)/taud2(kcCai);
-(0.001*(GT*p.^2.*q.*(Veff-fVCa(cCai))+...
GL*c.^2.*d1.*d2.*(Veff-fVCa(cCai))))/(2*Far*(10236*10^(-9)))-cCai/(tauCa)]; % Gain factor as in Kamaruvelu et al. (2016)
else
Veff = f1VeffPa(Q);
Out = [ESi(t)-10^(-3)*(Gl*(Veff-Vl)+Gna*m.^3.*h.*(Veff-Vna)+Gk*n.^4.*(Veff-Vk)+...
   GT*p.^2.*q.*(Veff-fVCa(cCai))+GCa.*r.^2.*(Veff-Vk)+GA*a.^2.*b.*(Veff-Vk)+...
   +GL*c.^2.*d1.*d2.*(Veff-fVCa(cCai)));
cellfun(@(X) f1rtPa.(['a_' X])(Q)-f1rtPa.(['apb_' X])(Q)*rate.(X),SONICgates);
(rinf(kcCai)-r)/taur(kcCai)
(d2inf(kcCai)-d2)/taud2(kcCai);
-(0.001*(GT*p.^2.*q.*(Veff-fVCa(cCai))+...
GL*c.^2.*d1.*d2.*(Veff-fVCa(cCai))))/(2*Far*(10236*10^(-9)))-cCai/(tauCa)]; % Gain factor as in Kamaruvelu et al. (2016)
end
end