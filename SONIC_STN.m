function Out = SONIC_STN(ESi,USPaT,DISPLAY,tNICE,t,Q,a,b,c,d1,h,m,n,p,q,r,d2,cCai,...
    Gna,Vna,Gk,Vk,Gl,Vl,GT,fVCa,GCa,GA,GL,Far,tauCa,f1Veff0,f1VeffPa,f1rt0,f1rtPa,SONICgates)
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
rate.('p') = p*(p<=1)+(p>1);
rate.('h') = h*(h<=1)+(h>1);
rate.('q') = q*(q<=1)+(q>1);
rate.('r') = r*(r<=1)+(r>1);
rate.('a') = a*(a<=1)+(a>1);
rate.('b') = b*(b<=1)+(b>1);
rate.('c') = c*(c<=1)+(c>1);
rate.('d1') = d1*(d1<=1)+(d1>1);
rate.('d2') = d2*(d2<=1)+(d2>1);

kcCai = 10^3*cCai;          % A trick to improve the Jacobian (and speed-up the solver)

if USPaT(t) == 0
Veff = f1Veff0(Q);
Out = [ESi(t)-10^(-3)*(Gl*(Veff-Vl)+Gna*m.^3.*h.*(Veff-Vna)+Gk*n.^4.*(Veff-Vk)+...
   GT*p.^2.*q.*(Veff-fVCa(cCai))+GCa.*r.^2.*(Veff-Vk)+GA*a.^2.*b.*(Veff-Vk)+...
   +GL*c.^2.*d1.*d2.*(Veff-fVCa(cCai)));
cellfun(@(X) f1rt0.(['a_' X])(Q)-f1rt0.(['apb_' X])(Q)*rate.(X),SONICgates);
(rinf(kcCai)-r)/taur(kcCai)
(d2inf(kcCai)-d2)/taud2(kcCai);
-(0.001*(GT*p.^2.*q.*(V-fVCa(cCai))+...
GL*c.^2.*d1.*d2.*(V-fVCa(cCai))))/(2*Far*(10236*10^(-9)))-cCai/(tauCa)]; % Gain factor as in Kamaruvelu et al. (2016)
else
Veff = f1VeffPa(Q);
Out = [ESi(t)-10^(-3)*(Gl*(Veff-Vl)+Gna*m.^3.*h.*(Veff-Vna)+Gk*n.^4.*(Veff-Vk)+...
   GT*p.^2.*q.*(Veff-fVCa(cCai))+GCa.*r.^2.*(Veff-Vk)+GA*a.^2.*b.*(Veff-Vk)+...
   +GL*c.^2.*d1.*d2.*(Veff-fVCa(cCai)));
cellfun(@(X) f1rtPa.(['a_' X])(Q)-f1rtPa.(['apb_' X])(Q)*rate.(X),SONICgates);
(rinf(kcCai)-r)/taur(kcCai)
(d2inf(kcCai)-d2)/taud2(kcCai);
-(0.001*(GT*p.^2.*q.*(V-fVCa(cCai))+...
GL*c.^2.*d1.*d2.*(Veff-fVCa(cCai))))/(2*Far*(10236*10^(-9)))-cCai/(tauCa)]; % Gain factor as in Kamaruvelu et al. (2016)
end
end