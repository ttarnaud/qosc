function Out = SONIC_STN_nanoMC(ESi,USPaT,DISPLAY,tNICE,t,Q1,a1,b1,c1,d11,h1,m1,n1,p1,q1,r1,d21,cCai1,Q2,a2,b2,c2,d12,h2,m2,n2,p2,q2,r2,d22,cCai2,...
    Gna,Vna,Gk,Vk,Gl,Vl,GT,fVCa,GCa,GA,GL,Far,tauCa,f1Veff0,f1VeffPa,f1rt0,f1rtPa,rinf,d2inf,taur,taud2,f1rtV,SONICgates,Cm0,aBLS,fBLS,RSI)
if DISPLAY == 1
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>)
if DISPLAY == 1
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));  
end
rate1 = struct; rate2 = struct;
m1 = m1*(m1<=1&m1>=0)+(m1>1); m2 = m2*(m2<=1&m2>=0)+(m2>1);
n1 = n1*(n1<=1&n1>=0)+(n1>1); n2 = n2*(n2<=1&n2>=0)+(n2>1);
p1 = p1*(p1<=1&p1>=0)+(p1>1); p2 = p2*(p2<=1&p2>=0)+(p2>1);
h1 = h1*(h1<=1&h1>=0)+(h1>1); h2 = h2*(h2<=1&h2>=0)+(h2>1);
q1 = q1*(q1<=1&q1>=0)+(q1>1); q2 = q2*(q2<=2&q2>=0)+(q2>1);
r1 = r1*(r1<=1&r1>=0)+(r1>1); r2 = r2*(r2<=1&r2>=0)+(r2>1);
a1 =  a1*(a1<=1&a1>=0)+(a1>1); a2 =  a2*(a2<=1&a2>=0)+(a2>1);
b1 = b1*(b1<=1&b1>=0)+(b1>1); b2 = b2*(b2<=1&b2>=0)+(b2>1);
c1 =  c1*(c1<=1&c1>=0)+(c1>1); c2 =  c2*(c2<=1&c2>=0)+(c2>1);
d11 = d11*(d11<=1&d11>=0)+(d11>1); d12 = d12*(d12<=1&d12>=0)+(d12>1);
d21 = d21*(d21<=1&d21>=0)+(d21>1); d22 = d22*(d22<=1&d22>=0)+(d22>1);
rate1.('m') = m1; rate1.('n') = n1; rate1.('p') = p1; rate1.('h') = h1;
rate1.('q') = q1; rate1.('r') = r1; rate1.('a') = a1; rate1.('b') = b1;
rate1.('c') = c1; rate1.('d1') = d11; rate1.('d2') = d21;
rate2.('m') = m2; rate2.('n') = n2; rate2.('p') = p2; rate2.('h') = h2;
rate2.('q') = q2; rate2.('r') = r2; rate2.('a') = a2; rate2.('b') = b2;
rate2.('c') = c2; rate2.('d1') = d12; rate2.('d2') = d22;

kcCai1 = 10^3*cCai1;          % A trick to improve the Jacobian (and speed-up the solver)
kcCai2 = 10^3*cCai2;

if USPaT(t) == 0
Veff1 = f1Veff0(Q1);
Veff2 = 1000*Q2/Cm0;

Out1 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*RSI))*(Veff2-Veff1)-10^(-3)*(Gl*(Veff1-Vl)+Gna*m1.^3.*h1.*(Veff1-Vna)+Gk*n1.^4.*(Veff1-Vk)+...
   GT*p1.^2.*q1.*(Veff1-fVCa(cCai1))+GCa.*r1.^2.*(Veff1-Vk)+GA*a1.^2.*b1.*(Veff1-Vk)+...
   +GL*c1.^2.*d11.*d21.*(Veff1-fVCa(cCai1)));
cellfun(@(X) f1rt0.(['a_' X])(Q1)-f1rt0.(['apb_' X])(Q1)*rate1.(X),SONICgates);
(rinf(kcCai1)-r1)/taur(kcCai1)
(d2inf(kcCai1)-d21)/taud2(kcCai1);
-(0.001*(GT*p1.^2.*q1.*(Veff1-fVCa(cCai1))+...
GL*c1.^2.*d11.*d21.*(Veff1-fVCa(cCai1))))/(2*Far*(10236*10^(-9)))-cCai1/(tauCa)]; % Gain factor as in Kamaruvelu et al. (2016)

Out2 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*(1/fBLS-1)*RSI))*(Veff1-Veff2)-10^(-3)*(Gl*(Veff2-Vl)+Gna*m2.^3.*h2.*(Veff2-Vna)+Gk*n2.^4.*(Veff2-Vk)+...
   GT*p2.^2.*q2.*(Veff2-fVCa(cCai2))+GCa.*r2.^2.*(Veff2-Vk)+GA*a2.^2.*b2.*(Veff2-Vk)+...
   +GL*c2.^2.*d12.*d22.*(Veff2-fVCa(cCai2)));
cellfun(@(X) f1rtV.(['a_' X])(Veff2)-f1rtV.(['apb_' X])(Veff2)*rate2.(X),SONICgates);
(rinf(kcCai2)-r2)/taur(kcCai2)
(d2inf(kcCai2)-d22)/taud2(kcCai2);
-(0.001*(GT*p2.^2.*q2.*(Veff2-fVCa(cCai2))+...
GL*c2.^2.*d12.*d22.*(Veff2-fVCa(cCai2))))/(2*Far*(10236*10^(-9)))-cCai2/(tauCa)]; % Gain factor as in Kamaruvelu et al. (2016)

Out = [Out1;Out2];
else
Veff1 = f1VeffPa(Q1);
Veff2 = 1000*Q2/Cm0;

Out1 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*RSI))*(Veff2-Veff1)-10^(-3)*(Gl*(Veff1-Vl)+Gna*m1.^3.*h1.*(Veff1-Vna)+Gk*n1.^4.*(Veff1-Vk)+...
   GT*p1.^2.*q1.*(Veff1-fVCa(cCai1))+GCa.*r1.^2.*(Veff1-Vk)+GA*a1.^2.*b1.*(Veff1-Vk)+...
   +GL*c1.^2.*d11.*d21.*(Veff1-fVCa(cCai1)));
cellfun(@(X) f1rtPa.(['a_' X])(Q1)-f1rtPa.(['apb_' X])(Q1)*rate1.(X),SONICgates);
(rinf(kcCai1)-r1)/taur(kcCai1)
(d2inf(kcCai1)-d21)/taud2(kcCai1);
-(0.001*(GT*p1.^2.*q1.*(Veff1-fVCa(cCai1))+...
GL*c1.^2.*d11.*d21.*(Veff1-fVCa(cCai1))))/(2*Far*(10236*10^(-9)))-cCai1/(tauCa)]; % Gain factor as in Kamaruvelu et al. (2016)

Out2 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*(1/fBLS-1)*RSI))*(Veff1-Veff2)-10^(-3)*(Gl*(Veff2-Vl)+Gna*m2.^3.*h2.*(Veff2-Vna)+Gk*n2.^4.*(Veff2-Vk)+...
   GT*p2.^2.*q2.*(Veff2-fVCa(cCai2))+GCa.*r2.^2.*(Veff2-Vk)+GA*a2.^2.*b2.*(Veff2-Vk)+...
   +GL*c2.^2.*d12.*d22.*(Veff2-fVCa(cCai2)));
cellfun(@(X) f1rtV.(['a_' X])(Veff2)-f1rtV.(['apb_' X])(Veff2)*rate2.(X),SONICgates);
(rinf(kcCai2)-r2)/taur(kcCai2)
(d2inf(kcCai2)-d22)/taud2(kcCai2);
-(0.001*(GT*p2.^2.*q2.*(Veff2-fVCa(cCai2))+...
GL*c2.^2.*d12.*d22.*(Veff2-fVCa(cCai2))))/(2*Far*(10236*10^(-9)))-cCai2/(tauCa)]; % Gain factor as in Kamaruvelu et al. (2016)

Out = [Out1;Out2];
end
end