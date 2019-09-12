function Out = SONIC_RE_nanoMC(ESi,USPaT,DISPLAY,tNICE,t,Q1,h1,m1,n1,s1,u1,cCai1,Q2,h2,m2,n2,s2,u2,cCai2,...
    Gna,Vna,Gk,Vk,GT,fVCa,Gl,Vl,Far,deffCa,tauCa,f1Veff0,f1VeffPa,f1rt0,f1rtPa,f1rtV,SONICgates,Cm0,aBLS,fBLS,RSI)
if DISPLAY == 1
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg)); 
end

rate1 = struct; rate2 = struct;
m1 = m1*(m1<=1&m1>=0)+(m1>1); m2 = m2*(m2<=1&m2>=0)+(m2>1);
n1 =  n1*(n1<=1&n1>=0)+(n1>1); n2 =  n2*(n2<=1&n2>=0)+(n2>1);
h1 = h1*(h1<=1&h1>=0)+(h1>1); h2 = h2*(h2<=1&h2>=0)+(h2>1);
s1 = s1*(s1<=1&s1>=0)+(s1>1); s2 = s2*(s2<=1&s2>=0)+(s2>1);
u1 = u1*(u1<=1&u1>=0)+(u1>1); u2 = u2*(u2<=1&u2>=0)+(u2>1);
rate1.('m') = m1; rate1.('n') = n1; rate1.('h') = h1;
rate1.('s') = s1; rate1.('u') = u1;
rate2.('m') = m2; rate2.('n') = n2; rate2.('h') = h2;
rate2.('s') = s2; rate2.('u') = u2;

if USPaT(t) == 0
Veff1 = f1Veff0(Q1);
Veff2 = 1000*Q2/Cm0;

Out1 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*RSI))*(Veff2-Veff1)-10^(-3)*(Gna*m1^3*h1*(Veff1-Vna)+Gk*n1^4*(Veff1-Vk)+GT*s1^2*u1*(Veff1-fVCa(cCai1))+Gl*(Veff1-Vl));
cellfun(@(X) f1rt0.(['a_' X])(Q1)-f1rt0.(['apb_' X])(Q1)*rate1.(X),SONICgates);
-(10^(-3)*GT*s1^2*u1*(Veff1-fVCa(cCai1)))/(2*Far*deffCa)-cCai1/tauCa];
Out2 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*(1/fBLS-1)*RSI))*(Veff1-Veff2)-10^(-3)*(Gna*m2^3*h2*(Veff2-Vna)+Gk*n2^4*(Veff2-Vk)+GT*s2^2*u2*(Veff2-fVCa(cCai2))+Gl*(Veff2-Vl));
cellfun(@(X) f1rtV.(['a_' X])(Veff2)-f1rtV.(['apb_' X])(Veff2)*rate2.(X),SONICgates);
-(10^(-3)*GT*s2^2*u2*(Veff2-fVCa(cCai2)))/(2*Far*deffCa)-cCai2/tauCa];

Out = [Out1;Out2];
else
Veff1 = f1VeffPa(Q1);
Veff2 = 1000*Q2/Cm0;

Out1 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*RSI))*(Veff2-Veff1)-10^(-3)*(Gna*m1^3*h1*(Veff1-Vna)+Gk*n1^4*(Veff1-Vk)+GT*s1^2*u1*(Veff1-fVCa(cCai1))+Gl*(Veff1-Vl));
cellfun(@(X) f1rtPa.(['a_' X])(Q1)-f1rtPa.(['apb_' X])(Q1)*rate1.(X),SONICgates);
-(10^(-3)*GT*s1^2*u1*(Veff1-fVCa(cCai1)))/(2*Far*deffCa)-cCai1/tauCa];
Out2 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*(1/fBLS-1)*RSI))*(Veff1-Veff2)-10^(-3)*(Gna*m2^3*h2*(Veff2-Vna)+Gk*n2^4*(Veff2-Vk)+GT*s2^2*u2*(Veff2-fVCa(cCai2))+Gl*(Veff2-Vl));
cellfun(@(X) f1rtV.(['a_' X])(Veff2)-f1rtV.(['apb_' X])(Veff2)*rate2.(X),SONICgates);
-(10^(-3)*GT*s2^2*u2*(Veff2-fVCa(cCai2)))/(2*Far*deffCa)-cCai2/tauCa];

Out = [Out1;Out2];
end
end