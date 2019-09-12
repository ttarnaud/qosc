function Out = SONIC_HH_nanoMC(ESi,USPaT,DISPLAY,tNICE,t,Q1,h1,m1,n1,Q2,h2,m2,n2,Gna,Vna,Gk,Vk,Gl,Vl,f1Veff0,f1VeffPa,f1rt0,f1rtPa,f1rtV,SONICgates,Cm0,aBLS,fBLS,RSI)
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
rate1.('m') = m1; rate1.('n') = n1; rate1.('h') = h1;
rate2.('m') = m1; rate2.('n') = n1; rate2.('h') = h2;

if USPaT(t) == 0
Veff1 = f1Veff0(Q1);
Veff2 = 1000*Q2/Cm0;

Out = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*RSI))*(Veff2-Veff1)-10^(-3)*(Gl*(Veff-Vl)+Gna*m1.^3.*h1.*(Veff1-Vna)+...
    Gk*n1.^4.*(Veff1-Vk));
    cellfun(@(X) f1rt0.(['a_' X])(Q1)-f1rt0.(['apb_' X])(Q1)*rate1.(X),SONICgates)];
Out = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*(1/fBLS-1)*RSI))*(Veff1-Veff2)-10^(-3)*(Gl*(Veff2-Vl)+Gna*m2.^3.*h2.*(Veff2-Vna)+...
    Gk*n2.^4.*(Veff2-Vk));
    cellfun(@(X) f1rtV.(['a_' X])(Veff2)-f1rtV.(['apb_' X])(Veff2)*rate2.(X),SONICgates)];

Out = [Out1;Out2];
else
Veff1 = f1VeffPa(Q1);
Veff2 = 1000*Q2/Cm0;

Out = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*RSI))*(Veff2-Veff1)-10^(-3)*(Gl*(Veff-Vl)+Gna*m1.^3.*h1.*(Veff1-Vna)+...
    Gk*n1.^4.*(Veff1-Vk));
    cellfun(@(X) f1rtPa.(['a_' X])(Q1)-f1rtPa.(['apb_' X])(Q1)*rate1.(X),SONICgates)];
Out = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*(1/fBLS-1)*RSI))*(Veff1-Veff2)-10^(-3)*(Gl*(Veff2-Vl)+Gna*m2.^3.*h2.*(Veff2-Vna)+...
    Gk*n2.^4.*(Veff2-Vk));
    cellfun(@(X) f1rtV.(['a_' X])(Veff2)-f1rtV.(['apb_' X])(Veff2)*rate2.(X),SONICgates)];

Out = [Out1;Out2];
end
end