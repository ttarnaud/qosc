function Out = SONIC_GPe_nanoMC(ESi,USPaT,DISPLAY,tNICE,t,Q1,h1,n1,r1,CA1,Q2,h2,n2,r2,CA2,Gna,Vna,Gk,Vk,Gl,Vl,...
    GT,VT,GCa,VCa,Gahp,minf,ainf,sinf,f1Veff0,f1VeffPa,f1rt0,f1rtPa,f1rtV,SONICgates,Cm0,aBLS,fBLS,RSI)
if DISPLAY == 1
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));  
end
rate1 = struct; rate2 = struct;
n1 = n1*(n1<=1&n1>=0)+(n1>1); n2 = n2*(n2<=1&n2>=0)+(n2>1);
h1 = h1*(h1<=1&h1>=0)+(h1>1); h2 = h2*(h2<=1&h2>=0)+(h2>1);
r1 = r1*(r1<=1&r1>=0)+(r1>1); r2 = r2*(r2<=1&r2>=0)+(r2>1);
rate1.('n') = n1; rate1.('h') = h1; rate1.('r') = r1;
rate2.('n') = n2; rate2.('h') = h2; rate2.('r') = r2;

if USPaT(t) == 0
Veff1 = f1Veff0(Q1);
Veff2 = 1000*Q2/Cm0;

Out1 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*RSI))*(Veff2-Veff1)-10^(-3)*(Gl*(Veff1-Vl)+Gna*minf(Veff1).^3.*h1.*(Veff1-Vna)+...
   Gk*n1.^4.*(Veff1-Vk)+GT*ainf(Veff1).^3.*r1.*(Veff1-VT)+...
   GCa.*sinf(Veff1).^2.*(Veff1-VCa)+Gahp*(Veff1-Vk).*(CA1./(CA1+10)));
cellfun(@(X) f1rt0.(['a_' X])(Q1)-f1rt0.(['apb_' X])(Q1)*rate1.(X),SONICgates)
-10^(-4)*(0.1*GT*ainf(Veff1).^3.*r1.*(Veff1-VT)+...
   0.1*GCa.*sinf(Veff1).^2.*(Veff1-VCa)+15*10^(3)*CA1)];
Out2 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*(1/fBLS-1)*RSI))*(Veff1-Veff2)-10^(-3)*(Gl*(Veff2-Vl)+Gna*minf(Veff2).^3.*h2.*(Veff2-Vna)+...
   Gk*n2.^4.*(Veff2-Vk)+GT*ainf(Veff2).^3.*r2.*(Veff2-VT)+...
   GCa.*sinf(Veff2).^2.*(Veff2-VCa)+Gahp*(Veff2-Vk).*(CA2./(CA2+10)));
cellfun(@(X) f1rtV.(['a_' X])(Veff2)-f1rtV.(['apb_' X])(Veff2)*rate2.(X),SONICgates)
-10^(-4)*(0.1*GT*ainf(Veff2).^3.*r.*(Veff2-VT)+...
   0.1*GCa.*sinf(Veff2).^2.*(Veff2-VCa)+15*10^(3)*CA2)];

Out = [Out1;Out2];
else
Veff1 = f1VeffPa(Q1);
Veff2 = 1000*Q2/Cm0;

Out1 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*RSI))*(Veff2-Veff1)-10^(-3)*(Gl*(Veff1-Vl)+Gna*minf(Veff1).^3.*h1.*(Veff1-Vna)+...
   Gk*n1.^4.*(Veff1-Vk)+GT*ainf(Veff1).^3.*r1.*(Veff1-VT)+...
   GCa.*sinf(Veff1).^2.*(Veff1-VCa)+Gahp*(Veff1-Vk).*(CA1./(CA1+10)));
cellfun(@(X) f1rtPa.(['a_' X])(Q1)-f1rtPa.(['apb_' X])(Q1)*rate1.(X),SONICgates)
-10^(-4)*(0.1*GT*ainf(Veff1).^3.*r1.*(Veff1-VT)+...
   0.1*GCa.*sinf(Veff1).^2.*(Veff1-VCa)+15*10^(3)*CA1)];
Out2 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*(1/fBLS-1)*RSI))*(Veff1-Veff2)-10^(-3)*(Gl*(Veff2-Vl)+Gna*minf(Veff2).^3.*h2.*(Veff2-Vna)+...
   Gk*n2.^4.*(Veff2-Vk)+GT*ainf(Veff2).^3.*r2.*(Veff2-VT)+...
   GCa.*sinf(Veff2).^2.*(Veff2-VCa)+Gahp*(Veff2-Vk).*(CA2./(CA2+10)));
cellfun(@(X) f1rtV.(['a_' X])(Veff2)-f1rtV.(['apb_' X])(Veff2)*rate2.(X),SONICgates)
-10^(-4)*(0.1*GT*ainf(Veff2).^3.*r.*(Veff2-VT)+...
   0.1*GCa.*sinf(Veff2).^2.*(Veff2-VCa)+15*10^(3)*CA2)];

Out = [Out1;Out2];
end
end