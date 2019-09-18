function Out = SONIC_ThRT_nanoMC(ESi,USPaT,DISPLAY,tNICE,t,Q1,h1,r1,Q2,h2,r2,Gna,Vna,Gk,Vk,Gl,Vl,...
    GT,VT,minf,pinf,f1Veff0,f1VeffPa,f1rt0,f1rtPa,f1rtV,SONICgates,Cm0,aBLS,fBLS,RSI,proteinMode)
if DISPLAY == 1
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));  
end
switch proteinMode
    case 0, xP = 1; xl = 1;             % (ratio of protein coverage (electrolytes, leak) in the BLS compartment)
    case 1, xP = 0; xl = 1;
    case 2, xP = 0; xl = 0;
end

rate1 = struct; rate2 = struct;
h1 = h1*(h1<=1&h1>=0)+(h1>1); h2 = h2*(h2<=1&h2>=0)+(h2>1);
r1 = r1*(r1<=1&r1>=0)+(r1>1); r2 = r2*(r2<=1&r2>=0)+(r2>1);
rate1.('h') = h1; rate1.('r') = r1;
rate2.('h') = h2; rate2.('r') = r2;

if USPaT(t) == 0
Veff1 = f1Veff0(Q1);
Veff2 = 1000*Q2/Cm0;

Out1 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*RSI))*(Veff2-Veff1)-10^(-3)*(xl*Gl*(Veff1-Vl)+xP*Gna*minf(Veff1).^3.*h1.*(Veff1-Vna)+...
    xP*Gk*(0.75*(1-h1)).^4.*(Veff1-Vk)+xP*GT*pinf(Veff1).^2.*r1.*(Veff1-VT));
cellfun(@(X) f1rt0.(['a_' X])(Q1)-f1rt0.(['apb_' X])(Q1)*rate1.(X),SONICgates)];

Out2 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*(1/fBLS-1)*RSI))*(Veff1-Veff2)-10^(-3)*(Gl*(Veff2-Vl)+Gna*minf(Veff2).^3.*h2.*(Veff2-Vna)+...
    Gk*(0.75*(1-h2)).^4.*(Veff2-Vk)+GT*pinf(Veff2).^2.*r2.*(Veff2-VT));
cellfun(@(X) f1rtV.(['a_' X])(Veff2)-f1rtV.(['apb_' X])(Veff2)*rate2.(X),SONICgates)];

Out = [Out1;Out2];
else
Veff1 = f1VeffPa(Q1);
Veff2 = 1000*Q2/Cm0;

Out1 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*RSI))*(Veff2-Veff1)-10^(-3)*(xl*Gl*(Veff1-Vl)+xP*Gna*minf(Veff1).^3.*h1.*(Veff1-Vna)+...
    xP*Gk*(0.75*(1-h1)).^4.*(Veff1-Vk)+xP*GT*pinf(Veff1).^2.*r1.*(Veff1-VT));
cellfun(@(X) f1rtPa.(['a_' X])(Q1)-f1rtPa.(['apb_' X])(Q1)*rate1.(X),SONICgates)];

Out2 = [ESi(t)+10^(-3)*(1/(pi*aBLS^2*(1/fBLS-1)*RSI))*(Veff1-Veff2)-10^(-3)*(Gl*(Veff2-Vl)+Gna*minf(Veff2).^3.*h2.*(Veff2-Vna)+...
    Gk*(0.75*(1-h2)).^4.*(Veff2-Vk)+GT*pinf(Veff2).^2.*r2.*(Veff2-VT));
cellfun(@(X) f1rtV.(['a_' X])(Veff2)-f1rtV.(['apb_' X])(Veff2)*rate2.(X),SONICgates)];

Out = [Out1;Out2];
end
end