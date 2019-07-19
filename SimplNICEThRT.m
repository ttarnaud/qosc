function Out = SimplNICEThRT(ESi,DISPLAY,tNICE,t,Q,h,r,CmR,Gna,Vna,Gk,Vk,Gl,Vl,...
    GT,VT,hinf,tauh,rinf,taur,minf,pinf)
if DISPLAY == 2
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));  
end
h = h*(h<=1)+(h>1);
r = r*(r<=1)+(r>1);

V = 10^(3)*Q/CmR(t);
Out = [ESi(t)-10^(-3)*(Gl*(V-Vl)+Gna*minf(V).^3.*h.*(V-Vna)+...
    Gk*(0.75*(1-h)).^4.*(V-Vk)+GT*pinf(V).^2.*r.*(V-VT));
(hinf(V)-h)/tauh(V);
(rinf(V)-r)/taur(V)];
end