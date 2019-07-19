function Out = SimplNICEGPe(ESi,DISPLAY,tNICE,t,Q,n,h,r,CA,CmR,Gna,Vna,Gk,Vk,Gl,Vl,...
    GT,VT,GCa,VCa,Gahp,minf,ainf,sinf,ninf,hinf,rinf,taun,tauh,taur)
if DISPLAY == 2
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));  
end
n = n*(n<=1)+(n>1);
h = h*(h<=1)+(h>1);
r = r*(r<=1)+(r>1);

V = 10^(3)*Q/CmR(t);
Out = [ESi(t)-10^(-3)*(Gl*(V-Vl)+Gna*minf(V).^3.*h.*(V-Vna)+...
   Gk*n.^4.*(V-Vk)+GT*ainf(V).^3.*r.*(V-VT)+...
   GCa.*sinf(V).^2.*(V-VCa)+Gahp*(V-Vk).*(CA./(CA+10)));
(ninf(V)-n)/taun(V);
(hinf(V)-h)/tauh(V);
(rinf(V)-r)/taur(V);
-10^(-4)*(0.1*GT*ainf(V).^3.*r.*(V-VT)+...
   0.1*GCa.*sinf(V).^2.*(V-VCa)+15*10^(3)*CA)];
end