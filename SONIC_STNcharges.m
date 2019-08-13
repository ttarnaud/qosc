function Out = SimplNICESTNLUcharges(ESi,DISPLAY,tNICE,t,Q,...
    m,n,p,h,q,r,a,b,c,d1,d2,cCai,~,~,~,~,~,~,~,CmR,...
    Gna,Vna,Gk,Vk,Gl,Vl,GT,fVCa,GCa,GA,GL,minf,ninf,pinf,hinf,qinf,rinf,ainf,...
    binf,cinf,d1inf,d2inf,taum,taun,taup,tauh,tauq,taur,taua,taub,tauc,taud1,taud2,...
    Far,tauCa,VLIMs,kcCaiLIMs)
% Also calculates the charges QNa,QK,QT,QL,QKCa,QA,Ql in exchange for computational cost
% Note: it is not possible to obtain the charges in post-processing, due to
% the complicated Adam-Bashforth-Moulton iteration sheme that is used
if DISPLAY == 2
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));  
end
m = m*(m<=1)+(m>1); 
n = n*(n<=1)+(n>1);
p = p*(p<=1)+(p>1);
h = h*(h<=1)+(h>1);
q = q*(q<=1)+(q>1);
r = r*(r<=1)+(r>1);
a = a*(a<=1)+(a>1);
b = b*(b<=1)+(b>1);
c = c*(c<=1)+(c>1);
d1 = d1*(d1<=1)+(d1>1);
d2 = d2*(d2<=1)+(d2>1);

kcCai = 10^3*cCai;          % A trick to improve the Jacobian (and speed-up the solver)
V = 10^(3)*Q/CmR(t);
if (V<VLIMs(1)||V>VLIMs(2)), error('Error: transmembrane potential outside interpolation interval!'); end
if (kcCai<kcCaiLIMs(1)||kcCai>kcCaiLIMs(2)), error('Error: calcium concentration outside interpolation interval!'); end

Out = [ESi(t)-10^(-3)*(Gl*(V-Vl)+Gna*m.^3.*h.*(V-Vna)+Gk*n.^4.*(V-Vk)+...
   GT*p.^2.*q.*(V-fVCa(cCai))+GCa.*r.^2.*(V-Vk)+GA*a.^2.*b.*(V-Vk)+...
   +GL*c.^2.*d1.*d2.*(V-fVCa(cCai)));
(minf(V)-m)/taum(V);
(ninf(V)-n)/taun(V);
(pinf(V)-p)/taup(V);
(hinf(V)-h)/tauh(V);
(qinf(V)-q)/tauq(V);
(rinf(kcCai)-r)/taur(kcCai);
(ainf(V)-a)/taua(V);
(binf(V)-b)/taub(V);
(cinf(V)-c)/tauc(V);
(d1inf(V)-d1)/taud1(V);
(d2inf(kcCai)-d2)/taud2(kcCai);
-(0.001*(GT*p.^2.*q.*(V-fVCa(cCai))+...
GL*c.^2.*d1.*d2.*(V-fVCa(cCai))))/(2*Far*(10236*10^(-9)))-cCai/(tauCa);
-10^(-3)*(Gna*m.^3.*h.*(V-Vna));
-10^(-3)*(Gk*n.^4.*(V-Vk));
-10^(-3)*(GT*p.^2.*q.*(V-fVCa(cCai)));
-10^(-3)*(GL*c.^2.*d1.*d2.*(V-fVCa(cCai)));
-10^(-3)*(GCa.*r.^2.*(V-Vk));   
-10^(-3)*(GA*a.^2.*b.*(V-Vk));
-10^(-3)*(Gl*(V-Vl))]; % Gain factor as in Kamaruvelu et al. (2016)
end