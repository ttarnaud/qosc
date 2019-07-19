function Out = SimplNICELTSLU(ESi,DISPLAY,tNICE,t,Q,m,n,p,h,s,u,CmR,Gna,Vna,Gk,Vk,Gm,GT,VCa,Gl,Vl,am,ampbm,an,anpbn,pinf,taup,ah,ahpbh,sinf,taus,uinf,tauu,VLIMs)
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
s = s*(s<=1)+(s>1);
u = u*(u<=1)+(u>1);
V = 10^(3)*Q/CmR(t);
if (V<VLIMs(1)||V>VLIMs(2)), error('Error: transmembrane potential outside interpolation interval!'); end
Out = [ESi(t)-10^(-3)*(Gna*m^3*h*(V-Vna)+Gk*n^4*(V-Vk)+Gm*p*(V-Vk)+GT*s^2*u*(V-VCa)+Gl*(V-Vl));
am(V)-(ampbm(V))*m;
an(V)-(anpbn(V))*n;
(pinf(V)-p)/taup(V);
ah(V)-(ahpbh(V))*h;
(sinf(V)-s)/taus(V);
(uinf(V)-u)/tauu(V)];
end