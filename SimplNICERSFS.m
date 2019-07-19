function Out = SimplNICERSFS(ESi,DISPLAY,tNICE,t,Q,m,n,p,h,CmR,Gna,Vna,Gk,Vk,Gm,Gl,Vl,am,bm,an,bn,pinf,taup,ah,bh)
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

V = 10^(3)*Q/CmR(t);
Out = [ESi(t)-10^(-3)*(Gna*m^3*h*(V-Vna)+Gk*n^4*(V-Vk)+Gm*p*(V-Vk)+Gl*(V-Vl));
am(V)-(am(V)+bm(V))*m;
an(V)-(an(V)+bn(V))*n;
(pinf(V)-p)/taup(V);
ah(V)-(ah(V)+bh(V))*h];
end