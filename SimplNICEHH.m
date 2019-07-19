function Out = SimplNICEHH(ESi,DISPLAY,tNICE,t,Q,m,n,h,CmR,Gna,Vna,Gk,Vk,Gl,Vl,...
    am,bm,an,bn,ah,bh)
if DISPLAY == 2
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tNICE(1))/(tNICE(2)-tNICE(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
m = m*(m<=1)+(m>1); 
n = n*(n<=1)+(n>1);
h = h*(h<=1)+(h>1);

V = 10^(3)*Q/CmR(t);
Out = [ESi(t)-10^(-3)*(Gl*(V-Vl)+Gna*m.^3.*h.*(V-Vna)+...
    Gk*n.^4.*(V-Vk));
am(V)-(am(V)+bm(V))*m;
an(V)-(an(V)+bn(V))*n;
ah(V)-(ah(V)+bh(V))*h];
end