function Out = BLS1Q(DISPLAY,tBLS,t,Z,dZ,na,R,rhol,PecQ,Q0,DeltaQ,Pin,Pm,Po,PaT,omega,PS,delta0,mus,mul,S,Da,Ci,ka,ksi) %#ok<INUSL>
if DISPLAY==2
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tBLS(1))/(tBLS(end)-tBLS(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg)); 
end
Q = Q0+DeltaQ*cos(omega*t);

Out = [dZ;
(-3/(2*R(Z)))*dZ^2+1/(rhol*abs(R(Z)))*(PecQ(Q,Z)+Pin(na,Z)+...
Pm(Z)-Po+PaT(t)*sin(omega*t)-PS(Z)-(4/abs(R(Z)))*dZ*...
    ((3*delta0*mus)/abs(R(Z))+mul));
    2*S(Z)*Da*(Ci-Pin(na,Z)/ka)/ksi];
end