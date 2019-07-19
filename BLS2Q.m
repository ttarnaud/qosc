function Out = BLS2Q(DISPLAY,tBLS,t,Z,dZ,na,Ca,R,rhol,PecQ,Q0,Pin,Pm,Po,PaT,omega,PS,delta0,mus,mul,S,Da,ka,A,B,deltaR)
if DISPLAY == 2
global reverseStr; %#ok<TLEV>
Progress = 100*(t-tBLS(1))/(tBLS(end)-tBLS(1));  %#ok<*NASGU>
msg = sprintf('Progress: %3.1f', Progress); 
fprintf([reverseStr, msg]);
reverseStr = repmat(sprintf('\b'), 1, length(msg));  
end
Out = [dZ;
(-3/(2*R(Z)))*dZ^2+1/(rhol*abs(R(Z)))*(Pin(na,Z)+PecQ(Q0,Z)+...
Pm(Z)-Po+PaT(t)*sin(omega*t)-PS(Z)-(4/abs(R(Z)))*dZ*...
    ((3*delta0*mus)/abs(R(Z))+mul));
    S(Z)*Da*(Ca(1)-Pin(na,Z)/ka)/deltaR;
    A*Ca+B(Pin(na,Z))];
end