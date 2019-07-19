function Out = FunCm(Z,Cm0,a,delta)
temp = (((Cm0*delta)/a^2)*(Z+((a^2-Z.^2-Z*delta)./(2*Z))...
    .*log((2*Z+delta)/delta)));
temp(isnan(temp)) = Cm0;
Out = temp;
end