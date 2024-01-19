function Out = FunPm(Z,delta,Pm)
if Z<=-delta/2
error('FunPm:UNPHYS','Unphysical solution: Z<=-delta/2!');
else
Out = Pm(Z);
end
end