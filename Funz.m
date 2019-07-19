function Out = Funz(r,Z,R)
if Z == 0
Out=0;
else
Out = sign(R)*sqrt(R^2-r.^2)-R+Z;
end
end