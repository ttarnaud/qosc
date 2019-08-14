c = 1515;				% Speed of sound surrounding medium (m/s)
rhol = 1028;			% Density surrounding medium (kg/m^3)
Pa2I = @(Pa) Pa^2/(2*rhol*c);
PaR = 1e3*logspace(log10(50),log10(150),4);

for i = 1:length(PaR)
SONICrun('0.300','2','0.050','0.150','0.5e6','1','0',num2str(Pa2I(PaR(i))),'0','0','1','0','0','1','1','0','0','0')
end