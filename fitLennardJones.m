function [Arx,deltaxx,xx,yx,gof,output] = fitLennardJones(a,delta)
Ar = 10^5;				% Attraction/repulsion pressure coefficient (Pa) 
deltax = 1.4*10^(-9); 	% Initial gap between leaflets (no charge) (m) 
x = 5;					% Repulsion exponent
y = 3.3;				% Attraction exponent


R = @(Z) (a^2+Z.^2)/(2.*Z); % Radius of curvature [m]
z = @(r,Z) Funz(r,Z,R(Z)); % Curvature [m]
f = @(r,Z) Ar*((deltax./(2*z(r,Z)+delta)).^x-(deltax./(2*z(r,Z)+delta)).^y); % [Pa]
Pm = @(Z) (2./(Z.^2+a^2)).*integral(@(r)(r.*f(r,Z)),0,a); % Attraction/repulsion moulecular pressure [Pa]

ll = 0; rl = a; NN = 1e3;
PmRange = linspace(ll,rl,NN);
PmI = zeros(length(PmRange),1);
for iR = 1:length(PmRange)
PmI(iR) = Pm(PmRange(iR));
end

[status,errmsg] = license('checkout','Curve_Fitting_Toolbox'); %#ok<ASGLU>
if (status == 1)
fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0 0 0 0],'StartPoint',[Ar*10^(-5),deltax*10^(9),x,y]);
ft = fittype('Arxem5*10^(5)*((deltaxxe9*10^(-9)./(2*Z+delta)).^xx-(deltaxxe9*10^(-9)./(2*Z+delta)).^yx)','problem',{'delta'},'Independent','Z','options',fo);
[cfit,gof,output] = fit(PmRange',PmI,ft,'problem',delta,'display','off');
Arx = cfit.Arxem5*10^5;
xx = cfit.xx; yx = cfit.yx;
deltaxx = cfit.deltaxxe9*10^(-9);
else
LenJonesFit = @(X,Xdata) X(1)*10^(5)*((X(4)*10^(-9)./(2*Xdata+delta)).^X(2)-(X(4)*10^(-9)./(2*Xdata+delta)).^X(3));
X0 = [Ar*10^(-5),x,y,deltax*10^(9)]; xdata = PmRange'; ydata = PmI; lb = [0 0 0 0]; ub = [];
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','Display','off');
[X,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(LenJonesFit,X0,xdata,ydata,lb,ub,options); %#ok<ASGLU>
Arx = X(1)*10^5; xx = X(2); yx = X(3); deltaxx = X(4)*10^(-9);
gof = resnorm;
end
end