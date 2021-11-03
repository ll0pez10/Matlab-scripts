syms lambda gama y x 
%metodo alternativo de calculo do campo eletrico

u0 = 4*pi*10^-7;
e0 = 8.854187*10^(-12);
e1 = 10;
w = 2*pi*60;
R = 42.5*10^(-3);
h = 1;
J = 500/(pi*R^2);
d = sqrt((h-y)^2 +x^2);
D = sqrt((h+y)^2 +x^2);
sigma1 = 10^-3;
e1 = 10;
sigma2 = 10^-15;

gama1 = sqrt(1i*w*u0*(sigma1 + 1i*w*e1*e0));
gama2 = sqrt(1i*w*u0*(sigma2 + 1i*w*e0));
eta1 = gama1^2- gama^2;
eta2 = gama2^2- gama^2;


u1 = sqrt(lambda^2 + gama1^2 - gama^2);
u2 = sqrt(lambda^2 + gama2^2 - gama^2);
n = gama2/gama1;

f1 = ((exp(-u1*(h+y)))/(u1 + u2))*cos(x*lambda);
f2 = ((exp(-u1*(h+y)))/(n*(u1 + u2)))*cos(x*lambda);

E1 = -((1i*w*u0*J)/(2*pi))*(besselk(0,eta1*d) - besselk(0,eta1*D) + 2*int(f1,lambda, [0 inf]) - (gama^2/gama1^2)*(besselk(0,eta1*d) - besselk(0,eta1*D) + 2*int(f2,lambda ,[0 inf]) ));