syms theta x k s
%metodo das imagens

%resistividade do nucleo central
pc = 3.365*10^(-8);
%resistividade da blindagem metalica
ps = 1.718*10^(-8);

w = 2*pi*60;
theta = 0;
u0 = 4*pi*10^(-7);
e0 = 8.854187*10^(-12); 
r1 = 19.5*10^(-3);
r2 = 37.75*10^(-3);
r3 = 37.97*10^(-3);
r4 = 42.5*10^(-3);
I_fase = 500;
xc = [0 0 0];
yc = [50*10^(-2) 70*10^(-2) 90*10^(-2)];
R = 42.5*10^(-3);
V_fase = 300*10^3/sqrt(3);

%coeficiente dieletrico da camada interna de isolante
e_isol = 2.85*e0;
%coeficiente dieletrico da camada externa de isolante
e_isol_ext = 2.51*e0;

%potenciais das camadas isolantes
Pc =1/(2*pi*e_isol)*log(r2/r1);
Ps =1/(2*pi*e_isol_ext)*log(r4/r3);

Pcc = Pc + Ps;

%matriz de coeficientes de potencial
P = [Pcc Ps;Ps Ps];

u = [V_fase*cos(theta); 0];
u2 = [V_fase*cos(theta-(2*pi)/3); 0];
u3 = [V_fase*cos(theta -(4*pi)/3); 0];

e1 = (1/R)*inv(P)*u;
e2 = (1/R)*inv(P)*u2;
e3 = (1/R)*inv(P)*u3;

q = [e1(1,1) ;e1(2,1)];
q2 = [e2(1,1); e2(2,1)];
q3 = [e3(1,1); e3(2,1)];

i = [ I_fase*cos(theta); -I_fase*cos(theta)];
i2 = [ I_fase*cos(theta-(2*pi)/3); -I_fase*cos(theta-(2*pi)/3)];
i3 = [ I_fase*cos(theta-(4*pi)/3); -I_fase*cos(theta-(4*pi)/3)];

xci = [xc(1) xc(1)];

yci = [yc(1) -yc(1)];

xci_2 = [xc(2) xc(2)];

yci_2 = [yc(2) -yc(2)];

xci_3 = [xc(3) xc(3)];

yci_3 = [yc(3) -yc(3)];

dummy = 0;
Esolo = 0;
for k = 1:2
    dummy = q(k)*[(x- xci(k))/((x - xci(k))^2 +(-yci(k))^2) -yci(k)/((x - xci(k))^2 +(-yci(k))^2)]; 
    Esolo = Esolo + dummy; 
end

for k = 1:2
    dummy = q2(k)*[(x- xci_2(k))/((x - xci_2(k))^2 +(-yci_2(k))^2) -yci_2(k)/((x - xci_2(k))^2 +(-yci_2(k))^2)]; 
    Esolo = Esolo + dummy; 
end

for k = 1:2
    dummy = q3(k)*[(x- xci_3(k))/((x - xci_3(k))^2 +(-yci_3(k))^2) -yci_3(k)/((x - xci_3(k))^2 +(-yci_3(k))^2)]; 
    Esolo = Esolo + dummy; 
end
E_eficaz = sqrt(Esolo(1).^2+Esolo(2).^2);


dummy = 0;
Hsolo = 0;
for k = 1:2
    dummy = i(k)*[(x- xci(k))/((x - xci(k))^2 +(-yci(k))^2) -yci(k)/((x - xci(k))^2 +(-yci(k))^2)]; 
    Hsolo = Hsolo + dummy; 
end
for k = 1:2
    dummy = i2(k)*[(x- xci_2(k))/((x - xci_2(k))^2 +(-yci_2(k))^2) -yci_2(k)/((x - xci_2(k))^2 +(-yci_2(k))^2)]; 
    Hsolo = Hsolo + dummy; 
end
for k = 1:2
    dummy = i3(k)*[(x- xci_3(k))/((x - xci_3(k))^2 +(-yci_3(k))^2) -yci_3(k)/((x - xci_3(k))^2 +(-yci_3(k))^2)]; 
    Hsolo = Hsolo + dummy; 
end
H_eficaz = sqrt(Hsolo(1)^2+Hsolo(2)^2);


%campo externo ao condutor (no solo)
sigma_solo = 10^(-3);
e_solo = 10;

%rô do solo
pg = 1/sigma_solo;
        
%cte de euler
gama = 0.577215665;
Ng = sqrt(1i*w*u0/pg);
%d = sqrt((hi - hk)^2 + x^2);
d = 30*10^-2;
h = 1;
R = 42.5*10^-3;

l = 2*h;          %soma das profunidades
        
Z11 = (1i*w*u0/(2*pi))*(-log((gama*Ng*R)/2)+0.5-(4/3)*Ng*h);
Z22 = (1i*w*u0/(2*pi))*(-log((gama*Ng*R)/2)+0.5-(4/3)*Ng*h);
Z33 = (1i*w*u0/(2*pi))*(-log((gama*Ng*R)/2)+0.5-(4/3)*Ng*h);
Z12 = (1i*w*u0/(2*pi))*(-log((gama*Ng*d)/2)+0.5-(2/3)*Ng*l);
Z13 = (1i*w*u0/(2*pi))*(-log((gama*Ng*2*d)/2)+0.5-(2/3)*Ng*l);
Z23 = (1i*w*u0/(2*pi))*(-log((gama*Ng*d)/2)+0.5-(2/3)*Ng*l);
Z21 = Z12;
Z31 = Z13;
Z32 = Z23;

        
Zew = [Z11 Z12 Z13;Z21 Z22 Z23;Z31 Z32 Z33]
        
I =[500*cos(theta); -500*cos(theta-(2*pi)/3); 500*cos(theta-(4*pi)/3)]
        
Eext = -Zew*I

%campo na superficie do condutor

Nc = sqrt(1i*w*u0/pc);
Ns = sqrt(1i*w*u0/ps);

Z1 = (pc*Nc/(2*pi*r1))*coth(0.777*Nc*r1) + 0.356*pc/(pi*r1^2);
Z2 = (1i*w*u0/(2*pi))*log(r2/r1);
Z6 = (1i*w*u0/(2*pi))*log(r4/r3);
Z3 = ((ps*Ns)/(2*pi*r2))*coth(Ns*(r3-r2)) + ps/(2*pi*r2*(r2+r3));
Z4 = ((ps*Ns)/(pi*(r2+r3)))*csch(Ns*(r3-r2));
Z5 = ((ps*Ns)/(2*pi*r3))*coth(Ns*(r3-r2)) + ps/(2*pi*r3*(r2+r3));

Zcabo = [Z1+Z2+Z3+Z5+Z6-2*Z4 Z5+Z6-Z4; Z5+Z6-Z4 Z5+Z6];

Eint = Zcabo*i

k = 600;
a = linspace (-10, 10, k);
d = linspace (-40, 40 , k);
cEsol = [ 0, 2905167679615873./(18446744073709551616*(a.^2 + 49/100)) - 8300479084616783./(36893488147419103232*(a.^2 + 1/4)) + 7470431176155111./(36893488147419103232*(a.^2 + 81/100))];

cEfic = (((9152857963535559*(a - 3/10))./(295147905179352825856*((a - 3/10).^2 + 1)) + (18305715927071123*(a + 3/10))./(590295810358705651712*((a + 3/10).^2 + 1)) + (73222863708284479*a)./(2361183241434822606848*(a.^2 + 1))).^2 + (87408601626032737./(2361183241434822606848*(a.^2 + 1)) + 10926075203254089./(295147905179352825856*((a - 3/10).^2 + 1)) + 32491453844819373./(590295810358705651712*((a + 3/10).^2 + 1))).^2).^(1/2);
plot (d, cEfic)
xlabel ('distancia[m]')
ylabel ('Campo elétrico [kV/m]')
title ('Campo elétrico eficaz')