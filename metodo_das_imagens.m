syms theta x k
%metodo das imagens

xc = [-30*10^-2 0 30*10^-2];
yc = [1 1 1];

theta = 0;
R = 42.5*10^(-3);
V_fase = 69*10^3/sqrt(3);
I_fase = 500;

e0 = 8.854187*10^(-12);
e_isol = 2.85*e0;
e_isol_ext = 2.51*e0;

r1 = 19.5*10^(-3);
r2 = 37.75*10^(-3);
r3 = 37.97*10^(-3);
r4 = 42.5*10^(-3);

Pc = 1/(2*pi*e_isol)*log(r2/r1);
Ps = 1/(2*pi*e_isol_ext)*log(r4/r3);

Pcc = Pc + Ps;

P = [Pcc Ps; Ps Ps];

u = [V_fase*cos(theta); -V_fase*cos(theta)];
u2 = [V_fase*cos(theta-(2*pi)/3); -V_fase*cos(theta)];
u3 = [V_fase*cos(theta -(4*pi)/3); -V_fase*cos(theta)];

e1 = (2*pi*e0)*inv(P)*u;
e2 = (2*pi*e0)*inv(P)*u2;
e3 = (2*pi*e0)*inv(P)*u3;

q = [e1(1,1) ;e1(2,1)];
q2 = [e2(1,1); e2(2,1)];
q3 = [e3(1,1); e3(2,1)];


i =  [I_fase*cos(theta); I_fase*cos(theta)];
i2 = [I_fase*cos(theta-(2*pi)/3); I_fase*cos(theta-(2*pi)/3)];
i3 = [I_fase*cos(theta-(4*pi)/3); I_fase*cos(theta-(4*pi)/3)];

xci = [xc(1) xc(1)];

yci = [yc(1) -yc(1)];

xci_2 = [xc(2) xc(2)];

yci_2 = [yc(2) -yc(2)];

xci_3 = [xc(3) xc(3)];

yci_3 = [yc(3) -yc(3)];

dummy = 0;
Esolo = 0;
for k = 1:2
    dummy = (q(k)/(10*e0))*[(x- xci(k))/((x - xci(k))^2 +(-yci(k))^2), -yci(k)/((x - xci(k))^2 +(-yci(k))^2)]; 
    Esolo = Esolo + dummy; 
end

for k = 1:2
    dummy = (q2(k)/(10*e0))*[(x- xci_2(k))/((x - xci_2(k))^2 +(-yci_2(k))^2), -yci_2(k)/((x - xci_2(k))^2 +(-yci_2(k))^2)]; 
    Esolo = Esolo + dummy; 
end

for k = 1:2
    dummy = (q3(k)/(10*e0))*[(x- xci_3(k))/((x - xci_3(k))^2 +(-yci_3(k))^2), -yci_3(k)/((x - xci_3(k))^2 +(-yci_3(k))^2)]; 
    Esolo = Esolo + dummy; 
end
E_eficaz = sqrt(Esolo(1)^2+Esolo(2)^2);

dummy = 0;
Hsolo = 0;
for k = 1:2
    dummy = (i(k)/(2*pi))*[(x- xci(k))/((x - xci(k))^2 +(-yci(k))^2), -yci(k)/((x - xci(k))^2 +(-yci(k))^2)]; 
    Hsolo = Hsolo + dummy; 
end
for k = 1:2
    dummy = (i2(k)/(2*pi))*[(x- xci_2(k))/((x - xci_2(k))^2 +(-yci_2(k))^2), -yci_2(k)/((x - xci_2(k))^2 +(-yci_2(k))^2)]; 
    Hsolo = Hsolo + dummy; 
end
for k = 1:2
    dummy = (i3(k)/(2*pi))*[(x- xci_3(k))/((x - xci_3(k))^2 +(-yci_3(k))^2), -yci_3(k)/((x - xci_3(k))^2 +(-yci_3(k))^2)]; 
    Hsolo = Hsolo + dummy; 
end
H_eficaz = sqrt(Hsolo(1)^2+Hsolo(2)^2);


%plotar
k = 600;
a = linspace (-5, 5, k);
d = linspace (-40, 40 , k);
campEsol = [ - (9152857963535559*(a - 3/10))./(295147905179352825856*((a - 3/10).^2 + 1)) - (18305715927071123*(a + 3/10))./(590295810358705651712*((a + 3/10).^2 + 1)) - (73222863708284479*a)./(2361183241434822606848*(a.^2 + 1)), - 87408601626032737./(2361183241434822606848*(a.^2 + 1)) - 10926075203254089./(295147905179352825856.*((a - 3/10).^2 + 1)) - 32491453844819373./(590295810358705651712.*((a + 3/10).^2 + 1))];

campEfic = (((9152857963535559*(a - 3/10))./(295147905179352825856*((a - 3/10).^2 + 1)) + (18305715927071123*(a + 3/10))./(590295810358705651712*((a + 3/10).^2 + 1)) + (73222863708284479*a)./(2361183241434822606848*(a.^2 + 1))).^2 + (87408601626032737./(2361183241434822606848*(a.^2 + 1)) + 10926075203254089./(295147905179352825856*((a - 3/10).^2 + 1)) + 32491453844819373./(590295810358705651712*((a + 3/10).^2 + 1))).^2).^(1/2);
plot (d, campEfic)
xlabel ('distancia[m]')
ylabel ('Campo elétrico [kV/m]')
title ('Campo elétrico eficaz')