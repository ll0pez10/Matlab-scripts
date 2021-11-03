syms s L C V(z,t) I(z,t) w G l R R1

%os valores de impedancia e indutancia são pelo comprimento
%indutancia do condutor
L = (u/2*pi)*ln(b/a)
%resistencia do condutor
R = (1/(2*pi*d*sigma_cond)*((1/a) + (1/b))
%capacitancia do condutor
C = (2*pi*e)/ln(b/a)
%condutancia do condutor
G = (2*pi*sigma)/ln(b/a)

%condutância
G = 1/R;

%velocidade da onda caso sem perdas
s = 1/(L*C);

%equacao telegráfica da tensão
pde1 = diff(V, z) == -(R * I + L*diff(I, t));

%equacao telegrafica da corrente
pde2 = diff(I, z) == -(G*V + C*diff(V, t));

%eq de onda geral p/ L.T da tensao
pde3 = diff(V,z,2) == L*C*diff(V,t,2) + (L*G + R*C)*diff(V,t) + R*G*V;

%eq de onda geral p/ L.T da corrente
pde4 = diff(I,z,2) == L*C*diff(I,t,2) + (L*G + R*C)*diff(I,t) + R*G*I;

%condicao de propagacao sem perdas
%R = 0;
%G = 0;

%impedancia caracteristica caso sem perdas
Z = (L/C)^(1/2);

%admitancia paralela
Yp = G + wCi;

%impedancia série equivalente da rede
Zs = R + wLi;

%condicao de heavside
R*C = G * L;

%coeficientes alpha + jbeta
alpha = (1/2)*(R*(C/L)^(1/2) + G*(L/C)^(1/2));

beta = w*(L/C)^(1/2)*(1 + (1/8)*( G/(w*C) - R/(w*L))^2);

%matriz R
L = [Vr; Ir];
%matriz S
S = [Vs; Is];

sigma = w*(L*C)^(1/2)i;

%n lembro o que representam
Vs = Vr*cos(theta) + Z*Ir*sin(theta);
Is = Ir*cos(theta) + (Vr/Z)*sin(theta)i;

L = [cosh(sigma*l) -Z*sinh(sigma*l); -(1/Z)*sinh(sigma*l) cosh(sigma*l)]*S;
S = [cosh(sigma*l) Z*sinh(sigma*l); (1/Z)*sinh(sigma*l) cosh(sigma*l)]*L;


%variaveis
h = 20;

%parametros iniciais
e0 = 8.854*10^(-12);
u0 = 1.257*10^(-6);
sigma_central = inf;
sigma_isolante = 0;
sigma_coroa = inf;
e_isolante = 3;

%primeiro item
R2 = 30*10^(-3);
R3 = 32*10^(-3);
I = 150;
V = 30*10^3;
    %formula vetor de poynting
    P = real(V* conj(I));
    %formula da densidade de energia elétrica
    
    %formula da desidade de energia magnetica
    
    %razao que minimiza o campo elétrico
    razao = R1/R2;

%segundo item (valores universais)

r1 = 19.5 * 10^(-3);
r2 = 37.75 * 10^(-3);
r3 = 37.97 * 10^(-3);
r4 = 42.5 * 10^(-3);
I = 150;
V = (69 * 10^3)/sqrt(3);

%nucleo
ro_central = 3.365 * 10^(-8);

%isolante
e_isolante_2 = 2.85 * e0;

%coroa
ro_coroa = 1.718 * 10^(-8);

%isolante externo
e_isolante_ext = 2.51 * e0;

    %formula vetor de poynting
    P = real(V* conj(I));
    
    %formula da densidade de energia elétrica
    
    %formula da desidade de energia magnetica
    
    %razao que minimiza o campo elétrico
    razao = R1/R2;
    
    %razao que minimiza o campo elétrico
    razao_linha = R1/R2;

%terceiro item

%distancia entre os cabos
d = 30 * 10^(-2);
%profundidade
h = 1;

    %valor campo eletromagnetico na superfície do solo
    eletrico_3 = 0;
    magnetico_3 = 0;
    
%quarto item

%distancia entre os cabos
d = 20 * 10^(-2);
%profundidade mínima
h_min = 50 * 10^(-2);

    %valor campo eletromagnetico na superfície do solo
    eletrico_4 = 0;
    magnetico_4 = 0;

%quinto item

%distancia entre os cabos
d = 2 * D;
%profundidade
h = 800 * 10^(-3);

l = 1.75 * d;

    %parametros da caixa de areia
    sigma_areia = 10^(-5);
    e_areia = 2;
    
    %parametros do solo
    sigma_solo = 10^(-3);
    e_solo = 10;
    
    %casamento de impedancias
    w_max =;
    



