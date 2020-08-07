%nome: Luan Lopes dos Santos
%DRE: 115173479
%disciplina: circuitos elétricos em CA
%professor: Heloi

%operador a formato retangular
a =  -0.5+0.866i;
a2 = a^2;

%matriz A
A = [1 1 1; a2 a 1; a a2 1];

%valores
Za = 15;
Zb = 10;
Zc = 20;
Rn = 5;
Va = 8.66+5i;
Vb = 10-17.32i;
Vc = -17.32+10i;

%matriz de impedâncias
Z = [Za 0 0; 0 Zb 0; 0 0 Zc];

%matriz de tensões original
Vabc = [Va; Vb; Vc];

% VALORES matriz de tensoes 120 
V = inv(A)* Vabc;

%VALORES matriz de correntes 120
I = A * inv(Z) * inv(A) * V;

%VALORES matriz de correntes original 
Iabc = A * I;

% %neutro solidamente aterrado, componentes 1 e 2 permanecem iguais
In = 3*I(3);

%verificacao dos resultados o somatorio das correntes deve ser igual a
%corrente no neutro
prova = Iabc(1) + Iabc(2) + Iabc(3) - In;

%potencia aparente da configuração
S120 = conj(V')*conj(I);
    disp (3*S120);
Sabc = conj(Vabc')*conj(Iabc);
    disp(Sabc);
