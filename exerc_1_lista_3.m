%nome: Luan Lopes dos Santos
%DRE: 115173479
%disciplina: circuitos elétricos em CA
%professor: Heloi
%exercicio 1 lista 3 

%operador a formato retangular
a =  -0.5+0.866i;
a2 = a^2;

%matriz A
A = [1 1 1; a2 a 1; a a2 1];
disp(A);
%valores
Za = 2;
Zb = 2;
Zc = 2;
Va = 10;
Vb = 20-17.32i;
Vc = 20+17.32i;
Iabc(1) = 2;

%matriz de impedâncias com a resistencia do neutro, todas acopladas
Z = [Za+3*Rn 0 0; 0 Zb+3*Rn 0; 0 0 Zc+3*Rn];
disp(inv(Z));

%matriz de tensões original
Vabc = [Va; Vb; Vc];

% VALORES matriz de tensoes 120 
V = inv(A)* Vabc;
disp(V);
disp(inv(A));

%VALORES matriz de correntes 120
I = A * inv(Z) * inv(A) * V;

%VALORES matriz de correntes original 
Iabc = A * I;
In = 3*I(3);

%verificacao dos resultados o somatorio das correntes deve ser igual a
%corrente no neutro
prova = Iabc(1) + Iabc(2) + Iabc(3) - In;

%potencia aparente da configuração
S120 = conj(V')*conj(I);
    
Sabc = conj(Vabc')*conj(Iabc);