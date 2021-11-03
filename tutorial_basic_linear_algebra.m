A = [ 3 2 -1; 2 -2 4; -1 0.5 -1]
b = [ 1 -2 -2]

%Ax = b
%resolver sistemas lineares
x = A\b

%determinante de uma matriz
det(A)

%autovalor de uma matriz
eig(A)

%autovetores de uma matriz
%D matriz de autovalores
%V matriz de autovetores
[V, D] = eig(A)

%comparacao de valores
1 == 0

%normalizar os vetores (euclidiano)
norm(A*V - V*D)

