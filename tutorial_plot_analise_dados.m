x = 1:10;
%grafico de barras
bar(x)

x = randn(1000,1);
%histograma, 2 argumento e o numero de intervalos
hist(x, 50)

%grafico de pizza
x = 1:10;
pie(x)

%plot scatter
x = linspace(0, 2*pi,1000);
y = 10*sin(x) + randn(1,1000);
plot(x,y)
scatter(x,y)

x = randn(1000,1)*2;
y = 5*sin(x) + randn(1000,1);
scatter(x,y)

