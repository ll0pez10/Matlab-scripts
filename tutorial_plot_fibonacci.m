% plot(x, y)
% y2 = cos(x);
% plot(x,y,x,y2);
% plot();
% plot();

%os tamanhos dos arrays devem ser iguais
y =[1,1,2,3,5,8,13,21];
plot(y)
x = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7];
plot(x,y)

%inicio, fim e quantidade de pontos funcao linspace
x = linspace(0,2*pi,200);
y = sin(x);
plot(x,y)
y2 = cos(x);
plot(x,y,'--',x,y2,'.');
openfig('tutorial_plot.fig')