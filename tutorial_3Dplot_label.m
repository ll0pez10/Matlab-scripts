X = linspace(0, 2*pi,100);
Y = sin(X);
Z = Y+ rand(1,100);
plot(X,Z,X,Y);

subplot(2,1,1);
plot(X,Y);
title('sin(x)')
subplot(2,1,2);
plot(X,Z);

Z = randn(100,100);
surf(Z);
contour(Z);

X = linspace(0, 2*pi,100);
Y = sin(X);
plot(X,Y)
title('A plot of sin(x)')
xlabel('my label for the x-axis')
ylabel('my label for the y-axis')




