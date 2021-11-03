
%supondo que o campo elétrico tenha apenas componente em x
E(x, t, w, alpha) = [cos(w*t)*exp(-alpha*x), 0, 0];
X = [x y z];
pretty(curl(E, X));

%supondo que o campo elétrico tenha apenas componente em y
E(x, t, w, alpha) = [0, cos(w*t)*exp(-alpha*x), 0];
dBdt = curl(E, X);
pretty(dBdt);

%campo magnetico tera apenas componente em Z
H = (-exp(-x*alpha)*alpha*cos(w*t))/ u;
result = int(H, t);
pretty(result);

%supondo que o campo elétrico tenha apenas componente em z
E(x, t, w, alpha) = [0, 0, cos(w*t)*exp(-alpha*x)];
pretty(curl(E, X));
