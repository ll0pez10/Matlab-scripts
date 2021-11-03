%a transformada de fourier  deixa ver a componente de frequencia do sinal
%matlab tem uma funcao para calcular a FFT

x = linspace(0,5*pi,5000);
y = sin(x);
Y = fft(y);
%real da a parte real de um numero imaginairo
plot(real(Y))

%a transformada de fourier é periodica com periodo = 2pi
%o sinal no DT tem o mesmo tamanho do sinal no DF
d = audioread('hello_world.wav');
D = fft(d(:,1));
plot(real(D));

%short time fourier transform == spectogram
y2 = sin(x.*x);
spectrogram(y2,500);
%principio de incerteza na transformada de fourier, se queremos uma
%frequencia muito acurada precisamos ter uma janela de tempo muito longa
%mas isso nos da frequencias diferentes em tempos diferentes se queremos
%saber a frequencia em um tempo especifico precisamos estreitar a janela
%mas assim nao temos uma frequencia muito acurada