d = audioread('hello_world.wav');
plot(d);
size(d);
[d, Fs] = audioread('hello_world.wav');
%sound(d, Fs);

%inverter um sinal
d2 = flipud(d);

%comparando os plots dos sinais
subplot(2,1,1);
plot(d);
title('sinal original')
subplot(2,1,2);
plot(d2);
title('sinal invertido')

%salvando o arquivo invertido
audiowrite('helloworld_invertido.mp4', d2, Fs);
%lendo o arquivo invertido em novas variaveis
[d3,fs] = audioread('helloworld_invertido.mp4');
%tocando o arquivo invertido
%sound(d3, fs);

%acelerando o som em 2x
d4 = downsample(d, 2);
sound(d4, Fs);
sound(d4, Fs/2);



%d representa os dados armazenados no matlab
%Fs representa a frequencia da amostra
%audioread pode ler qualquer tipo de arquivo de som