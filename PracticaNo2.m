clc;
clear all;
close all;

%Sistemas Celulares
%Practica 2
%Elaborado por Ulices Yarit Mora Alavez
%Parte 1. Características de la señal OFDM
%Primero se generan los 8 bits aleatorios
bkt = randi([0,1], 1, 8);%16 bits

%Convertir los bits generados aleatoriamente a la representacion en
%sistemas celulares 1 pasa a -1 y 0 a 1.
bkt(bkt == 0) =-1;

% Agregar un valor adicional al final igual al último valor para
% representar correctamente el ultimo pulso
bktAux = [bkt, bkt(end)];

% Crear un vector de tiempo para representar los pulsos
sepPort = 1600;
sepPort2 = 2000;
sepPort3 = 3200;

%Tiempo útil de cada simbolo
Tu=1/sepPort;
Tu2=1/sepPort2;
Tu3=1/sepPort3;

% Tiempo de guarda
Ts=1/16000;

% Vector de tiempo de Ts de cada subportadora
T=0:Ts:Tu-Ts;
T2=0:Ts:Tu2-Ts;
T3=0:Ts:Tu3-Ts;

% Numero de Subportadoras
Nc=8;


t = 0:Tu:(length(bktAux)-1)*Tu;
t1 = Tu:Tu:(length(bktAux)-1)*Tu;
t2 = 0:Tu:(length(bktAux)-2)*Tu;
L = length(T)*160;
Fs=1/Ts;
f = Fs*(-(L/2):(L/2)-1)/L;
stairs(t, bktAux, 'LineWidth', 2);

% Ajustar los ejes para que cada pulso dure una unidad de tiempo
xlim([0, (length(bktAux)-1)*Tu]);
ylim([-1.5, 1.5]);

% Etiquetas de los ejes
xlabel('Tu');
ylabel('Magnitud');

% Título del gráfico
title('Representación de los bits Aleatorios');

% Añadir grid
grid on;

modsig = zeros(1,length(T));
modsig2 = zeros(1,length(T2));
modsig3 = zeros(1,length(T3));
modsig1 = zeros(1,length(T));
modsig12 = zeros(1,length(T2));
modsig13 = zeros(1,length(T3));
espectroS = zeros(1,length(f));
espectroS2 = zeros(1,length(f));
espectroS3 = zeros(1,length(f));

%Señal modulada
for c=1:Nc
    modsig = (modsig+bkt(c)*exp(1i*2*pi*sepPort));
    modsig2 = (modsig2+bkt(c)*exp(1i*2*pi*sepPort2));
    modsig3 = modsig3+bkt(c)*exp(1i*2*pi*T3(c-1)*sepPort3);
    modsig1(c,:) = bkt(c)*exp(1i*2*pi*T(c-1)*sepPort);
    modsig12(c,:) = bkt(c)*exp(1i*2*pi*T2(c-1)*sepPort2);
    modsig13(c,:) = bkt(c)*exp(1i*2*pi*T3(c-1)*sepPort3);
    figure
    subplot(2,1,1)
    plot(T,abs(modsig1(c,:)));
    grid on;
    title(['Magnitud de la señal modulada ' num2str(c) ' a 1600 Hz']);
    xlabel('\Deltat');
    ylabel('Magnitud');

    espectroS=fftshift(abs(fft(modsig1(c,:), L)));
    subplot(2,1,2)
    plot(f,abs(espectroS)/max(abs(espectroS)))
    grid on;
    title(['Magnitud del espectro de la señal modulada ' num2str(c) ' a 1600Hz']);
    xlabel('Hz');
    ylabel('Magnitud');


    figure
    subplot(2,1,1)
    plot(T2,abs(modsig12(c,:)));
    grid on;
    title(['Magnitud de la señal modulada ' num2str(c) ' a 2000 Hz']);
    xlabel('\Deltat');
    ylabel('Magnitud');

    espectroS2=fftshift(abs(fft(modsig12(c,:), L)));
    subplot(2,1,2)
    plot(f,abs(espectroS2)/max(abs(espectroS2)))
    grid on;
    title(['Magnitud del espectro de la señal modulada ' num2str(c) ' a 2000Hz']);
    xlabel('Hz');
    ylabel('Magnitud');



    figure
    subplot(2,1,1)
    plot(T3,abs(modsig13(c,:)));
    grid on;
    title(['Magnitud de la señal modulada ' num2str(c) ' a 3200 Hz']);
    xlabel('\Deltat');
    ylabel('Magnitud');

    espectroS3=fftshift(abs(fft(modsig13(c,:), L)));
    subplot(2,1,2)
    plot(f,abs(espectroS3)/max(abs(espectroS3)))
    grid on;
    title(['Magnitud del espectro de la señal modulada ' num2str(c) ' a 3200Hz']);
    xlabel('Hz');
    ylabel('Magnitud');
end


figure
plot(T, abs(modsig));
title('Magnitud de la señal multiplexada s(t) a 1600 Hz');
xlabel('\Deltat');
ylabel('Magnitud');
grid on

figure
plot(T2, abs(modsig2));
title('Magnitud de la señal multiplexada s(t) a 2000 Hz');
xlabel('\Deltat');
ylabel('Magnitud');
grid on

figure
plot(T3, abs(modsig3));
title('Magnitud de la señal multiplexada s(t) a 3200 Hz');
xlabel('\Deltat');
ylabel('Magnitud');
grid on


%calculando el espectro de s(t)
SF=fftshift(abs(fft(modsig, L)));
SF2=fftshift(abs(fft(modsig2, L)));
SF3=fftshift(abs(fft(modsig3, L)));
figure
subplot(3,1,1)
plot(f,abs(SF)/max(abs(SF)))
xlabel('Hz')
ylabel('Magnitud')
title('Magnitud de S(f) a 1600 Hz')
grid on

subplot(3,1,2)
plot(f,abs(SF2)/max(abs(SF2)))
xlabel('Hz')
ylabel('Magnitud')
title('Magnitud de S(f) a 2000 Hz')

subplot(3,1,3)
plot(f,abs(SF3)/max(abs(SF3)))
%stem(f,abs(espectroS1),'Color','r');
xlabel('Hz')
ylabel('Magnitud')
title('Magnitud de S(f) a 3200 Hz')
grid on


%%%%%%%%%%%%%%%%%%%%
%Segunda parte de la práctica
N=16;
N2=6;
N3=8;
N4=32;
nTs=Tu/N;
nTs2=Tu/N2;
nTs3=Tu/N3;
nTs4=Tu/N4;
n = 0:nTs:Tu-nTs;
n2 = 0:nTs2:Tu-nTs2;
n3 = 0:nTs3:Tu-nTs3;
n4 = 0:nTs4:Tu-nTs4;
sn = ifft(bkt,N);
sn2 = ifft(bkt,N2);
sn3 = ifft(bkt,N3);
sn4 = ifft(bkt,N4);
figure
stem(n,N*abs(sn),'Color','r');
plot(n,N*abs(sn), 'color', 'r', 'LineWidth',2)
xlabel('Tiempo')
ylabel('Magnitud')
title(['Comparación de la magnitud de s_n y s(t), N=' num2str(N)])
grid on
hold on

plot(T, abs(modsig), 'color', 'b', 'LineWidth',2);
legend('s_n','s(t)')


figure
subplot(2,2,1)
plot(n,N*abs(sn), 'color', 'r', 'LineWidth',2)
xlabel('Tiempo')
ylabel('Magnitud')
title(['Comparación de la magnitud de s_n y s(t), N=' num2str(N)])
grid on
hold on

plot(T, abs(modsig), 'color', 'b', 'LineWidth',2);
legend('s_n','s(t)')

subplot(2,2,2)
plot(n2,N2*abs(sn2), 'color', 'r', 'LineWidth',2)
xlabel('Tiempo')
ylabel('Magnitud')
title(['Comparación de la magnitud de s_n y s(t), N=' num2str(N2)])
grid on
hold on

plot(T, abs(modsig), 'color', 'b', 'LineWidth',2);
legend('s_n','s(t)')

subplot(2,2,3)
plot(n3,N3*abs(sn3), 'color', 'r', 'LineWidth',2)
xlabel('Tiempo')
ylabel('Magnitud')
title(['Comparación de la magnitud de s_n y s(t), N=' num2str(N3)])
grid on
hold on

plot(T, abs(modsig), 'color', 'b', 'LineWidth',2);
legend('s_n','s(t)')

subplot(2,2,4)
plot(n4,N4*abs(sn4), 'color', 'r', 'LineWidth',2)
xlabel('Tiempo')
ylabel('Magnitud')
title(['Comparación de la magnitud de s_n y s(t), N=' num2str(N4)])
grid on
hold on

plot(T, abs(modsig), 'color', 'b', 'LineWidth',2);
legend('s_n','s(t)')