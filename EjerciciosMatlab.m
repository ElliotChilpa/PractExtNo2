clc;
clear all;
close all;

%% Practica Numero 2
% Sistemas Celulares
%% Creación de 8 Bits Aleatorios y Grafica

% Definir 8 bits (Secuencia aleatoria de 0s y 1s)
bits = randi([0, 1], 1, 8); % Genera 8 bits aleatorios

% Cambiar de 0s y 1s a señalización (-1s y 1s)
bits(bits == 0) = -1;

% Parámetros de la señal de pulsos
amplitud = 1;          % Amplitud de los pulsos
duracion_pulso = 1;   % Duración de cada pulso
tiempo_muestreo = 0.01; % Intervalo de muestreo

% Crear un vector de tiempo para la señal de pulsos eje x
t8bit = 0:tiempo_muestreo:duracion_pulso*length(bits);

pulsos = zeros(1, length(t8bit));

% Generar la señal de pulsos en base a la secuencia de bits
for i = 1:length(bits)
    inicio = (i - 1) * duracion_pulso / tiempo_muestreo + 1;
    fin = i * duracion_pulso / tiempo_muestreo;
    pulsos(inicio:fin) = bits(i) * amplitud;
end

% Graficar la señal de pulsos
figure()
plot(t8bit, pulsos, "r","LineWidth",2);
xlabel('Tiempo');
ylabel('Amplitud');
title('Bits Aleatorios enre (1 y -1)');
grid on;
axis([0, t8bit(end), -1.5, 1.5]);

%% Realizar señales
% Agreagamos a nuestra señal el valor de nuetro ultimo bit Cyclic prefix CP
% o tiempo de guarga Tg
bitsAux = [bits, bits(end)];

% Separación entre portadoras DeltaF en Frecuencias
sepPort = 1600;
sepPort1 = sepPort;
sepPort2 = sepPort * 2;
sepPort3 = sepPort * 3;
sepPort4 = sepPort * 4;
sepPort5 = sepPort * 5;
sepPort6 = sepPort * 6;
sepPort7 = sepPort * 8;
sepPort8 = sepPort * 9;

%Tiempo útil de cada simbolo
Tu1 = 1/sepPort;
Tu2 = 1/sepPort2;
Tu3 = 1/sepPort3;
Tu4 = 1/sepPort4;
Tu5 = 1/sepPort5;
Tu6 = 1/sepPort6;
Tu7 = 1/sepPort7;
Tu8 = 1/sepPort8;

% Tiempo de guarda
Ts = 1/16000;

% Vector de tiempo de Ts de cada subportadora
T1 = 0:Ts:Tu1; % La primera subportadora inicia en 0, con tiempo de muestreo ts y el final es el tiempo util
T2 = 0:Ts:Tu2;
T3 = 0:Ts:Tu3;
T4 = 0:Ts:Tu4;
T5 = 0:Ts:Tu5;
T6 = 0:Ts:Tu6;
T7 = 0:Ts:Tu7;
T8 = 0:Ts:Tu8;

% Numero de Subportadoras
Nc = 8;

% 
%t = 0:Tu:(length(bitsAux)-1)*Tu;
%t1 = Tu:Tu:(length(bitsAux)-1)*Tu;
%t2 = 0:Tu:(length(bitsAux)-2)*Tu;

L = length(T1)*160;

Fs = 1 / Ts;

f = Fs*(-(L/2):(L/2)-1)/L;


modsig = zeros(1,length(T1));
modsig2 = zeros(1,length(T2));
modsig3 = zeros(1,length(T3));

%modsig1 = zeros(1,length(T));
%modsig12 = zeros(1,length(T2));
%modsig13 = zeros(1,length(T3));

%espectroS = zeros(1,length(f));
%espectroS2 = zeros(1,length(f));
%espectroS3 = zeros(1,length(f));

%Señal modulada
%exp(1j * 2 * pi * fc * t);




%%%

% Define number of subcarriers
N = 8; % 8 subcarriers

% Define sampling frequency (Fs)
Fs = 1e6; % Example: 1 MHz sampling frequency

% Define subcarrier spacing (Delta_f)
Delta_f = Fs/N; % Subcarrier spacing equal to Fs/N

% Generate time vector
t = 0 : 1/Fs : 1/Fs; % Time vector with sampling period 1/Fs

% Generate subcarriers
subcarriers = zeros(N, length(t));

for k = 0:N-1
    % Generate kth subcarrier
    subcarriers(k+1,:) = exp(1j*2*pi*k*Delta_f*t);
end

% Define frequency selective channel (example)
H = zeros(N, length(t)); % Channel matrix

for k = 0:N-1
    % Introduce frequency-selective fading
    H(k+1,:) = exp(1j*2*pi*rand*k*Delta_f*t);
end

% Apply channel to subcarriers
A_reshaped = reshape(H, [, 2]);

received_signal = H * subcarriers;

% Plot subcarriers in time domain (optional)
figure;
for k = 1:N
    plot(t, real(received_signal(k,:)));
    hold on;
end
xlabel('Time (s)');
ylabel('Amplitude');
title('Received Signal in Time Domain');

% Calculate FFT of received signal
fft_received_signal = fft(received_signal, [], 2);

% Plot magnitude spectrum of subcarriers
figure;
for k = 1:N
    f = Fs/length(fft_received_signal(k,:))*[-length(fft_received_signal(k,:))/2+1:length(fft_received_signal(k,:))/2];
    plot(f, abs(fft_received_signal(k,:)));
    hold on;
end
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title('Magnitude Spectrum of Subcarriers');
