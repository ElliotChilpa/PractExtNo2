clear all
close all
%% PARTE 1. CARACTERÍSTICAS DE LA SEÑAL OFDM.

% Genere 8 bits aleatorios, 𝑏𝑘(𝑡), y represéntelos con pulsos rectangulares de duración 
% 𝑇𝑢 y cuyas magnitudes pueden ser +1 y -1.

Tu = 1/1600; % Duración del símbolo 
Nc = 8; % Número de subportadoras
deltaF = [1600 2000 3200]; % Separación entre subportadoras 1600 Parte 1. y (2000 y 3200)
deltaT = 1/160000; % Resolución para el vector tiempo t, Recomentación del profe 1/16000

% Generar 8 bits aleatorios de magnitud +1 y -1
bk = (2 * randi([0,1], 1, 8) - 1).'; 

% Cada bit representa un simbolo, se replica el Número de muestras por símbolo (Tu/Ts) que en este caso es 10
bk_t = kron(bk, ones(1, Tu/deltaT)); % Matriz de 8x10

% Transforma la matriz de 8x10 en un vector de 1x80 concatenando filas solo
% para representar tren de bits
bk_t_Tren = reshape(bk_t.', 1, []); 

% Obtenemos el vector tiempo por simbolo
t = 0:deltaT:Tu-deltaT;

% Vector tiempo para el tren de pulsos y graficarlo
t_total = 0:deltaT:(length(bk_t_Tren)-1)*deltaT;

% Graficar el tren de pulsos 
figure(1)
plot(t_total, bk_t_Tren, 'LineWidth', 2, 'Color', 'b'); 
title('Tren de Pulsos de 8 Bits Generados Aleatoriamente');
xlabel('Tiempo (s)');
ylabel('Amplitud');
%ylim([-1.5, 1.5]); 
grid on;

% Graficar cada símbolo individualmente en subplots
figure(2)
for i = 1:Nc
    subplot(2, 4, i)
    %subplot(8, 1, i)
    plot(t,bk_t(i,:), 'LineWidth', 2, 'Color', 'r');
    title(['Simbolo ', num2str(i)]); % Título de cada subplot
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    ylim([-1.5, 1.5]); 
    grid on;
end
sgtitle('Simbolos con 10 Muestras*Simbolo'); % Título principal de la figura


% Cada símbolo se multiplíca por una portadora de la forma 𝑠𝑘(𝑡) = 𝑒^{𝑗∙2𝜋∙𝑘∙Δf∙t} , 𝑘 ∈ [0,1, . .7].
%inciso 1
sk_t{1} = zeros(8,1);

%inciso 3
sk_t{2} = zeros(8,1);
sk_t{3} = zeros(8,1);

%Bucle Externo ( k = 0:7 )
%Itera sobre cada uno de los 8 símbolos

% Bucle Interno (tiempo=1:length(t)):
% Itera sobre cada una de las muestras de tiempo OSEA la duración del
% simbolo que son 10 muestras
for k=0:7
    for tiempo=1: length(t)
        % Para inciso 1
        % Cada simbolo tiene 10 muestras por lo que sk_t debemos de tener
        % 8 (sk_t) ya que son 8 simbolos y cada simbolo tiene 10 muestras
        sk_t{1}(k+1,tiempo) = exp(j*2*pi*deltaF(1)*t(tiempo)*k);

        %Para inciso 3
        sk_t{2}(k+1,tiempo) = exp(j*2*pi*deltaF(2)*t(tiempo)*k);
        sk_t{3}(k+1,tiempo) = exp(j*2*pi*deltaF(3)*t(tiempo)*k);
    end
end

%Para inciso 1
sk_bk_t{1} = bk_t.* sk_t{1}; % Multiplicamos las señales Bk por la portadoras Sk

%Para inciso 3
sk_bk_t{2} = bk_t.* sk_t{2}; % Multiplicamos las señales Bk por la portadoras Sk
sk_bk_t{3} = bk_t.* sk_t{3}; % Multiplicamos las señales Bk por la portadoras Sk


% Graficar las señales moduladas (inciso 1)
figure(3)
for i = 1:Nc
    subplot(2, 4, i)
    plot(t, real(sk_bk_t{1}(i, :)));
    title(['Símbolo ', num2str(i)]);
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    ylim([-1.5, 1.5]); 
    grid on;
end
sgtitle('Símbolos Modulados con ΔF = 1600 '); 

% Cálculo de la magnitud de las señales moduladas (inciso 1)
mag_sk_bk_t{1} = abs(sk_bk_t{1});

% Magnitud de las señales moduladas inciso 3 ( deltaF(2) y deltaF(3) )
mag_sk_bk_t{2} = abs(sk_bk_t{2});
mag_sk_bk_t{3} = abs(sk_bk_t{3});

figure(4) %Graficamos la señal multiplicada
 for i=1:Nc % las 8 subportadoras
    subplot(2,4,i)
    plot(t, mag_sk_bk_t{1}(i, :),"LineWidth", 2,"Color", "b");
    title(['| Símbolo ', num2str(i), ' |']);
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    ylim([-1.5, 1.5]); % Establecer límites del eje y para visualizar mejor
    grid on;
end
sgtitle('Magnitud de los Símbolos Modulados con F(1600)');

N = 1000000; % Número de puntos para la FFT

% Creamos un vetor para graficar todas las sub-portadoras
% Espaciamos de manera lineal dos valores que nos dan un intervalo
w = linspace(-1/deltaT/2,1/deltaT/2,N)*2*pi;

% Para los espectros de cada subportadora
%Para inciso 1
sk_bk_f{1} = zeros(8,N);

%Para inciso 3
sk_bk_f{2} = zeros(8,N);
sk_bk_f{3} = zeros(8,N);

figure(5) % Graficamos los espectros de cada subportadora
% Para graficar el espectro calculamos la fft de los 8 simbolos
for i=1:Nc
    %Para inciso 1
    subplot(2,4,i)
    %subplot(8,1,i)
    sk_bk_f{1}(i,:) = fftshift(fft(sk_bk_t{1}(i,:),N))*deltaT; % fft
    plot(w/(2*pi),abs(sk_bk_f{1}(i,:))) % raices cada w/2pi
    axis([-4000+1600*(i-1) 4000+1600*(i-1) 0 1.1*max(abs(sk_bk_f{1}(i,:)))])
    hold on;

    % Esto es para marcar el punto medio del espectro con una línea punteada
    mid_freq = 1600*(i-1); % Frecuencia central de las portadoras
    xline(mid_freq, '--r', 'LineWidth', 1.5);

    % Nota de texto para indicar la frecuencia central
    text(mid_freq, 0.9*max(abs(sk_bk_f{1}(i,:))), ['f = ', num2str(mid_freq), ' Hz'], 'Color', 'red', 'HorizontalAlignment', 'right');
    hold off;
    title(['Espectro Sk(f) ', num2str(i)]);
    grid on;

    %Para inciso 3
    sk_bk_f{2}(i,:) = fftshift(fft(sk_bk_t{2}(i,:),N))*deltaT;
    sk_bk_f{3}(i,:) = fftshift(fft(sk_bk_t{3}(i,:),N))*deltaT;
end


figure(6) %Graficamos todos los espectros traslapados
for i = 1:Nc % Nc numero de sub-portadoras
    plot(w/(2*pi),(abs(sk_bk_f{1}(i,:))), "LineWidth", 0.75)
    hold on
    grid on;
end
hold off;
title('Espectros de las señales OFDM')
xlabel('Frecuencia (Hz)');
ylabel('|bk * Sk|');
axis([-5000 15000 0 1.1*max(abs(sk_bk_f{1}(1,:)))]);

figure(7)
% Graficar el espectro en la columna 1
for i = 1:Nc
    subplot(Nc, 2, 2*i-1) % Columna 1 para el espectro
    plot(w/(2*pi), abs(sk_bk_f{1}(i,:))) % Raices cada w/2pi
    axis([-4000+1600*(i-1) 4000+1600*(i-1) 0 1.1*max(abs(sk_bk_f{1}(i,:)))])
    hold on;
    % Marcar el punto medio del espectro con una línea punteada
    mid_freq = 1600 * (i-1); % Frecuencia central de las portadoras
    xline(mid_freq, '--r', 'LineWidth', 1.5);
    % Nota de texto para indicar la frecuencia central
    text(mid_freq, 0.9 * max(abs(sk_bk_f{1}(i,:))), ['f = ', num2str(mid_freq), ' Hz'], 'Color', 'red', 'HorizontalAlignment', 'right');
    hold off;
    title(['Espectro Sk(f) ', num2str(i)]);
    grid on;
end

% Graficar la magnitud en la columna 2
for i = 1:Nc
    subplot(Nc, 2, 2*i) % Columna 2 para la magnitud
    plot(t, mag_sk_bk_t{1}(i, :), "LineWidth", 2, "Color", "b");
    title(['| Símbolo ', num2str(i), ' |']);
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
    ylim([-1.5, 1.5]); % Establecer límites del eje y para visualizar mejor
    grid on;
end

% Título general para la figura
sgtitle('Espectro y Magnitud de los Símbolos Modulados con F(1600)');

%% ---------------------------------- INCISO 2 ----------------------------
%Obtenga la señal multiplexada 𝑠(𝑡) = ∑ 𝑏𝑘(𝑡) ∙ 𝑠𝑘(𝑡) 𝑘=0 y grafique su magnitud. 
%Además, obtenga el espectro de esta señal, al que se le denominará S(f), y grafique la magnitud de éste

% Obtenemos s(t) de cada simbolo que tiene duración tu y tiene 10 muestras
% Para inciso 2
s_t{1} = zeros(1,length(t));

% Para inciso 3
s_t{2} = zeros(1,length(t));
s_t{3} = zeros(1,length(t));

% Para cada simbolo que tiene una duracion de Tu llamada t y tenemos 8 simb
for tiempo = 1: length(t)
    for k = 0:7 % Para las 8 sub-portadoras
        %Para inciso 2
        s_t{1}(1,tiempo) = s_t{1}(1,tiempo) + sk_bk_t{1}(k+1,tiempo);

        %Para inciso 3
        s_t{2}(1,tiempo) = s_t{2}(1,tiempo) + sk_bk_t{2}(k+1,tiempo);
        s_t{3}(1,tiempo) = s_t{3}(1,tiempo) + sk_bk_t{3}(k+1,tiempo);
    end
end

figure(8) 
subplot(2,1,1) % Graficamos la señal multiplexada en función del tiempo
plot(t, abs(s_t{1}),"LineWidth", 1.25, "Color", "r")
ylim([0 1.1*max(abs(s_t{1}))])
title('Señal OFDM Multiplexada Tiempo S(t)')
xlabel('Tiempo (s)');
ylabel('Amplitud');
grid on

% Obtenemos la señal multiplexado en frecuencia
% Para inciso 2
S_f{1} = fftshift(fft(s_t{1},N))*deltaT;

% Para inciso 3
S_f{2} = fftshift(fft(s_t{2},N))*deltaT;
S_f{3} = fftshift(fft(s_t{3},N))*deltaT;


subplot(2,1,2) % Graficamos la señal multiplexada en función de la frecuencia
plot(w/(2*pi),abs(S_f{1}),"LineWidth", 1.25, "Color", "b")
axis([-15000 25000 0 1.1*max(abs(S_f{1}))]);
title('Señal OFDM Multiplexada Frec S(f)')
xlabel('Frecuencia (Hz)');
ylabel('|S(f)|');
grid on


%% Repita las actividades 1 y 2 para Δf=2,000 y 3,200 Hz.

figure(9) % Graficamos las señales multiplexadas para f=1600, f=2000, f=3200 

subplot(3,1,1) % Graficamos la señal multiplexada f = 1600
plot(w/(2*pi),abs(S_f{1}), "LineWidth", 1.25,"Color", 'r') 
axis([-15000 25000 0 1.1*max(abs(S_f{1}))]);
title('Señal OFDM Multiplexada Frec F( 1600 )')
xlabel('Frecuencia (Hz)');
ylabel('|S(f)|');
grid on

subplot(3,1,2)%Graficamos la señal multiplexada f = 2000
plot(w/(2*pi),abs(S_f{2}), "LineWidth", 1.25,"Color", 'b')
axis([-15000 25000 0 1.1*max(abs(S_f{2}))]);
title('Señal OFDM Multiplexada Frec F( 2000 )')
xlabel('Frecuencia (Hz)');
ylabel('|S(f)|');
grid on

subplot(3,1,3)%Graficamos la señal multiplexada f = 3200
plot(w/(2*pi),abs(S_f{3}), "LineWidth", 1.25,"Color", 'magenta')
axis([-15000 30000 0 1.1*max(abs(S_f{3}))]);
title('Señal OFDM Multiplexada Frec F( 3200 )')
xlabel('Frecuencia (Hz)');
ylabel('|S(f)|');
grid on

%% PARTE 2. IMPLEMENTACIÓN DIGITAL DE LA SEÑAL OFDM.
% Genere la señal Sn. Para esta actividad use inicialmente N = 16.
% Usamos la magnitud de S(t)

modulo = [6 8 16 32];

%Obtenemos todas las frecuencias de N

for i=1:length(modulo)
    fsn{i} = modulo(i)*1/Tu; % Aquí calculamos las muestras por simbolo = frecuencia de muestreo
    Tsn{i} = 1/fsn{i}; % Tiempo de muestreo porque nuestra magnitud de S(t) esta en tiempo
    tn{i} = 0:Tsn{i}:Tu-Tsn{i}; % Creamos un vector tiempo de todo lo que dura el simbolo

    s_tn{i} = modulo(i)*ifft(bk, modulo(i)); % Calculamos la transformada inversa
end

% Graficamos la magnitud de nuestra señal en el tiempo Y Muestreamos la
% señla con N = 16
figure(10)
plot(t, abs(s_t{1}),"LineWidth", 1.25,"Color","b");
hold on;
stem(tn{3}, abs(s_tn{3}), 'LineWidth', 1.25, 'Color', 'magenta');
title('Muestreo de S(t) con N = 16');
%xlabel('Tiempo (s)');
%ylabel('Amplitud');

grid on;
    
figure(11)
for i = 1:length(modulo)
    subplot(length(modulo),1,i)
    plot(t, abs(s_t{1}),"LineWidth", 1.25,"Color","b")
    hold on;
    stem(tn{i},abs(s_tn{i}),"LineWidth", 1.25,"Color",'magenta')
    %xlabel('Tiempo (s)');
    %ylabel('Amplitud');
    title(['Muestreo de S(t) con N = ', num2str(modulo(i))]);
    grid on;
    hold off
end

