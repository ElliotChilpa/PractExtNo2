clear all
close all
%%  ********************************* PARTE 1 ******************************
%% ---------------------------------- INCISO 1 ----------------------------
% Genere 8 bits aleatorios, 𝑏𝑘(𝑡), y represéntelos con pulsos rectangulares de duración 
% 𝑇𝑢 y cuyas magnitudes pueden ser +1 y -1.
Tu = 1/1600;
deltaT = 1/160000;
Nc = 8;
deltaF = [1600 2000 3200]; %Los otros dos valores serviran para el inciso 3

bk = (2*randi([0,1],1,8)-1).'; %8 bits generados aleatorios de magnitud +1 y -1
bk_t = kron(bk, ones(1, Tu/deltaT));%Repetira los bits dependiendo deltaT en este caso, por cada simbolo seran 10 muestras

%Obtenemos el vector tiempo
t = 0:deltaT:Tu-deltaT;

figure(1)
 for i=1:Nc
    subplot(2,4,i)
    plot(t,bk_t(i,:));
    title(['Bit ', num2str(i)]);
    ylim([-1.5, 1.5]); % Establecer límites del eje y para visualizar mejor
    grid on;
end

% A cada uno de estos símbolos multiplíquelo por una portadora de la forma 𝑠𝑘(𝑡) = 𝑒^{𝑗∙2𝜋∙𝑘∙Δf∙t} , 𝑘 ∈ [0,1, . .7].
%inciso 1
sk_t{1} = zeros(8,1);

%inciso 3
sk_t{2} = zeros(8,1);
sk_t{3} = zeros(8,1);

for k=0:7
    for tiempo=1: length(t)
        %Para inciso 1
        sk_t{1}(k+1,tiempo) = exp(j*2*pi*deltaF(1)*t(tiempo)*k);

        %Para inciso 3
        sk_t{2}(k+1,tiempo) = exp(j*2*pi*deltaF(2)*t(tiempo)*k);
        sk_t{3}(k+1,tiempo) = exp(j*2*pi*deltaF(3)*t(tiempo)*k);
    end
end

%Para inciso 1
sk_bk_t{1} = bk_t.* sk_t{1}; % ---> Multiplicamos la señal

%Para inciso 3
sk_bk_t{2} = bk_t.* sk_t{2}; % ---> Multiplicamos la señal
sk_bk_t{3} = bk_t.* sk_t{3}; % ---> Multiplicamos la señal

figure(2) %Graficamos la señal multiplicada
 for i=1:Nc
    subplot(2,4,i)
    plot(t,abs(sk_bk_t{1}(i,:)));
    title(['Bit ', num2str(i)]);
    ylim([-1.5, 1.5]); % Establecer límites del eje y para visualizar mejor
    grid on;
end

N=1000000;
%Espaciamos de manera lineal dos valores que nos dan un intervalo
w=linspace(-1/deltaT/2,1/deltaT/2,N)*2*pi;

figure(3) %Graficamos los espectros de cada bit
%Para inciso 1
sk_bk_f{1} = zeros(8,N);

%Para inciso 3
sk_bk_f{2} = zeros(8,N);
sk_bk_f{3} = zeros(8,N);

for i=1:Nc
    %Para inciso 1
    subplot(2,4,i)
    sk_bk_f{1}(i,:) = fftshift(fft(sk_bk_t{1}(i,:),N))*deltaT;
    plot(w/(2*pi),abs(sk_bk_f{1}(i,:)))
    axis([-4000+1600*(i-1) 4000+1600*(i-1) 0 1.1*max(abs(sk_bk_f{1}(i,:)))])
    title('Espectro Sk(f)', i)
    grid on;

    %Para inciso 3
    sk_bk_f{2}(i,:) = fftshift(fft(sk_bk_t{2}(i,:),N))*deltaT;
    sk_bk_f{3}(i,:) = fftshift(fft(sk_bk_t{3}(i,:),N))*deltaT;
end

figure(4) %Graficamos todos los espectros traslapados
for i=1:Nc
    plot(w/(2*pi),(abs(sk_bk_f{1}(i,:))))
    hold on
    grid on;
end
hold off;
title('Espectros de las señales bk(t)*sk(t)')
xlabel('Hz');
ylabel('|Bk(f)*Sk(f)|');
axis([-5000 15000 0 1.1*max(abs(sk_bk_f{1}(1,:)))]);

%% ---------------------------------- INCISO 2 ----------------------------
%Obtenga la señal multiplexada 𝑠(𝑡) = ∑ 𝑏𝑘(𝑡) ∙ 𝑠𝑘(𝑡) 𝑘=0 y grafique su magnitud. 
%Además, obtenga el espectro de esta señal, al que se le denominará S(f), y grafique la magnitud de éste

%-------------  Obtenemos s(t) ------------------------
%Para inciso 2
s_t{1} = zeros(1,length(t));

%Para inciso 3
s_t{2} = zeros(1,length(t));
s_t{3} = zeros(1,length(t));

for tiempo=1: length(t)
    for k=0:7
        %Para inciso 2
        s_t{1}(1,tiempo) = s_t{1}(1,tiempo) + sk_bk_t{1}(k+1,tiempo);

        %Para inciso 3
        s_t{2}(1,tiempo) = s_t{2}(1,tiempo) + sk_bk_t{2}(k+1,tiempo);
        s_t{3}(1,tiempo) = s_t{3}(1,tiempo) + sk_bk_t{3}(k+1,tiempo);
    end
end

figure(5) 
subplot(2,1,1) %Graficamos la señal multiplexada en función del tiempo
plot(t, abs(s_t{1}))
ylim([0 1.1*max(abs(s_t{1}))])
title('Señal Multiplexada en tiempo s(t)')
xlabel('s');
ylabel('|s(t)|');
grid on

%-------------  Obtenemos S(f) -----------------------
%Para inciso 2
S_f{1} = fftshift(fft(s_t{1},N))*deltaT;

%Para inciso 3
S_f{2} = fftshift(fft(s_t{2},N))*deltaT;
S_f{3} = fftshift(fft(s_t{3},N))*deltaT;


subplot(2,1,2)%Graficamos la señal multiplexada en función de la frecuencia
plot(w/(2*pi),abs(S_f{1}))
axis([-15000 25000 0 1.1*max(abs(S_f{1}))]);
title('Señal Multiplexada en frecuencia S(f)')
xlabel('Hz');
ylabel('|S(f)|');
grid on

%% ---------------------------------- INCISO 3.1 con f=2000 ----------------------------
%Repita las actividades 1 y 2 para Δf=2,000 y 3,200 Hz.


figure(6) %Graficamos las señales multiplexadas para f=1600, f=2000, f=3200 

subplot(3,1,1)%Graficamos la señal multiplexada f=1600
plot(w/(2*pi),abs(S_f{1}))
axis([-15000 25000 0 1.1*max(abs(S_f{1}))]);
title('Señal Multiplexada en frecuencia S(f) f=1600')
xlabel('Hz');
ylabel('|S(f)|');
grid on

subplot(3,1,2)%Graficamos la señal multiplexada f=2000
plot(w/(2*pi),abs(S_f{2}))
axis([-15000 25000 0 1.1*max(abs(S_f{2}))]);
title('Señal Multiplexada en frecuencia S(f) f=2000')
xlabel('Hz');
ylabel('|S(f)|');
grid on

subplot(3,1,3)%Graficamos la señal multiplexada f=3200
plot(w/(2*pi),abs(S_f{3}))
axis([-15000 30000 0 1.1*max(abs(S_f{3}))]);
title('Señal Multiplexada en frecuencia S(f) f=3200')
xlabel('Hz');
ylabel('|S(f)|');
grid on

%%  ********************************* PARTE 2 ******************************
%Utilice la IFFT para crear las muestras de s(t) a partir de los bits utilizados en la 
%Actividad 1, es decir, genere la señal sn. Para esta actividad use inicialmente N=16.
%Determine la magnitud de sn.

modulo = [6 8 16 32];

%Obtenemos todas las frecuencias
figure(7)
for i=1:length(modulo)
    fsn{i} = modulo(i)*1/Tu;
    Tsn{i} = 1/fsn{i};
    tn{i} = 0:Tsn{i}:Tu-Tsn{i};

    s_tn{i} = modulo(i)*ifft(bk, modulo(i));

    subplot(length(modulo),1,i)
    plot(t, abs(s_t{1}))
    hold on
    stem(tn{i},abs(s_tn{i}))
    title(['N =', num2str(modulo(i))])
    hold off
end
% fsn{1} = modulo() 
% fs1=modulo(4)*1/Tu;
% Ts1=1/fs1; 
% t1=0:Ts1:Tu-Ts1;
% figure(7)
% s_tIFFT = modulo(4)*ifft(bk,modulo(4));
% stem(t1,abs(s_tIFFT))
% hold on
% plot(t, abs(s_t{1}))
% hold off
