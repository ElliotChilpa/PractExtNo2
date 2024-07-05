clc
close all

% Generamos Bits Aleatorios 
% 
bits = 2 * randi([0, 1], 1, 8) - 1;
tu = 1/1600;
delta_f=1600;
fs=1/16000;
t = 0:fs:tu-fs;
Nc=8;
T=0:fs:tu-fs;
L = length(T)*160;
f_total = (1/fs)*(-(L/2):(L/2)-1)/L;

bk = [];
aux = [];

for i = 1:Nc
    if bits(i) == -1
        aux = -ones(size(t)); 
    else
        aux = ones(size(t)); 
    end
    bk = [bk aux];
end
clear aux

% Vector de tiempo para todos los pulsos
t_total = linspace(0,tu-fs,length(bk));

figure(1)
stairs(t_total, bk, 'b', 'LineWidth', 2);
xlabel('Tiempo [s]');
ylabel('Amplitud');
title('Señal b_k(t)');
ylim([-1.5, 1.5]);
%figure_style(hfig,'fig1')

s_t = zeros(Nc,length(bk));

%DOMINIO DEL TIEMPO
figure(2)
tiledlayout(4,2)
for j=1:Nc
    nexttile
    sk = exp(1i*2*pi*(j-1)*delta_f*t_total);
    mod=  conv(bits(j),sk);
    s_t(j,:) = mod; 
    plot(t_total,abs(mod),'MarkerFaceColor','#0072BD','MarkerEdgeColor','#0072BD')
    title(['Símbolo = ' num2str(j)])
    ylim([-1.5, 1.5]);
    xlabel('Tiempo [s]');
    ylabel('Amplitud');
end
sgtitle('Modulación de los símbolos en el dominio del tiempo f(t)')

%DOMINIO DE LA FRECUENCIA

figure(3)
tiledlayout(4,2)
for j=1:Nc
    nexttile
    sk = exp(1i*2*pi*j*delta_f*t_total);
    mod= abs(fftshift(fft( bits(j).*sk,16000)));
    % Vector de frecuencia para todos los pulsos
    f_total = linspace(-delta_f*8,delta_f*16,length(mod));

    plot(f_total,mod,'MarkerFaceColor','#D95319','MarkerEdgeColor','#D95319')
    title(['Símbolo = ' num2str(j)])
    xlabel('Frecuencia [Hz]');
    xlim([-10 (delta_f*9)+10 ])
end
sgtitle('Modulación de los símbolos en el dominio de la frecuencia F(\omega)')


figure (4)
tiledlayout(2,1)
nexttile
suma_st= sum(s_t);
plot(t_total,abs(suma_st));
title('Señal multiplexada s(t) en tiempo')
xlabel('Tiempo [s]');
nexttile
plot(f_total,abs(fftshift(fft(suma_st,16000))));
title('Señal multiplexada S(f) en frecuencia')
xlabel('Frecuencia [Hz]');

%-------- df = 1600, 2000, 3,2000 -----


%DOMINIO DEL TIEMPO
frecuencias = cell(1, 3);
freq=[1600,2000,3200];
for k=1:3
    figure(5+k)
    tiledlayout(8,1)
    for j=1:Nc
        nexttile
        sk = exp(1i*2*pi*(j-1)*freq(k)*t_total);
        mod= bits(j)*sk;
        frecuencias{k}(j,:) = mod; 
        plot(t_total,abs(mod))
        title(['Símbolo = ' num2str(j)])
        ylim([-1.5, 1.5]);
        xlabel('Tiempo [s]');
        ylabel('Amplitud');
    end
    sgtitle(['Modulación de los símbolos en el dominio del tiempo f(t) para \Deltaf = ' num2str(freq(k)) ])
end
    
%DOMINIO DE LA FRECUENCIA
for k=1:3
    figure(8+k)
    tiledlayout(4,2)
    for j=1:Nc
        nexttile
        sk = exp(1i*2*pi*(j-1)*freq(k)*t_total);
        mod= abs( fftshift(fft( bits(j)*sk,16000)));
        % Vector de frecuencia para todos los pulsos
        f_total = linspace(0,freq(k)*8,length(mod));
    
        plot(f_total,mod)
        title(['Símbolo = ' num2str(j)])
        xlabel('Frecuencia [Hz]');
    end
    sgtitle(['Modulación de los símbolos en el dominio de la frecuencia F(\omega) para \Deltaf = ' num2str(freq(k)) ])
end

%------- GRÁFICA DE s(t) y S(f)  para las frecuencias  ------------

figure(12)
tiledlayout(2,1)
nexttile
for k=1:3
    suma_st =sum(frecuencias{k});
    plot(t_total,abs(suma_st));
   
    hold on
end
title('Señal multiplexada s(t) en tiempo')
xlabel('Tiempo [s]');
legend('\Deltaf=1600','\Deltaf=2000','\Deltaf=3200')

nexttile
for k=1:3
    suma_st =sum(frecuencias{k});
    f_total = linspace(0,freq(k)*8,length(mod));
    plot(f_total,abs(fftshift(fft(suma_st,16000))));
    hold on
end
title('Señal multiplexada S(f) en frecuencia')
xlabel('Frecuencia [Hz]');
legend('\Deltaf=1600','\Deltaf=2000','\Deltaf=3200')

%--------------------  PARTE 2  ---------------------------------

n = [6,8,16,32,1600];
% figure (20) ifft(Y,n,2);

figure(13)
tiledlayout(5,1)
for j=1:5
    st_reconstruida = ifftshift(ifft(sum(frecuencias{1}),n(j)));
    nexttile
    f_total = linspace(0,8*delta_f,n(j));
    plot(f_total,abs(st_reconstruida));
    title( ['Señal reconstruida s(t) con N= ' num2str(n(j)) ] ) 

end
