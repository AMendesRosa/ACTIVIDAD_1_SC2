% Actividad 1 - Ítem 1
close all; clear all; clc
pkg load control
pkg load signal

% Parámetros del circuito
R = 220; % Resistencia [Ω]
L = 500e-3; % Inductancia [Hy]
C = 2.2e-6; % Capacitancia [F]

% Matrices del sistema en variables de estado
A = [-R/L, -1/L; 1/C, 0];
B = [1/L; 0];
C = [R, 0];
D = [0];
% La matriz C se da así dado que el vector de estados es [I VC]

% Ahora sacamos la funcion de transferencia según las matrices
[num, den] = ss2tf(A, B, C, D);
G = tf(num, den)

% Calculamos los polos de la FT y tomamos su parte real
polos = pole(G)
realPolos = real(polos);

% Definimos el polo más cercano y más lejano al origen
polo_cercano = min(realPolos)
polo_lejano = max(realPolos)

% Estimamos el tiempo de muestreo que permita observar las dinámicas rápidas
% según el polo más lejano del origen
tR = -log(0.95)/abs(polo_lejano)

% Y para estimar el tiempo de simulación utilizamos el polo más cercano
% al origen, y así permita abarcar toda la dinámica lenta
tL = -log(0.05)/abs(polo_cercano)

% Para ver varias repeticiones de la respuesta del sistema, tomamos
% un valor mayor de tiempo de simulación. Mientras que el tiempo de muestreo
% será al menos 10 veces menor al tR
t_muestreo = 1e-5;
t_sim = 0.2;

% Tiempo de simulación con 10 ms de entrada 0V
t = linspace(0, t_sim, t_sim/t_muestreo);
u = linspace(0, 0, t_sim/t_muestreo); % Inicializar la entrada con ceros

% Se define "u" a partir de t = 10ms con saltos entre 12 V y -12 V cada 10 ms
u(t > 0.01) = 12*(-1).^(floor((t(t > 0.01) - 0.01)/0.01));

% Los cambios de "u" se deben a que el exponente va haciendo un recuento de los
% periodos de 10 ms que ocurren y de esa manera el (-1) se vuelve una secuencia
% de "1" y "-1", ya que el floor convierte el exponente en valores enteros.

I(1) = 0;
VC(1) = 0;
VR(1) = 0; % Esta es la salida "y" del sistema
x = [I(1) VC(1)]'; % Planteamos las condiciones iniciales
x0 = [0 0]'; % Planteamos el punto de operación


% Le damos valores a la salida y las variables de estado en cada punto
for i = 1:(t_sim/t_muestreo) - 1;
  % Ecuaciones de estado y de salida
  x_punto = A*(x - x0) + B*u(i);
  x = x + x_punto*t_muestreo; % Integración de Euler
  y = C*x;

  % Salidas y variables de estado
  VR(i+1) = y;
  I(i+1) = x(1);
  VC(i+1) = x(2);
end

%Finalmente grafico la entrada, salida y variables de estado
figure;
subplot(4,1,1);
plot(t, u);
title('Tension de entrada');
xlabel('Tiempo (s)');
ylabel('Voltaje (V)');
grid on

subplot(4,1,2);
plot(t, VC);
title('Tension en el capacitor');
xlabel('Tiempo (s)');
ylabel('Voltaje (V)');
grid on

subplot(4,1,3);
plot(t, I);
title('Corriente en la malla');
xlabel('Tiempo (s)');
ylabel('Corriente (A)');
grid on;

subplot(4,1,4);
plot(t, VR);
title('Tension en la resistencia');
xlabel('Tiempo (s)');
ylabel('Voltaje (V)');
grid on;
