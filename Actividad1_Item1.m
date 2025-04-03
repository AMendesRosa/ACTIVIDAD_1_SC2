% Actividad 1 - Ítem 1
close all; clear all; clc
pkg load control

% Parámetros del circuito
R = 220; % Resistencia [Ω]
L = 500e-3; % Inductancia [Hy]
C = 2.2e-6; % Capacitancia [F]

% Matrices del sistema en variables de estado
A = [-R/L, -1/L; 1/C, 0];
B = [1/L; 0];
C_il = [1, 0]; % Salida: IL (primera variable de estado)
C_vc = [0, 1]; % Salida: VC (segunda variable de estado)
C_vr = [R, 0]; % Salida: VR (IL*R)
D = [0];

% La matriz C se da así dado que el vector de estados es [IL VC]

% Tiempo de simulación con 10 ms de entrada 0V
t = 0:1e-4:0.2;
u = zeros(size(t)); % Inicializar la entrada con ceros
u(t > 0.01) = 12*(-1).^(floor((t(t > 0.01) - 0.01)/0.01));

% Los cambios de signo de "u" son gracias al factor (-1) y su exponente, éste
% hace un recuento de los periodos de 10 ms que van ocurriendo y de esa manera
% el factor (-1), debido a que floor convierte el exponente en valores enteros,
% se transforma en una secuencia de "1" y "-1" que le dan el signo a "u".

% Almacenamos valores de cada salida en función del tiempo y la excitación
[y_vc, t_out] = lsim(ss(A, B, C_vc, D), u, t); % VC
[y_il, ~] = lsim(ss(A, B, C_il, D), u, t); % IL
[y_vr, ~] = lsim(ss(A, B, C_vr, D), u, t); % VR

% La función "ss" crea un objeto de espacio de estados basado en las
% matrices de su argumento.
% Se utiliza "~" ya que los mismos valores se guardan en "t_out"

% Graficar resultados
figure;

subplot(4, 1, 1);
plot(t_out, u); % Grafica la entrada
title('Entrada: Tensión aplicada');
xlabel('Tiempo (s)');
ylabel('Voltaje (V)');
grid on;

subplot(4, 1, 2);
plot(t_out, y_vc); % Grafica el voltaje sobre el capacitor
title('Voltaje en el capacitor (VC)');
xlabel('Tiempo (s)');
ylabel('Voltaje (V)');
grid on;

subplot(4, 1, 3);
plot(t_out, y_il); % Grafica la corriente sobre el inductor
title('Corriente en el inductor (i_L)');
xlabel('Tiempo (s)');
ylabel('Corriente (A)');
grid on;

subplot(4, 1, 4);
plot(t_out, y_vr); % Grafica la tensión sobre la resistencia
title('Caída de tensión en la resistencia (v_R)');
xlabel('Tiempo (s)');
ylabel('Voltaje (V)');
grid on;
