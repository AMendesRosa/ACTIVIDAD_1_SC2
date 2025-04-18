% Actividad 1 - Ítem 3
close all; clear all; clc;
pkg load signal
pkg load control

% Importar datos desde el archivo Excel
tabla = 'Curvas_Medidas_RLC_2025.xlsx';
data = xlsread(tabla, 1);
t = data(501:end, 1); % Tiempo
I_original = data(501:end, 2); % Corriente en el circuito
Vin_original = data(501:end, 4);

% Ploteo los datos originales de las tablas
subplot(2, 1, 1);
plot(t, Vin_original); % Grafica la tensión de entrada
title('Entrada: Tensión aplicada');
xlabel('Tiempo (s)');
ylabel('Voltaje (V)');
grid on;

subplot(2, 1, 2);
plot(t, I_original); % Grafica la corriente en la malla
title('Corriente en la malla');
xlabel('Tiempo (s)');
ylabel('Corriente (A)');
grid on;

% Definimos tiempo de muestreo igual al paso de la tabla:
t_muestreo = 0.00001;

% La función de transferencia del sistema con los parámetros calculados,
% se desarrolla a continuación:
R = 220; % [Ω]
L = 1.9591e-3; % [Hy]
C = 2.2306e-6; % [F]

% Definimos la FT según las matrices del sistema
A = [-R/L, -1/L; 1/C, 0];
B = [1/L; 0];
C = [R, 0];
D = [0];

[num, den] = ss2tf(A, B, C, D);
G = tf(num, den)

x0 = [0 0]'; % Planteamos el punto de operación

I(1) = 0;
VC(1) = 0;
x = [I(1) VC(1)]'; % Planteamos las condiciones iniciales

% Le damos valores a las variables de estado en cada punto
for i = 1:19500;
  % Ecuación
  x_punto = A*(x - x0) + B*Vin_original(i);
  x = x + x_punto*t_muestreo; % Integración de Euler

  % Variable de estado de interés
  I(i+1) = x(1);
end

% Se plotean los datos de la tabla junto con la corriente del modelo inferido
% para comparar la similitud entre ambos
figure;
plot(t, I_original, 'b-', t, I, "r-");
xlabel('Tiempo (s)');
ylabel('Corriente (A)');
title('Comparación de las respuestas de la corriente en la malla');
legend('Respuesta Tabulada', 'Respuesta del Modelo')
grid on;

disp("Terminado");
