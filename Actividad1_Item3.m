% Actividad 1 - Ítem 3
close all; clear all; clc;
pkg load io
pkg load control
s = tf('s');

% Importar datos desde el archivo Excel
tabla = 'Curvas_Medidas_RLC_2025.xlsx';
data = xlsread(tabla, 1);
t = data(:, 1); % Tiempo
I_original = data(:, 2); % Corriente en el circuito
Vin_original = data(:, 4); % Tensión de excitación

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

% Definimos la función de excitación
u = zeros(size(t)); % Inicializar la entrada con ceros
indices_u = find(t > 0.01);
u(indices_u) = 12*(-1).^(floor((t(indices_u)/0.05)));

% Interpolar la señal de excitación u(t)
u_interp = interp1(t, u, t, 'previous');

% Esta ultima linea toma los puntos de "u" en los valores de "t" encontrados
% en la tabla, de manera que se puedan simular y graficar conjuntamente las
% respuestas de corriente y así comparar ambas curvas

% La función de transferencia del sistema con los parámetros calculados,
% tomando la corriente como salida, se desarrolla a continuación:
R = 220; % [Ω]
L = 4.155*10^(-3); % [Hy]
C = 2.2787*10^(-6); % [F]

%               V_in = VR + VL + VC
%             V_in = (R + sL + 1/(sC))*I
%         FT = I/V_in = 1/(sL + R + 1/(sC))

Sys_Model = 1/(s*L + R + 1/(s*C))
[I_model, t] = lsim(Sys_Model, u_interp, t);

% Se plotean los datos de la tabla junto con la corriente del modelo inferido
% para comparar la similitud entre ambos
figure;
plot(t, I_original, 'b', 'DisplayName', 'Respuesta Tabulada');
hold on;
plot(t, I_model, 'r', 'DisplayName', 'Respuesta del Modelo');
xlabel('Tiempo (s)');
ylabel('Corriente (A)');
title('Comparación de la Respuesta Tabulada y del Modelo');
xlim([0.05, max(t)]);
grid on;


disp("Terminado");
