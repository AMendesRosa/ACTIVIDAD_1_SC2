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
VC_original = data(:, 3);
Vin_original = data(:, 4);

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

% Definimos paso y tiempo de simulación iguales a los de la tabla:
t_sim = 0.2;
paso = 0.00001;

% Definimos la función de excitación
u = linspace(0, 0, t_sim/paso); % Inicializar la entrada con ceros
u(t > 0.01) = 12*(-1).^(floor((t(t > 0.01)/0.05)));

% La función de transferencia del sistema con los parámetros calculados,
% tomando la corriente como salida, se desarrolla a continuación:
R = 220; % [Ω]
L = 4.155*10^(-3); % [Hy]
C = 2.2787*10^(-6); % [F]

%               V_in = VR + VL + VC
%             V_in = (R + sL + 1/(sC))*I
%         FT = I/V_in = 1/(sL + R + 1/(sC))

Sys_Model = 1/(s*L + R + 1/(s*C))
[I_model, t] = lsim(Sys_Model, u, t);

% Se plotean los datos de la tabla junto con la corriente del modelo inferido
% para comparar la similitud entre ambos
figure;
plot(t, I_original, 'b', 'DisplayName', 'Respuesta Tabulada');
hold on;
plot(t, I_model, 'r', 'DisplayName', 'Respuesta del Modelo');
xlabel('Tiempo (s)');
ylabel('Corriente (A)');
title('Comparación de la Respuesta Tabulada y del Modelo');
legend('Respuesta Tabulada', 'Respuesta del Modelo')
xlim([0.05, max(t)]);
grid on;

disp("Terminado");
