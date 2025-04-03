% Actividad 1 - Ítem 2
close all; clear all; clc;
pkg load io
pkg load control
s = tf('s');

% Importar datos desde el archivo Excel
tabla = 'Curvas_Medidas_RLC_2025.xlsx';
data = xlsread(tabla, 1);
t = data(:, 1); % Tiempo
I = data(:, 2); % Corriente en el circuito
VC_original = data(:, 3); % Voltaje en el capacitor

% Ploteo los datos de la tabla (sólo para VC) para ver el retardo
figure;
plot(t, VC_original);
title('Voltaje en el capacitor');
xlabel('Tiempo (s)');
ylabel('Voltaje (V)');
grid on;

delay = 0.01;

% Extraigo los valores de "t" y "VC" para t = 0.01 hasta t = 0.025 (donde
% ya se estabilizó la respuesta) para así eliminar el retardo del gráfico
indices = find(t >= 0.01 & t <= 0.025);
t_resp = t(indices) - delay;
VC_resp = VC_original(indices);

% Gráfico de la respuesta de VC al escalón de 12V sin el retardo
figure;
plot(t_resp, VC_resp);
title('Voltaje en el capacitor (respuesta al escalón de 12V)');
xlabel('Tiempo (s)');
ylabel('Voltaje (V)');
grid on;

% La función de transferencia del sistema VC/V_in es:
%             V_in = (R + sL + 1/(sC))*I
%           VC = 1/(sC)*V_in/(R+sL+1/(sC))
%         FT = VC/V_in = 1/(s^(2)*LC + sRC + 1)

% Por la forma de la curva tomaremos el caso que tiene polos reales y
% distintos, y para buscar la forma de la FT usaremos el método de Chen,
% asemejándola a la forma FT = K*(T3*s + 1)/((T1*s + 1)*(T2*s + 1)).

% Para eso debemos fijar 3 puntos:
t_aux = 5.5e-4; % Será nuestro t1

% Se alineará "t_aux" con los valores de "t" originales sumando "delay", y se
% busca el índice con mínimo valor de diferencia, devolviendo aquel donde
% ocurre "t_aux", y usarlo para así poder tanto conocer el valor de salida
% en la curva, como también obtener los puntos necesarios para el método:
[~, lugar] = min(abs(t_aux + delay - t));
y_t1 = VC_original(lugar);
t1 = t(lugar) - delay;

[~, lugar] = min(abs(2*t1 + delay - t));
y_t2 = VC_original(lugar);
t2 = t(lugar) - delay;

[~, lugar] = min(abs(3*t1 + delay - t));
y_t3 = VC_original(lugar);
t3 = t(lugar) - delay;

StepAmplitude = 12;
K = 1;
k1 = (1/StepAmplitude)*y_t1/K - 1;
k2 = (1/StepAmplitude)*y_t2/K - 1;
k3 = (1/StepAmplitude)*y_t3/K - 1;

b = 4*k1^3*k3 - 3*k1^2*k2^2 - 4*k2^3 + k3^2 + 6*k1*k2*k3
a1 = (k1*k2 + k3 - sqrt(b))/(2*(k1^2 + k2))
a2 = (k1*k2 + k3 + sqrt(b))/(2*(k1^2 + k2))
beta = (2*k1^3 + 3*k1*k2 + k3 - sqrt(b))/(sqrt(b))

T1 = (-t1/log(a1))
T2 = (-t1/log(a2))
T3 = (beta*(T1 - T2) + T1)
T1 = real(T1) % Se tomaría la parte real ya que vimos que la respuesta
T2 = real(T2) % no es subamortiguada, sin embargo los polos según los
T3 = real(T3) % puntos elegidos ya son reales.

% Despreciando el cero, la FT inferida nos queda de la siguiente manera:
Sys_Model_Aux = K/((s*T1 + 1)*(s*T2 + 1))
Sys_Model = step(12*Sys_Model_Aux, t_resp);

% Comparamos la respuesta inferida con los valores de la tabla
plot(t_resp, Sys_Model, 'r-', t_resp, VC_resp, 'b-');
title('Voltaje en el capacitor');
xlabel('Tiempo (s)');
ylabel('Voltaje (V)');
legend('Inferida por el método de Chen', 'Original');
grid on;

% Ahora se supondrá un valor de R = 220Ω, y se calculará L y C en base
% a la función de transferencia que se encontró por el método de Chen
pkg load symbolic
syms s_s T1_s T2_s C_s R_s L_s real

Sys_Model_Symb = K/((s_s*T1_s + 1)*(s_s*T2_s + 1))
RLC_Model_Symb = 1/(s_s^(2)*L_s*C_s + s_s*R_s*C_s + 1)

% Igualamos denominadores y asignamos valores a las variables para calcular
eq1 = L_s*C_s == T1_s*T2_s;
eq2 = T1_s + T2_s == C_s*R_s;
R_s = 220;
T1_s = T1;
T2_s = T2;

C_s = double(solve(eval(eq2), C_s));
L = double(solve(eval(eq1), L_s))
R = R_s
C = C_s

% Armamos el modelo de FT con los valores de R, L y C calculados
RLC_Model = step(StepAmplitude*1/(s^(2)*L*C + s*R*C + 1), t_resp);

% Graficamos la respuesta a un escalon de 12V en los 3 modelos para comparar
figure;
plot(t_resp, Sys_Model, 'r-', t_resp, RLC_Model, 'g-', t_resp, VC_resp, 'b-');
title('Voltaje en el capacitor');
xlabel('Tiempo (s)');
ylabel('Voltaje (V)');
legend('Inferida por método de Chen', 'Modelo con RLC', 'Original');
grid on;

disp('Terminado');

