% Actividad 1 - Ítem 5
close all; clear all; clc;
pkg load control
pkg load io
pkg load signal
s = tf("s");

% Importar datos desde el archivo Excel
tabla = 'Curvas_Medidas_Motor_2025_v.xlsx';
data = xlsread(tabla, 1);
t = data(:, 1); % Tiempo
Wr_original = data(:, 2); % Velocidad angular
Ia_original = data(:, 3); % Corriente en la armadura
Vin_original = data(:, 4); % Tensión de entrada
TL_original = data(:, 5); % Torque de carga

% Ploteo los datos de la tabla para ver el comportamiento
figure;
subplot(4, 1, 1);
plot(t, Wr_original);
title('Velocidad Angular');
xlabel('Tiempo (s)');
ylabel('(Rad/s)');
grid on;

subplot(4, 1, 2);
plot(t, Ia_original);
title('Corriente en la armadura');
xlabel('Tiempo (s)');
ylabel('(A)');
grid on;

subplot(4, 1, 3);
plot(t, Vin_original);
title('Tensión de entrada');
xlabel('Tiempo (s)');
ylabel('Voltaje (V)');
grid on;

subplot(4, 1, 4);
plot(t, TL_original);
title('Torque de carga');
xlabel('Tiempo (s)');
ylabel('(N.m)');
grid on;

delay1 = 0.1; % Defino el retardo observado en la curva de Vin

% Extraigo los valores de t y Wr desde t = delay hasta t = 0.7
% (donde comienza el torque) para así eliminar el retardo del gráfico
indices1 = find(t >= 0.1 & t <= 0.7);
t_resp1 = t(indices1) - delay1;
Wr_resp1 = Wr_original(indices1);

% Ploteo la respuesta de Wr a Vin
figure;
plot(t_resp1, Wr_resp1);
title('Velocidad angular (respuesta al escalón de 12V)');
xlabel('Tiempo (s)');
ylabel('Vel. Angular (Rad/s)');
grid on;

% Sabiendo que es de segundo orden, por la forma de la curva tomaremos
% el caso que tiene polos reales y distintos, y para calcular los polos
% de la FT usaremos el método de Chen, con la forma:

%         FT = K*(T3*s + 1)/((T1*s + 1)*(T2*s + 1)).

% Para eso debemos fijar 3 puntos:
t_aux1 = 0.01; % Será nuestro t1

% Como ese valor está sin el retardo, se alineará con los valores de "t" al
% sumar el delay, y se busca el índice con mínimo valor de diferencia con "t",
% devolviendo el índice donde ocurre "t_aux", para poder conocer el valor de
% salida en la curva, como también definir los puntos necesarios para el método:
[~, lugar] = min(abs(t_aux1 + delay1 - t));
y_t1_w = Wr_original(lugar);
t1_w = t(lugar) - delay1

[~, lugar] = min(abs(2*t1_w+ delay1 - t));
y_t2_w = Wr_original(lugar);
t2_w = t(lugar) - delay1

[~, lugar] = min(abs(3*t1_w + delay1 - t));
y_t3_w = Wr_original(lugar);
t3_w = t(lugar) - delay1

% Luego de obtener los puntos y sus valores en verdadera magnitud, tomamos un
% último punto que donde la salida esté en estado de régimen, para así con este
% normalizar la altura de la entrada y calcular los valores de ki correctos,
% ya que el método de Chen aplica para un escalón unitario.
yss_w = 7.63;

% Se normaliza K dividiendolo por la altura del escalón, como así
% también a los valores de salidas
K_w = yss_w/2;
k1_w = (y_t1_w/2)/K_w - 1;
k2_w = (y_t2_w/2)/K_w - 1;
k3_w = (y_t3_w/2)/K_w - 1;

b_w = 4*(k1_w^3)*k3_w-3*(k1_w^2)*(k2_w^2)-4*(k2_w^3)+k3_w^2+6*k1_w*k2_w*k3_w
a1_w = (k1_w*k2_w + k3_w - sqrt(b_w))/(2*(k1_w^2 + k2_w))
a2_w = (k1_w*k2_w + k3_w + sqrt(b_w))/(2*(k1_w^2 + k2_w))
beta_w = (2*k1_w^3 + 3*k1_w*k2_w + k3_w - sqrt(b_w))/sqrt(b_w)

T1_w = real((-t1_w/log(a1_w)))
T2_w = real((-t1_w/log(a2_w)))
T3_w = (beta_w*(T1_w - T2_w) + T1_w)

% Por lo tanto, la FT nos queda de la siguiente manera
Sys_Model_Aux_w = K_w/((s*T1_w + 1)*(s*T2_w + 1))
Sys_Model_w = step(2*Sys_Model_Aux_w, t_resp1);

% Comparamos la respuesta inferida con los valores de la tabla
figure;
plot(t_resp1, Sys_Model_w, 'r-', t_resp1, Wr_resp1, 'b-');
title('Velocidad Angular respecto a Vin');
xlabel('Tiempo (s)');
ylabel('Vel. Angular (Rad/s)');
legend('Inferida por el método de Chen', 'Original');
grid on;

% Ahora aplicaremos el mismo método para el comportamiento de Ia respecto
% de Vin

delay1 = 0.101; % Defino el retardo observado en la curva de Ia

% Extraigo los valores de t y Wr desde t = delay hasta t = 0.7
% (donde comienza el torque) para así eliminar el retardo del gráfico
indices1 = find(t >= 0.101 & t <= 0.7);
t_resp1 = t(indices1) - delay1;
Ia_resp1 = Ia_original(indices1);

% Ploteo la respuesta de Ia a Vin
figure;
plot(t_resp1, Ia_resp1);
title('Velocidad angular (respuesta al escalón de 12V)');
xlabel('Tiempo (s)');
ylabel('Vel. Angular (Rad/s)');
grid on;

% Tomamos la misma forma de ft:
%         FT = K*(T3*s + 1)/((T1*s + 1)*(T2*s + 1)).

t_aux2 = 0.005; % Será nuestro t1

[~, lugar] = min(abs(t_aux2 + delay1 - t));
y_t1_2 = Ia_original(lugar);
t1_2 = t(lugar) - delay1

[~, lugar] = min(abs(2*t1_2 + delay1 - t));
y_t2_2 = Ia_original(lugar);
t2_2 = t(lugar) - delay1

[~, lugar] = min(abs(3*t1_2 + delay1 - t));
y_t3_2 = Ia_original(lugar);
t3_2 = t(lugar) - delay1

yss_2 = 0.04;

% Se normaliza K y los valores de salidas
K_2 = yss_2/2;
k1_2 = (y_t1_2/2)/K_2 - 1;
k2_2 = (y_t2_2/2)/K_2 - 1;
k3_2 = (y_t3_2/2)/K_2 - 1;

b_2 = 4*(k1_2^3)*k3_2-3*(k1_2^2)*(k2_2^2)-4*(k2_2^3)+k3_2^2+6*k1_2*k2_2*k3_2
a1_2 = (k1_2*k2_2 + k3_2 - sqrt(b_2))/(2*(k1_2^2 + k2_2))
a2_2 = (k1_2*k2_2 + k3_2 + sqrt(b_2))/(2*(k1_2^2 + k2_2))
beta_2 = (k1_2+a1_2)/(a1_2-a2_2)

T1_2 = real((-t1_2/log(a1_2)))
T2_2 = real((-t1_2/log(a2_2)))
T3_2 = (beta_2*(T1_2 - T2_2) + T1_2)

% Por lo tanto, la FT nos queda de la siguiente manera
Sys_Model_Aux_2 = K_2*(s*T3_2 + 1)/((s*T1_2 + 1)*(s*T2_2 + 1))
Sys_Model_2 = step(2*Sys_Model_Aux_2, t_resp1);

% Comparamos la respuesta inferida con los valores de la tabla
figure;
plot(t_resp1, Sys_Model_2, 'r-', t_resp1, Ia_resp1, 'b-');
title('Corriente de armadura respecto a Vin');
xlabel('Tiempo (s)');
ylabel('Corriente (A)');
legend('Inferida por el método de Chen', 'Original');
grid on;

% Establecemos los parámetros calculados según lo explicado en el informe
Ra = 2.4136
Ki = 3.815
J = 0.04009
Bm = 0.02
Km = 0.24946
Laa = 5.2382e-3

% Armo las matrices del espacio de estados según las ecuaciones de estado
% y la ecuación de salida:
A = [-Ra/Laa, -Km/Laa, 0; Ki/J, -Bm/J, 0; 0, 1, 0];
B = [1/Laa, 0; 0, -1/J; 0, 0];
C = [0, 1, 0];
D = [0, 0];

% Calculo los autovalores de A para evaluar estabilidad (polos de FT), y
% tomamos su módulo para cálculos posteriores (omitiendo el 0):
polos = eig(A)
mod_polos = abs(real(polos));
mod_polos_sin_0 = mod_polos(mod_polos ~= 0);

% Definimos el polo más cercano y más lejano al origen
polo_cercano = min(mod_polos_sin_0)
polo_lejano = max(mod_polos_sin_0)

% Tiempo asociado al polo más rápido (mayor abs de parte real)
tR = -log(0.95)/polo_lejano

% Tiempo asociado al polo más lento (menor abs de parte real)
tL = -log(0.05)/polo_cercano

% Definimos el tiempo de simulación y de muestreo
t_sim = 1.5;
t_muestreo = 1e-5;
num_puntos = round(t_sim/t_muestreo) + 1;

% Inicializamos valores
t2 = linspace(0, t_sim, num_puntos);
u = zeros(1, num_puntos);
TL = zeros(1, num_puntos);
Ia = zeros(1, num_puntos);
Wr = zeros(1, num_puntos);
Tita = zeros(1, num_puntos);

% Definimos el vector de estados
x = [0 0 0]'; % Condiciones iniciales
x0 = [0 0 0]'; % Punto de operación

% Le damos valores a la salida y las variables de estado en cada punto
for i = 1:(t_sim/t_muestreo) - 1;
  if(t2(i) > 0.1)
    u(i) = 2;
  end

  if (t2(i) > 0.7)
    TL(i) = 0.12  ;
  end
  if (t2(i) > 1)
    TL(i) = 0;

    end
  % Ecuaciones de estado y de salida
  x_punto = A*(x - x0) + B*[u(i) TL(i)]';
  x = x + x_punto*t_muestreo; % Integración de Euler
  y = C*x;

  % Salidas y variables de estado
  Wr(i+1) = y;
  Ia(i+1) = x(1);
  Tita(i+1) = x(3);
end

%Finalmente grafico la entrada, salida y variables de estado
figure;
subplot(4,1,1);
plot(t2, Wr, "linewidth", 1.2);
title('Velocidad angular');
xlabel('Tiempo (s)');
ylabel('Vel. angular (rad/s)');
grid on

subplot(4,1,2);
plot(t2, Ia, "linewidth", 1.2);
title('Corriente en la armadura');
xlabel('Tiempo (s)');
ylabel('Corriente (A)');
grid on;

subplot(4,1,3);
plot(t2, u, "linewidth", 1.2);
title('Tensión de entrada');
xlabel('Tiempo (s)');
ylabel("Voltaje (V)");
grid on;

subplot(4,1,4);
plot(t2, TL, "linewidth", 1.2);
title('Torque de carga');
ylabel("Torque (N.m)");
xlabel('Tiempo (s)');
grid on;

