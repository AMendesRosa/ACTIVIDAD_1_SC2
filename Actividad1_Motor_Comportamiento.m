% Actividad 1 - Ítem 4
close all; clear all; clc;

% Definición de parámetros
Laa = 366e-6;
J = 5e-9;
Ra = 55.6;
Bm = 0;
Ki = 6.49e-3;
Km = 6.53e-3;

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
t_sim = 0.2;
t_muestreo = 1e-7;

% Tiempo de simulación con 2 ms de entrada 0V y 10 ms de TL = 0
t = linspace(0, t_sim, t_sim/t_muestreo);
u = linspace(0, 0, t_sim/t_muestreo); % Inicializar la entrada con ceros
TL = linspace(0, 0, t_sim/t_muestreo); % Inicializar TL con ceros

% Definimos el vector de estados
Ia(1) = 0;
Wr(1) = 0;
Tita(1) = 0;
x = [Ia(1) Wr(1) Tita(1)]'; % Condiciones iniciales
x0 = [0 0 0]'; % Punto de operación

% Le damos valores a la salida y las variables de estado en cada punto
for i = 1:(t_sim/t_muestreo) - 1;
  if(t(i) > 0.01)
    u(i) = 12;
  end

  if (t(i) > 0.06)
    TL(i) = 5e-4;
  end
  if (t(i) > 0.12)
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
plot(t, u, "linewidth", 1.2);
title('Tension de entrada');
xlabel('Tiempo (s)');
ylabel('Voltaje (V)');
grid on

subplot(4,1,2);
plot(t, Wr, "linewidth", 1.2);
title('Velocidad angular');
xlabel('Tiempo (s)');
ylabel('Radianes por segundo (rad/s)');
grid on

subplot(4,1,3);
plot(t, Ia, "linewidth", 1.2);
title('Corriente en la armadura');
xlabel('Tiempo (s)');
ylabel('Corriente (A)');
grid on;

subplot(4,1,4);
plot(t, TL, "linewidth", 1.2);
title('Torque de carga');
xlabel('Tiempo (s)');
grid on;







