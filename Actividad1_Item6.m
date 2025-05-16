% Actividad 1 - Ítem 6
clear all; close all; clc

t_sim = 35;
t_muestreo = 2e-4;
t = 0:t_muestreo:t_sim;

% Inicialización de vectores
X = [0; 0; 0; 0]; 
u = 2;
TL = 0;

% Parámetros del motor
Ra = 2.4136;
Ki = 0.274521;
J  = 2.8221615e-3;
Km = 0.262123;  
Laa = 5.283097e-3;

N = length(t);
x1 = zeros(N,1);
x2 = zeros(N,1);
x3 = zeros(N,1);
x4 = zeros(N,1);

function [X] = modmotor_fast(dt, xant, u, TL, Ra, Ki, J, Km, Laa);
  
    h = 1e-5;
    N = round(dt / h);

    omega = xant(1);
    ia    = xant(3);
    tita  = xant(4);

    for i = 1:N
        domega = (Ki * ia - TL) / J;
        dia    = (u - Ra * ia - Km * omega) / Laa;

        omega = omega + h * domega;
        ia    = ia    + h * dia;
        tita  = tita  + h * omega;
    end

    X = [omega, domega, ia, tita];
end

% Referencia de velocidad
theta_ref = ones(N, 1);  % 1 rad constante

% PID parámetros (a tunear)
Kp = 1;
Ki = 1;
Kd = 0.00001;

e  = zeros(N,1);
u_pid = zeros(N,1);


for i = 1:N
    tiempo = t(i);

    if i == 1
    e_k = theta_ref(i) - X(4);  % Error de ángulo
    e_k_1 = 0;
    e_k_2 = 0;
    u_pid(i) = 0;
elseif i == 2
    e_k = theta_ref(i) - X(4);
    e_k_1 = e(i-1);
    e_k_2 = 0;
    u_pid(i) = u_pid(i-1) + Kp*(e_k - e_k_1) + Ki*t_muestreo*e_k + (Kd/t_muestreo)*(e_k - 2*e_k_1);
else
    e_k = theta_ref(i) - X(4);
    e_k_1 = e(i-1);
    e_k_2 = e(i-2);
    u_pid(i) = u_pid(i-1) + Kp*(e_k - e_k_1) + Ki*t_muestreo*e_k + (Kd/t_muestreo)*(e_k - 2*e_k_1 + e_k_2);
end
e(i) = e_k;


    % Aplica TL perturbador a partir del valor indicado de tiempo
    if ((tiempo >= 9.5) && (tiempo <= 2.15))
        TL = 0.12;
    else
        TL = 0;
    end

    X = modmotor_fast(t_muestreo, X, u_pid(i), TL, Ra, Ki, J, Km, Laa);
    x1(i) = X(1);
    x2(i) = X(2);
    x3(i) = X(3);
    x4(i) = X(4);
end

% Gráficos
figure;
subplot(4,1,1); plot(t, x1, 'g'); title('\omega'); grid on;
subplot(4,1,2); plot(t, x4, 'r'); title('\theta'); grid on;
subplot(4,1,3); plot(t, x3, 'b'); title('ia'); grid on;
subplot(5,1,5); plot(t, u_pid, 'k'); title('u_{PID}'); grid on;
xlabel('Tiempo [s]');

