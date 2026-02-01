%% Compare RK4, RKF, and ode45
clear; clc; close all;

%% Problem definition (simple harmonic oscillator)
func = @(t, u) [ ...
    u(2);        % x' = v
   -u(1)         % v' = -x
];

u0 = [1; 0];          % x(0)=1, v(0)=0
t_span = [0 20];

%% Parameters
para.ode_rank  = 2;
para.delta     = 0.1;     % step size for RK4
para.epsilon   = 1e-6;    % tolerance for RKF
para.h_init    = 0.05;
para.h_max     = 0.2;
para.max_steps = 20000;
para.safety_factor = 0.9;

%% Run RK4 (fixed step)
[t_rk4, U_rk4] = RK4_solver(func, t_span, u0, para);

%% Run RKF (adaptive)
[U_rkf, t_rkf] = RKF_solver(func, t_span, u0, para);

%% Run ode45
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t_ode45, U_ode45] = ode45(func, t_span, u0, opts);
U_ode45 = U_ode45.';   % states × time

%% Exact solution
t_exact = linspace(t_span(1), t_span(2), 4000);
x_exact = cos(t_exact);
v_exact = -sin(t_exact);

%% ------------------------------------------------------------
%% 1. Time history comparison (x(t))
figure;
plot(t_exact, x_exact, 'k--', 'LineWidth',1.5); hold on;
plot(t_rk4,   U_rk4(1,:), 'b');
plot(t_rkf,   U_rkf(1,:), 'r.');
plot(t_ode45, U_ode45(1,:), 'g');

legend('Exact','RK4','RKF','ode45','Location','best');
xlabel('t'); ylabel('x(t)');
title('Position vs Time');
grid on;

%% ------------------------------------------------------------
%% 2. Phase space (x vs v)
figure;
plot(x_exact, v_exact, 'k--', 'LineWidth',1.5); hold on;
plot(U_rk4(1,:),   U_rk4(2,:),   'b');
plot(U_rkf(1,:),   U_rkf(2,:),   'r.');
plot(U_ode45(1,:), U_ode45(2,:), 'g');

legend('Exact','RK4','RKF','ode45','Location','best');
xlabel('x'); ylabel('v');
title('Phase Space');
axis equal; grid on;

%% ------------------------------------------------------------
%% 3. Step-size diagnostics
h_rkf   = diff(t_rkf);
h_ode45 = diff(t_ode45);

figure;
subplot(2,1,1)
plot(h_rkf, '.-');
ylabel('h');
title('RKF adaptive step sizes');
grid on;

subplot(2,1,2)
plot(h_ode45, '.-');
ylabel('h');
xlabel('Step index');
title('ode45 adaptive step sizes');
grid on;

%% ------------------------------------------------------------
%% 4. Accuracy vs cost
x_exact_final = cos(t_span(end));

err_rk4   = abs(U_rk4(1,end)   - x_exact_final);
err_rkf   = abs(U_rkf(1,end)   - x_exact_final);
err_ode45 = abs(U_ode45(1,end) - x_exact_final);

fprintf('\nFinal time error:\n');
fprintf('RK4   : %.3e   (%d steps)\n', err_rk4,   length(t_rk4)-1);
fprintf('RKF   : %.3e   (%d steps)\n', err_rkf,   length(t_rkf)-1);
fprintf('ode45 : %.3e   (%d steps)\n', err_ode45, length(t_ode45)-1);

fprintf('\nAverage step size:\n');
fprintf('RK4   : %.3e\n', mean(diff(t_rk4)));
fprintf('RKF   : %.3e\n', mean(diff(t_rkf)));
fprintf('ode45 : %.3e\n', mean(diff(t_ode45)));
