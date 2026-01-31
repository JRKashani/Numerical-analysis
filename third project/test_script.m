%% Compare RK4, RKF, and ode45
clear; clc; close all;

%% Problem definition
func = @(t, u) [ ...
    u(2);        % x' = v
   -u(1)         % v' = -x
];

u0 = [1; 0];          % x(0)=1, v(0)=0
t_span = [0 20];

%% Parameters
para.ode_rank = 2;
para.delta    = 0.1;     % step size for RK4
para.tol      = 1e-6;     % tolerance for RKF (if used)
para.h_min    = 1e-5;
para.h_max    = 0.2;

%% Run RK4
[t_rk4, U_rk4] = RK4_solver(func, t_span, u0, para);

%% Run RKF (adaptive)
%[~, t_rkf, U_rkf] = RKF_solver(func, t_span, u0, para);

%% Run ode45
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[t_ode45, U_ode45] = ode45(func, t_span, u0, opts);
U_ode45 = U_ode45.';   % make it states × time

%% Exact solution
t_exact = linspace(t_span(1), t_span(2), 2000);
x_exact = cos(t_exact);

%% Plot position x(t)
figure;
plot(t_exact, x_exact, 'k--', 'LineWidth',1.5); hold on;
plot(t_rk4,   U_rk4(1,:), 'b');
%plot(t_rkf,   U_rkf(1,:), 'r.');
plot(t_ode45, U_ode45(1,:), 'g');

%legend('Exact','RK4','RKF','ode45','Location','best');
legend('Exact','RK4','ode45','Location','best');
xlabel('t'); ylabel('x(t)');
title('ODE Solver Comparison');
grid on;

%% Error comparison at final time
x_exact_final = cos(t_span(end));

err_rk4   = abs(U_rk4(1,end)   - x_exact_final);
%err_rkf   = abs(U_rkf(1,end)   - x_exact_final);
err_ode45 = abs(U_ode45(1,end) - x_exact_final);

fprintf('\nFinal time error:\n');
fprintf('RK4   : %.3e\n', err_rk4);
%fprintf('RKF   : %.3e\n', err_rkf);
fprintf('ode45 : %.3e\n', err_ode45);
