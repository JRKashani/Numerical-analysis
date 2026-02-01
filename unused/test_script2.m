%% Comparison Script: RK4 vs RK45 vs ode45 for Blasius Equation
% Tests performance, accuracy, and efficiency

clear; clc; close all;

%% Problem Setup: Blasius Boundary Layer Equation
% f''' + 0.5*f*f'' = 0
% f(0) = 0, f'(0) = 0, f''(0) = ? (shooting parameter)

% ODE system: [f, f', f''] -> [f', f'', -0.5*f*f'']
Blasius_derivs = @(t, y) [y(2); y(3); -0.5*y(1)*y(3)];

% Domain and initial conditions
eta_max = 10;
f_pp_0 = 0.33206;  % Known shooting value for Blasius
u0 = [0; 0; f_pp_0];

% Tolerances to test
tolerances = [1e-4, 1e-6, 1e-8];
n_tol = length(tolerances);

%% Storage for results
results = struct();

%% Run solvers for each tolerance
for idx = 1:n_tol
    tol = tolerances(idx);
    fprintf('\n=== Testing with tolerance = %.0e ===\n', tol);
    
    %% 1. Fixed-step RK4 (RKF_solver)
    fprintf('Running RK4 (fixed step)...\n');
    para_rk4.delta = 0.01;  % Fixed step size
    para_rk4.ode_rank = 3;
    
    tic;
    [t_rk4, U_rk4] = RKF_solver(Blasius_derivs, eta_max, u0, para_rk4);
    time_rk4 = toc;
    
    results(idx).rk4.t = t_rk4;
    results(idx).rk4.U = U_rk4;
    results(idx).rk4.time = time_rk4;
    results(idx).rk4.steps = length(t_rk4);
    results(idx).rk4.fevals = 4 * (length(t_rk4) - 1);  % 4 evals per step
    
    fprintf('  Time: %.4f s, Steps: %d\n', time_rk4, results(idx).rk4.steps);
    
    %% 2. Adaptive RK45 (RK45_solver)
    fprintf('Running RK45 (adaptive)...\n');
    para_rk45.epsilon = tol;
    para_rk45.h_max = 0.5;
    para_rk45.h_init = 0.01;
    para_rk45.ode_rank = 3;
    para_rk45.max_steps = 1e6;
    
    tic;
    [U_rk45, t_rk45] = RK45_solver(Blasius_derivs, eta_max, u0, para_rk45);
    time_rk45 = toc;
    
    results(idx).rk45.t = t_rk45;
    results(idx).rk45.U = U_rk45;
    results(idx).rk45.time = time_rk45;
    results(idx).rk45.steps = length(t_rk45);
    results(idx).rk45.fevals = 6 * results(idx).rk45.steps;  % Approx (includes rejections)
    
    fprintf('  Time: %.4f s, Steps: %d\n', time_rk45, results(idx).rk45.steps);
    
    %% 3. MATLAB's ode45
    fprintf('Running ode45...\n');
    options = odeset('RelTol', tol, 'AbsTol', tol);
    
    tic;
    [t_ode45, U_ode45] = ode45(Blasius_derivs, [0, eta_max], u0, options);
    time_ode45 = toc;
    
    results(idx).ode45.t = t_ode45;
    results(idx).ode45.U = U_ode45';  % Transpose to match column-major
    results(idx).ode45.time = time_ode45;
    results(idx).ode45.steps = length(t_ode45);
    
    fprintf('  Time: %.4f s, Steps: %d\n', time_ode45, results(idx).ode45.steps);
    
    %% 4. Compute reference solution (high accuracy ode45)
    if idx == n_tol  % Only once with tightest tolerance
        fprintf('Computing reference solution...\n');
        options_ref = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
        [t_ref, U_ref] = ode45(Blasius_derivs, [0, eta_max], u0, options_ref);
        
        % Store reference
        results(idx).ref.t = t_ref;
        results(idx).ref.U = U_ref';
    end
end

%% Use tightest tolerance result as reference
t_ref = results(end).ref.t;
U_ref = results(end).ref.U;

%% Compute errors (interpolate to common grid)
eta_common = linspace(0, eta_max, 1000);

for idx = 1:n_tol
    % Interpolate reference
    U_ref_interp = interp1(t_ref, U_ref', eta_common)';
    
    % RK4 errors
    U_rk4_interp = interp1(results(idx).rk4.t, results(idx).rk4.U', eta_common)';
    results(idx).rk4.error = max(abs(U_rk4_interp - U_ref_interp), [], 2);
    results(idx).rk4.error_max = max(results(idx).rk4.error);
    
    % RK45 errors
    U_rk45_interp = interp1(results(idx).rk45.t, results(idx).rk45.U', eta_common)';
    results(idx).rk45.error = max(abs(U_rk45_interp - U_ref_interp), [], 2);
    results(idx).rk45.error_max = max(results(idx).rk45.error);
    
    % ode45 errors
    U_ode45_interp = interp1(results(idx).ode45.t, results(idx).ode45.U', eta_common)';
    results(idx).ode45.error = max(abs(U_ode45_interp - U_ref_interp), [], 2);
    results(idx).ode45.error_max = max(results(idx).ode45.error);
end

%% Plotting

%% Figure 1: Solution Comparison (for middle tolerance)
idx_plot = 2;  % Use middle tolerance for plotting
figure('Position', [100, 100, 1200, 800]);

subplot(2,2,1)
plot(results(idx_plot).rk4.t, results(idx_plot).rk4.U(1,:), 'b-', 'LineWidth', 1.5); hold on;
plot(results(idx_plot).rk45.t, results(idx_plot).rk45.U(1,:), 'r--', 'LineWidth', 1.5);
plot(results(idx_plot).ode45.t, results(idx_plot).ode45.U(1,:), 'k:', 'LineWidth', 2);
xlabel('\eta'); ylabel('f(\eta)'); title('f(\eta)');
legend('RK4 (fixed)', 'RK45 (adaptive)', 'ode45', 'Location', 'best');
grid on;

subplot(2,2,2)
plot(results(idx_plot).rk4.t, results(idx_plot).rk4.U(2,:), 'b-', 'LineWidth', 1.5); hold on;
plot(results(idx_plot).rk45.t, results(idx_plot).rk45.U(2,:), 'r--', 'LineWidth', 1.5);
plot(results(idx_plot).ode45.t, results(idx_plot).ode45.U(2,:), 'k:', 'LineWidth', 2);
xlabel('\eta'); ylabel('f''(\eta)'); title('f''(\eta) - Velocity Profile');
legend('RK4 (fixed)', 'RK45 (adaptive)', 'ode45', 'Location', 'best');
grid on;
yline(1.0, '--', 'f''(\infty) = 1', 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);

subplot(2,2,3)
plot(results(idx_plot).rk4.t, results(idx_plot).rk4.U(3,:), 'b-', 'LineWidth', 1.5); hold on;
plot(results(idx_plot).rk45.t, results(idx_plot).rk45.U(3,:), 'r--', 'LineWidth', 1.5);
plot(results(idx_plot).ode45.t, results(idx_plot).ode45.U(3,:), 'k:', 'LineWidth', 2);
xlabel('\eta'); ylabel('f''''(\eta)'); title('f''''(\eta)');
legend('RK4 (fixed)', 'RK45 (adaptive)', 'ode45', 'Location', 'best');
grid on;

subplot(2,2,4)
% Step size visualization for adaptive methods
plot(results(idx_plot).rk45.t(1:end-1), diff(results(idx_plot).rk45.t), 'r-', 'LineWidth', 1.5); hold on;
plot(results(idx_plot).ode45.t(1:end-1), diff(results(idx_plot).ode45.t), 'k:', 'LineWidth', 1.5);
xlabel('\eta'); ylabel('Step size h'); title('Adaptive Step Size');
legend('RK45', 'ode45', 'Location', 'best');
grid on;

sgtitle(sprintf('Blasius Solution Comparison (tol = %.0e)', tolerances(idx_plot)));

%% Figure 2: Error Analysis
figure('Position', [150, 150, 1200, 400]);

subplot(1,3,1)
for idx = 1:n_tol
    semilogy(eta_common, results(idx).rk4.error(1,:), 'DisplayName', ...
        sprintf('RK4, tol=%.0e', tolerances(idx)), 'LineWidth', 1.5); hold on;
end
xlabel('\eta'); ylabel('Error in f(\eta)');
title('RK4 Error vs Reference');
legend('Location', 'best'); grid on;

subplot(1,3,2)
for idx = 1:n_tol
    semilogy(eta_common, results(idx).rk45.error(1,:), 'DisplayName', ...
        sprintf('RK45, tol=%.0e', tolerances(idx)), 'LineWidth', 1.5); hold on;
end
xlabel('\eta'); ylabel('Error in f(\eta)');
title('RK45 Error vs Reference');
legend('Location', 'best'); grid on;

subplot(1,3,3)
for idx = 1:n_tol
    semilogy(eta_common, results(idx).ode45.error(1,:), 'DisplayName', ...
        sprintf('ode45, tol=%.0e', tolerances(idx)), 'LineWidth', 1.5); hold on;
end
xlabel('\eta'); ylabel('Error in f(\eta)');
title('ode45 Error vs Reference');
legend('Location', 'best'); grid on;

%% Figure 3: Performance Comparison
figure('Position', [200, 200, 1200, 500]);

subplot(1,3,1)
% Runtime comparison
times = zeros(3, n_tol);
for idx = 1:n_tol
    times(1, idx) = results(idx).rk4.time;
    times(2, idx) = results(idx).rk45.time;
    times(3, idx) = results(idx).ode45.time;
end
bar(times');
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.0e', x), tolerances, 'UniformOutput', false));
xlabel('Tolerance'); ylabel('Runtime (seconds)');
title('Runtime Comparison');
legend('RK4', 'RK45', 'ode45', 'Location', 'best');
grid on;

subplot(1,3,2)
% Steps comparison
steps = zeros(3, n_tol);
for idx = 1:n_tol
    steps(1, idx) = results(idx).rk4.steps;
    steps(2, idx) = results(idx).rk45.steps;
    steps(3, idx) = results(idx).ode45.steps;
end
bar(steps');
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('%.0e', x), tolerances, 'UniformOutput', false));
xlabel('Tolerance'); ylabel('Number of Steps');
title('Number of Steps');
legend('RK4', 'RK45', 'ode45', 'Location', 'best');
grid on;

subplot(1,3,3)
% Efficiency: Error vs Time
loglog([results.rk4.time], [results.rk4.error_max], 'b-o', 'LineWidth', 2, 'MarkerSize', 8); hold on;
loglog([results.rk45.time], [results.rk45.error_max], 'r-s', 'LineWidth', 2, 'MarkerSize', 8);
loglog([results.ode45.time], [results.ode45.error_max], 'k-d', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Runtime (seconds)'); ylabel('Maximum Error');
title('Efficiency: Error vs Runtime');
legend('RK4', 'RK45', 'ode45', 'Location', 'best');
grid on;

%% Print Summary Table
fprintf('\n\n========== PERFORMANCE SUMMARY ==========\n\n');
fprintf('%-10s | %-8s | %-8s | %-10s | %-12s\n', ...
    'Tolerance', 'Solver', 'Time(s)', 'Steps', 'Max Error');
fprintf('%s\n', repmat('-', 1, 65));

for idx = 1:n_tol
    fprintf('%-10.0e | %-8s | %8.4f | %10d | %12.2e\n', ...
        tolerances(idx), 'RK4', results(idx).rk4.time, results(idx).rk4.steps, results(idx).rk4.error_max);
    fprintf('%-10s | %-8s | %8.4f | %10d | %12.2e\n', ...
        '', 'RK45', results(idx).rk45.time, results(idx).rk45.steps, results(idx).rk45.error_max);
    fprintf('%-10s | %-8s | %8.4f | %10d | %12.2e\n', ...
        '', 'ode45', results(idx).ode45.time, results(idx).ode45.steps, results(idx).ode45.error_max);
    fprintf('%s\n', repmat('-', 1, 65));
end

%% Key Metrics at boundary condition
fprintf('\n========== BOUNDARY CONDITION CHECK ==========\n');
fprintf('Target: f''(eta_max) should approach 1.0\n\n');
fprintf('%-10s | %-8s | f''(%.1f) = %.8f\n', 'Tolerance', 'Solver', eta_max, '');
fprintf('%s\n', repmat('-', 1, 50));

for idx = 1:n_tol
    fprintf('%-10.0e | %-8s | f''(%.1f) = %.8f\n', ...
        tolerances(idx), 'RK4', eta_max, results(idx).rk4.U(2, end));
    fprintf('%-10s | %-8s | f''(%.1f) = %.8f\n', ...
        '', 'RK45', eta_max, results(idx).rk45.U(2, end));
    fprintf('%-10s | %-8s | f''(%.1f) = %.8f\n', ...
        '', 'ode45', eta_max, results(idx).ode45.U(2, end));
    fprintf('%s\n', repmat('-', 1, 50));
end

fprintf('\nReference: f''(∞) ≈ 0.99999... (should be very close to 1.0)\n');