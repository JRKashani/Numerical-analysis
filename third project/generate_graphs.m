% --- Setup & Parameters ---
clear; clc; close all;

% The correct shooting parameter found (approx 0.33206)
f_double_tag_0 = 0.332057336215; 

% Generate the data using your wrapper
[eta, U] = Blasius_Wrapping(f_double_tag_0);

% Extract variables for clarity
f       = U(1, :);
f_tag   = U(2, :); % Velocity profile (u/U_inf)
f_tag_tag = U(3, :); % Shear stress

% --- Figure 1: The Blasius Profiles ---
figure('Name', 'Blasius Profiles', 'Color', 'w');
hold on;
grid on;

% Plotting with distinct line styles
plot(eta, f,       '-k', 'LineWidth', 1.5, 'DisplayName', '$f(\eta)$');
plot(eta, f_tag,   '--b', 'LineWidth', 1.5, 'DisplayName', '$f''(\eta)$ (Velocity)');
plot(eta, f_tag_tag,':r', 'LineWidth', 1.5, 'DisplayName', '$f''''(\eta)$ (Shear)');

% Formatting for LaTeX look
xlabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Function Value', 'Interpreter', 'latex', 'FontSize', 14);
title('Solution of the Blasius Equation', 'Interpreter', 'latex', 'FontSize', 16);
legend('show', 'Interpreter', 'latex', 'Location', 'NorthWest', 'FontSize', 12);
ylim([0, 5]); % Focus on the relevant range
xlim([0, 8]);

% Save the figure
exportgraphics(gcf, 'blasius_profiles.png', 'Resolution', 300);

% --- Figure 2: Adaptive Step Size Analysis ---
% Calculate the step size vector h
h_steps = diff(eta);
eta_mid = eta(1:end-1); % Corresponds to the step locations

figure('Name', 'Adaptive Step Size', 'Color', 'w');
plot(eta_mid, h_steps, '.-m', 'LineWidth', 1, 'MarkerSize', 8);
grid on;

% Formatting
xlabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Step Size ($h$)', 'Interpreter', 'latex', 'FontSize', 14);
title('RKF45 Adaptive Step Size', 'Interpreter', 'latex', 'FontSize', 16);

% Annotation to explain the graph
dim = [0.15 0.6 0.3 0.3];


% Save the figure
exportgraphics(gcf, 'step_size_analysis.png', 'Resolution', 300);

fprintf('Graphs generated and saved as PNG files.\n');