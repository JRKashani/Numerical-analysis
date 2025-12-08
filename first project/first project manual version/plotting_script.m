% Stable function (no singularities)
g = @(x, Bi) x .* besselj(1, x) - Bi .* besselj(0, x);

% Domain between first two roots of J0
x = 2.41:0.0001:5.505;
x1 = 2:0.001:6;

y1 = f(x, 0.1);
y2 = f(x, 60);
y3 = 500*besselj(0,x1);

figure('Color', 'w');
plot(x, y1, 'LineWidth', 1.6, 'DisplayName', 'Bi = 0.1');
hold on
plot(x, y2, 'LineWidth', 1.6, 'DisplayName', 'Bi = 60');
hold on
plot(x1, y3, 'LineWidth', 1.6, 'DisplayName', '500 * J_0(x)');

yline(0, 'k--', 'LineWidth', 1, 'DisplayName', 'y = 0'); % Show root line

xlabel('x');
ylabel('x * (J_1(x)/J_0(x)) - Bi');
%title('Stable Function Form (No Singularities)');
legend('Location','northwest');
grid on
hold off
