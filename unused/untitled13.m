clear all

u0 = [0; 0; 2*0.332057336];

h = 0.001;
eta_span = [0:h:50];

U(:, 1) = u0;
u_in = u0;

for i = 1:eta_span(end)/h
    eta = i*h;
    u_out = RK4_solver(@(eta, u)simp_Blasius_dervs(eta, u), h, eta, u_in);
    U = [U u_out];
    u_in = u_out;
end

%plot(eta_span, U(1, :))
% plot3(U(1, :), U(2, :), U(3, :), 'r', 'LineWidth', 3);
% hold on
% set(gca, 'FontSize', 24)
% view(20, 40)
% 
% [t,y] = ode45(@(eta, u)simp_Blasius(eta, u), eta_span, u0);
% y = y';
% plot3(y(1, :), y(2, :), y(3, :), 'c--', 'LineWidth', 3);
% set(gcf, 'Position', [1500 500 800 600])
% legend('My RK4', 'MATLAB ode45');
plot(eta_span, U(2, :), 'b', eta_span, U(3, :), 'g')