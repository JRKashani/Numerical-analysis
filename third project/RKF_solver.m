function [t, U] = RKF_solver(func, t_span, u0, para)
% RKF_solver implementing the Runge-Kutta-Fehlberg from "Numerical
% Analysis" of Burden & Faires with a Butcher tableau stored in a different
% function.
% The function accept a function, a time span, initial values and some
% parameters: ODE rank = the number of equations connecting the function and
% it derivitives, epsilon = the tolerance - function as relative tolerance
% and absolute one at once, h_init = the initial step size, h_max = 
    
    if isfield(para, 'h_init')
        h = para.h_init;
    else
        h = 0.01;
    end
    if isfield(para, 'h_max')
        h_max = para.h_max;
    else
        h_max = 0.1;
    end
    if isfield(para, 'safety_factor')
        safety_factor = para.safety_factor;
    else
        safety_factor = 0.9;
    end
    
    eps = para.epsilon;
    [A, err, b5, c] = RKF_Butcher_Tableau();
    solver_ks = length(b5);
    n_vars = para.ode_rank;
    steps_limit = para.max_steps;
    n_steps = 0;

    if isscalar(t_span)
        t_span = [0, t_span];
    end
    
    min_est    = max(ceil((t_span(end)-t_span(1))/h_max), 1000);

    t          = zeros(10*min_est, 1); %data point legend
    U          = zeros(n_vars, 10*min_est); %variables stats
    K          = zeros(n_vars, solver_ks); %slopes

    t_length = length(t);
    curr_step  = 1;
    U(:, 1) = u0(:);
    t(1) = t_span(1);
    
    while (t(curr_step) < t_span(end))
        if n_steps > steps_limit
            error("Too many iterations");
        end
        uk = U(:, curr_step);
        for j = 1:solver_ks
            if j == 1
                temp_prob = uk;
            else
                temp_prob = uk + (h * (K(:, 1:j-1) * A(j, 1:j-1).'));
            end
            K(:, j) = func(t(curr_step) + c(j)*h, temp_prob);
        end

        u5     = uk + h * (K * b5);
        
        du_err = h * (K * err);
        %R     = (1/eps)*(1/h)*max(abs(du_err)); %Legacy error estimate (from 
        % "Numerical Analysis", was replaced by the new R, derived from
        % "Solving ODE I"
        scale_factor = eps * (1 + max(abs(uk), abs(u5)));
        R = max(abs(du_err) ./ scale_factor);
        delta = safety_factor * (1 / (R)) ^0.2; %The reason for 0.2 is that the
        % local error is of size O(h^5)        
        delta = max(0.1, min(4, delta));        
        
        if R <= 1 %step accepted, procced
            U(:, curr_step + 1) = u5;
            t(curr_step + 1) = t(curr_step) + h;            
            
            if curr_step > t_length %partial preallocation, 
                % 2000 elements each
                t = [t; zeros(500, 1)];
                U = [U,  zeros(n_vars, 500)];
                t_length = t_length + 500;
            end

            curr_step = curr_step + 1;
        end

        h       = min(delta*h, h_max);
        n_steps = n_steps + 1;
    end
    
    if t(curr_step) > t_span(2)
        curr_step = curr_step - 1;
        h = t_span(2) - t(curr_step);
        uk = U(:, curr_step);

        k1 = func(t(curr_step),      uk); %dervs at the original point, u_k
        k2 = func(t(curr_step)+ h/2, uk + ((h/2) * k1)); %dervs at the point 
                %which will get to after a half step in f1 direction
        k3 = func(t(curr_step)+ h/2, uk + ((h/2) * k2)); %dervs at the point 
                %which will get to after a half step in f2 direction
        k4 = func(t(curr_step)+ h,   uk +  (h    * k3)); %dervs at the point 

        U(:, curr_step + 1) = uk + ((h/6) * (k1 + 2*k2 + 2*k3 + k4));
        %n_steps = n_steps + 1;
    end
    t = t(1:curr_step + 1);
    U = U(:, 1:curr_step + 1);
end