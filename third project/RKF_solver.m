function [t, U] = RKF_solver(func, t_span, u0, para)
% RKF_solver implementing the Runge-Kutta-Fehlberg from "Numerical
% Analysis" of Burden & Faires with a Butcher tableau stored in a different
% function.
% The function accept a function, a time span, initial values and some
% parameters: ODE rank = the number of equations connecting the function and
% it derivitives, epsilon = the tolerance - function as relative tolerance
% and absolute one at once, h_init = the initial step size, h_max = largest
% step allowed, safety_factor = factor to slow down the step-size changing
    
    %unpacking the parameters struct
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
    n_vars = para.ode_rank;
    steps_limit = para.max_steps;

    %initialising the in-function parameters
    [A, err, b5, c] = RKF_Butcher_Tableau();
    solver_ks = length(b5);
    
    %in case the span isn't two variables as demanded, the t will begin in
    %zero
    if isscalar(t_span)
        t_span = [0, t_span];
    end
    
    %The pre-allocation of memory
    min_est    = max(ceil((t_span(end)-t_span(1))/h_max), 1000);
    t          = zeros(min_est, 1); %data point legend
    U          = zeros(n_vars, min_est); %variables stats
    K          = zeros(n_vars, solver_ks); %slopes

    %initial "condition"
    t_length = length(t);
    curr_step  = 1;
    U(:, 1) = u0(:);
    t(1) = t_span(1);
    n_steps = 0;
    

    while (t(curr_step) < t_span(end)) %there is a solution for the
        % t_span(end) data point
        if n_steps > steps_limit
            error("Too many iterations");
        end

        uk = U(:, curr_step);

        for j = 1:solver_ks %k1 to k6, for each variable
            if j == 1
                temp_prob = uk;
            else
                %muliplying the K's by the coefficients from A (which 
                % comes from the Bucther Tableau)
                temp_prob = uk + (h * (K(:, 1:j-1) * A(j, 1:j-1).'));
            end
            %evaluating the function with the current K at time t + ch when
            %c is the node from the bucther tableau
            K(:, j) = func(t(curr_step) + c(j)*h, temp_prob);
        end
        
        %using the K matrix (which now include the values of the function
        %and its derivitives at this point in time, t) and the weights for
        % 5th order evaluation to evaluate the function value at t+h.
        u5     = uk + h * (K * b5);
        
        du_err = h * (K * err); %same, but figuring the errors estimates.
        %R = (1/eps)*(1/h)*max(abs(du_err)); %Legacy error estimate (from 
        % "Numerical Analysis", was replaced by the new R, derived from
        % "Solving ODE I"
        scale_factor = eps * (1 + max(abs(uk), abs(u5))); %The epsilon is 
        % used as relative error and absolute one simoultansouly 
        R = max(abs(du_err) ./ scale_factor); %choosing the worst error
        % from the vector
        delta = safety_factor * (1 / (R)) ^0.25; %The reason for 0.25 is
        % that the gloabal  error is of size O(h^4) - this is the adaptive
        % part of the solver, it meant to optimize step-size; the safety
        % factor is there to dampen the change further.
        delta = max(0.1, min(4, delta)); %To prevent drastic changes to
        % step size
        
        if R <= 1 %step accepted, procced
            U(:, curr_step + 1) = u5;
            t(curr_step + 1) = t(curr_step) + h;            
            
            if curr_step > t_length %partial preallocation, 
                % 2000 elements each, in case the previous ones werent
                % enough
                t = [t; zeros(500, 1)];
                U = [U,  zeros(n_vars, 500)];
                t_length = t_length + 500;
            end
            curr_step = curr_step + 1;
        end
        
        %updating the step size and the step counter - to prevent endless
        %loops
        h       = min(delta*h, h_max);
        n_steps = n_steps + 1;
    end
    
    if t(curr_step) > t_span(2) %in case the last data point isn't t_span(2)
        curr_step = curr_step - 1; %returning to the privious data point
        h = t_span(2) - t(curr_step); %adjusting the step size to perfectly
        % align with t_span(2)
        uk = U(:, curr_step);
        
        %classic RK4, with no error handeling or adaptive step sizing.
        % it might be an overkill, but it will keep the global error at
        % O(h^4)
        k1 = func(t(curr_step),      uk); %dervs at the original point, u_k
        k2 = func(t(curr_step)+ h/2, uk + ((h/2) * k1)); %dervs at the point 
                %which will get to after a half step in f1 direction
        k3 = func(t(curr_step)+ h/2, uk + ((h/2) * k2)); %dervs at the point 
                %which will get to after a half step in f2 direction
        k4 = func(t(curr_step)+ h,   uk +  (h    * k3)); %dervs at the point 

        U(:, curr_step + 1) = uk + ((h/6) * (k1 + 2*k2 + 2*k3 + k4));
        %n_steps = n_steps + 1;
    end
    %trimming the redundant parts of the vectors
    t = t(1:curr_step + 1);
    U = U(:, 1:curr_step + 1);
end