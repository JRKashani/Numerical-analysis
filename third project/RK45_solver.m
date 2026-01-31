function [U, t, ] = RK45_solver(func, t_span, u0, para)
    
    eps = para.epsilon;

    [A, err, b5, c] = RKF_Butcher_Tableau();
    solver_ks = 6;
    U = zeros(length(t), para.ode_rank);

    rel_slopes = zeros(para.ode_rank, solver_ks);
    curr_step = 0;
    while (t(curr_step) < t_span)
        
        for j = 1:solver_ks
            if j == 1
                temp_prob = uk;
            else
                temp_prob = uk+(h* (rel_slopes(:, 1:j-1) * A(j, 1:j-1).'));
            end
            rel_slopes(:, j) = func(t(curr_step) + c(j)*h, temp_prob);
        end
        u5 = uk + h * (rel_slopes * b5);
        u_err = uk + h * (rel_slopes * err);
        R = err_coe*(1/h)*max(abs(u_err));
        if R <= eps
            U(curr_step + 1, :) = u5;
            curr_step = curr_step + 1;
        else
            delta = ()
        end
        
        
        rel_slopes(:,:) = 0;
    end

end