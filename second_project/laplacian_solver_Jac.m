function [p] = laplacian_solver_Jac(p)

    linear_equations_matrix = p.laplacian_mat;
    x_old = p.temp_mat(:);
    b = x_old .* p.the_boundary(:);
    counter = 0;
    flag = 0;
    
    L = tril(linear_equations_matrix, -1);
    D = spdiags(diag(linear_equations_matrix), 0, p.n_parameters,...
        p.n_parameters);
    U = triu(linear_equations_matrix, 1);
    LU = L+U;
    D1 = inverse_this_diag_matrix(D);

    while flag == 0
        counter = counter + 1;
        if counter > 10000
            error("The loop had too many iterations");
        end

        x_new = D1 * (b - LU * x_old);

        if max(abs(x_new - x_old)) < p.epsilon
            flag = 1;
        end

        x_old = x_new;
    end
    p.temp_mat = reshape(x_old, p.mat_size_1, p.mat_size_2);
end