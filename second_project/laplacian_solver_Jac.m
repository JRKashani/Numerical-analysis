function [p] = laplacian_solver_Jac(p)
% LAPLACIAN_SOLVER_JAC receives a structure of parameters (including the
% initial temperature guess and the Laplacian matrix) and returns the same
% structure with the stabilized (converged) heat matrix.

    %unpacking the struct
    linear_equations_matrix = p.laplacian_mat;
    x_old = p.temp_mat(:);

    % Construct the RHS vector 'b'.
    % The interior equations equal 0, but boundary rows equal the fixed
    % temperatures. We use element-wise multiplication with the boundary
    % mask to ensure only boundary values are non-zero.
    b = x_old .* p.the_boundary(:);

    counter = 0;
    flag = 0;
    
    % For the Jacobi method, the matrix A is decomposed into three parts:
    % Lower triangular (L), Upper triangular (U), and Diagonal (D).
    % The iteration requires the inverse of D.
    L = tril(linear_equations_matrix, -1);
    D = spdiags(diag(linear_equations_matrix), 0, p.n_parameters,...
        p.n_parameters);
    U = triu(linear_equations_matrix, 1);
    LU = L+U; % Pre-calculated to avoid re-computing inside the loop
    D1 = inverse_this_diag_matrix(D);

    while flag == 0
        counter = counter + 1;
        if counter > 10000
            error("exceeded the maximum number of iterations");
        end

        % The Jacobi Iteration Step: x_new = D_inv * (b - (L+U)*x_old)
        x_new = D1 * (b - LU * x_old);

        % Stopping Condition:
        % If the maximum change between iterations is less
        % than epsilon, the solution has converged.
        if max(abs(x_new - x_old)) < p.epsilon
            flag = 1;
        end
        
        x_old = x_new;
    end
    p.temp_mat = reshape(x_old, p.mat_size_1, p.mat_size_2);
end