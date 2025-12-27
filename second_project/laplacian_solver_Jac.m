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
    % The iteration requires the inverse of D. Given it will be faster, I
    % avoided using the L and U matrises and preferred subtracting the main
    % diagonal from the main matrix.

    %L = tril(linear_equations_matrix, -1);
     D = spdiags(diag(linear_equations_matrix), 0, p.n_parameters,...
         p.n_parameters);
    %U = triu(linear_equations_matrix, 1);
    %LU = L+U; 
    LU = linear_equations_matrix - D; % Pre-calculated to avoid
                                      % re-computing inside the loop
    D1 = p.inv_D; %from the laplacian_mat function

    while flag == 0
        counter = counter + 1;
        % The Jacobi Iteration Step: x_new = D_inv * (b - (L+U)*x_old)
        x_new = D1 * (b - LU * x_old);

        % Stopping Condition:
        % If the maximum change between iterations is less
        % than epsilon, the solution has converged.
        % Given the run-time price, will be checked only once every 10
        % iterations
        if mod(counter, 10) == 0
            if counter > 100000
                error("exceeded the maximum number of iterations");
            end
            
            diff_vect = abs(x_new - x_old);
            if max(diff_vect) < p.epsilon
                flag = 1;
            elseif max(diff_vect) > 100
                error("starting to diverge, no point of continuing");
            end
        end
        
        x_old = x_new;
    end
    p.temp_mat = reshape(x_old, p.mat_size_1, p.mat_size_2);
end