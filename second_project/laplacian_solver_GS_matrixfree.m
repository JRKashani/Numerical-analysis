function p = laplacian_solver_GS_matrixfree(p)
% Solves the steady-state 2D heat equation using matrix-free Gauss–Seidel

    T = p.temp_mat;
    eps = p.epsilon;

    for iter = 1:10000
        T_old = T;

        % Gauss–Seidel update (interior only)
        T(2:end-1, 2:end-1) = 0.25 * ( ...
            T(1:end-2, 2:end-1) + ...
            T(3:end,   2:end-1) + ...
            T(2:end-1, 1:end-2) + ...
            T(2:end-1, 3:end) );

        % Convergence check (infinity norm)
        if norm(T - T_old, inf) < eps
            break
        end

        % Divergence safety check
        if norm(T, inf) > 1e4
            error("Solution diverging");
        end
    end

    if iter == 10000
        error("Maximum iterations exceeded");
    end

    p.temp_mat = T;
end
