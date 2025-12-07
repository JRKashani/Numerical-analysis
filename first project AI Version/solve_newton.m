function [root, iter] = solve_newton(func_handle, x0, tol, max_iter)
% Inputs:
%   func_handle: Handle to the function returning [f, df]
%   x0:          Initial guess
%   tol:         Tolerance for stop condition
%   max_iter:    Maximum allowed iterations

    x_current = x0;

    for iter = 1:max_iter
        [f, df] = func_handle(x_current);

        % Stop Condition 1: Function value is close enough to zero
        if abs(f) < tol
            root = x_current;
            return;
        end

        % Safety: Avoid division by zero if derivative vanishes
        if abs(df) < 1e-10 
            error('Derivative is too close to zero. Newton method failed.');
        end

        % Newton Step
        x_next = x_current - (f / df);

        % Stop Condition 2: Step size is smaller than tolerance (convergence)
        if abs(x_next - x_current) < tol
            root = x_next;
            return;
        end

        x_current = x_next;
    end

    error('Newton-Raphson method did not converge within max_iter.');
end