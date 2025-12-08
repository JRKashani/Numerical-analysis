% --- Helper: Constrained Newton-Raphson ---
function [x, iter] = run_newton(func, x, tol, max_iter, xmin, xmax)
    for iter = 1:max_iter
        [f, df] = func(x);
        
        % Success Check
        if abs(f) < tol, return; end
        
        % Stability Check: Relative derivative magnitude
        if abs(df) < 1e-12 * (1 + abs(f))
            error('Derivative too small relative to f(x)');
        end
        
        x_new = x - f/df;
        
        % Constraint Check: Must stay inside bracket
        if x_new < xmin || x_new > xmax
            error('Newton iterate left the valid bracket');
        end
        if x_new > xmax - (xmax-xmin)/10
            x_new = (x_new + xmax)/2;
        end

        
        % Convergence Check
        if abs(x_new - x) < tol
            x = x_new;
            return;
        end
        x = x_new;
    end
    error('Newton did not converge within max_iter');
end